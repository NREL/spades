
template <typename PType, int NStructReal, int NStructInt>
SpadesParticleContainer<PType, NStructReal, NStructInt>::
    SpadesParticleContainer(amrex::AmrParGDB* par_gdb, int ngrow)
    : amrex::NeighborParticleContainer<NStructReal, NStructInt>(par_gdb, ngrow)
    , m_ngrow(ngrow)
{
    const int nlevs_max = par_gdb->maxLevel() + 1;

    if (nlevs_max > 1) {
        amrex::Abort(
            "spades::SPADES::SpadesParticleContainer::"
            "SpadesParticleContainer(): not supporting multilevel right "
            "now");
    }
}

template <typename PType, int NStructReal, int NStructInt>
SpadesParticleContainer<PType, NStructReal, NStructInt>::
    SpadesParticleContainer(
        const amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<amrex::DistributionMapping>& dmap,
        const amrex::Vector<amrex::BoxArray>& ba,
        int ngrow)
    : amrex::NeighborParticleContainer<NStructReal, NStructInt>(
          geom, dmap, ba, {2}, ngrow)
    , m_ngrow(ngrow)
{
    if (geom.size() > 1) {
        amrex::Abort(
            "spades::SPADES::SpadesParticleContainer::"
            "SpadesParticleContainer(): not supporting multilevel right "
            "now");
    }
}

template <typename PType, int NStructReal, int NStructInt>
void SpadesParticleContainer<PType, NStructReal, NStructInt>::initialize_state()
{
    BL_PROFILE("spades::SpadesParticleContainer::initialize_state()");

    m_counts.define(
        this->ParticleBoxArray(LEV), this->ParticleDistributionMap(LEV),
        PType::NTYPES, m_ngrow, amrex::MFInfo());

    m_offsets.define(
        this->ParticleBoxArray(LEV), this->ParticleDistributionMap(LEV),
        PType::NTYPES, m_ngrow, amrex::MFInfo());

    m_counts.setVal(0);
    m_offsets.setVal(0);
}

template <typename PType, int NStructReal, int NStructInt>
void SpadesParticleContainer<PType, NStructReal, NStructInt>::clear_state()
{
    BL_PROFILE("spades::SpadesParticleContainer::clear_state()");

    m_counts.clear();
    m_offsets.clear();
}
template <typename PType, int NStructReal, int NStructInt>
void SpadesParticleContainer<PType, NStructReal, NStructInt>::count_particles()
{
    BL_PROFILE("spades::SpadesParticleContainer::count_particles()");

    m_counts.setVal(0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = this->MakeMFIter(LEV); mfi.isValid(); ++mfi) {

        const amrex::Box& box = mfi.tilebox();
        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        const auto& cnt_arr = m_counts.array(mfi);
        auto& pti = this->GetParticles(LEV)[std::make_pair(gid, tid)];
        const auto& particles = pti.GetArrayOfStructs();
        const auto* pstruct = particles().dataPtr();
        const int np = pti.numParticles();

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pindex) noexcept {
            const auto& p = pstruct[pindex];
            const amrex::IntVect iv(AMREX_D_DECL(
                p.idata(CommonIntData::i), p.idata(CommonIntData::j),
                p.idata(CommonIntData::k)));

            if (box.contains(iv)) {
                amrex::Gpu::Atomic::AddNoRet(
                    &cnt_arr(iv, p.idata(CommonIntData::type_id)), 1);
            }
        });
    }
}

template <typename PType, int NStructReal, int NStructInt>
void SpadesParticleContainer<PType, NStructReal, NStructInt>::count_offsets()
{
    BL_PROFILE("spades::SpadesParticleContainer::count_offsets()");

    m_offsets.setVal(0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = this->MakeMFIter(LEV); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        const auto ncell = box.numPts();
        const auto& cnt_arr = m_counts.const_array(mfi);
        const auto& offsets_arr = m_offsets.array(mfi);
        int* p_offsets = offsets_arr.dataPtr();
        amrex::Scan::PrefixSum<int>(
            ncell,
            [=] AMREX_GPU_DEVICE(int i) -> int {
                const auto iv = box.atOffset(i);
                int total_particles = 0;
                for (int typ = 0; typ < PType::NTYPES; typ++) {
                    total_particles += cnt_arr(iv, typ);
                }
                return total_particles;
            },
            [=] AMREX_GPU_DEVICE(int i, const int& xi) { p_offsets[i] = xi; },
            amrex::Scan::Type::exclusive, amrex::Scan::noRetSum);

        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(
                     int i, int j, int AMREX_D_PICK(, , k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                for (int typ = 1; typ < PType::NTYPES; typ++) {
                    offsets_arr(iv, typ) =
                        offsets_arr(iv, typ - 1) + cnt_arr(iv, typ - 1);
                }
            });
    }
}

template <typename PType, int NStructReal, int NStructInt>
void SpadesParticleContainer<PType, NStructReal, NStructInt>::update_counts()
{
    BL_PROFILE("spades::SpadesParticleContainer::update_counts()");
    count_particles();
    count_offsets();
}

template <typename PType, int NStructReal, int NStructInt>
template <typename CompareFunc>
void SpadesParticleContainer<PType, NStructReal, NStructInt>::sort_impl(
    const CompareFunc& Compare)
{
    // Taking inspiration from AMReX's SortParticlesByBin
    BL_PROFILE("spades::SpadesParticleContainer::sort_impl()");

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = this->MakeMFIter(LEV); mfi.isValid(); ++mfi) {
        auto& particle_tile = this->ParticlesAt(LEV, mfi);
        const size_t np = particle_tile.numParticles();

        if (np == 0) {
            continue;
        }

        // BL_PROFILE_VAR(
        //     "spades::SpadesParticleContainer::sort::sort_prep",
        //     prep);
        amrex::Gpu::DeviceVector<amrex::Long> cell_list(np);
        auto* p_cell_list = cell_list.data();
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pindex) noexcept {
            p_cell_list[pindex] = pindex;
        });
        // amrex::Gpu::Device::synchronize();
        // BL_PROFILE_VAR_STOP(prep);

        // Sort particle indices based on the cell index
        // BL_PROFILE_VAR(
        //     "spades::SpadesParticleContainer::sort::sort",
        //     sort);
        const auto& particles = particle_tile.GetArrayOfStructs();
        const auto* pstruct = particles().dataPtr();
#ifdef AMREX_USE_GPU
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
        thrust::sort(
            thrust::device, cell_list.begin(), cell_list.end(),
            [=] AMREX_GPU_DEVICE(
                const amrex::Long xi, const amrex::Long yi) noexcept {
                const auto& p1 = pstruct[xi];
                const auto& p2 = pstruct[yi];
                return Compare(p1, p2);
            });
#else
        // Perform sort on CPU, then copy back to device (not good)
        amrex::Vector<amrex::Long> h_cell_list(np, 0);
        amrex::Gpu::copy(
            amrex::Gpu::deviceToHost, cell_list.begin(), cell_list.end(),
            h_cell_list.begin());
        std::sort(
            h_cell_list.begin(), h_cell_list.end(),
            [=](const amrex::Long xi, const amrex::Long yi) {
                const auto& p1 = pstruct[xi];
                const auto& p2 = pstruct[yi];
                return Compare(p1, p2);
            });
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, h_cell_list.begin(), h_cell_list.end(),
            cell_list.begin());
#endif
#else
        std::sort(
            cell_list.begin(), cell_list.end(),
            [=](const amrex::Long xi, const amrex::Long yi) {
                const auto& p1 = pstruct[xi];
                const auto& p2 = pstruct[yi];
                return Compare(p1, p2);
            });
#endif
        // amrex::Gpu::Device::synchronize();
        // BL_PROFILE_VAR_STOP(sort);

        // Reorder the particles in memory
        // BL_PROFILE_VAR(
        //     "spades::SpadesParticleContainer::sort::"
        //     "ReorderParticles",
        //     reorder);
        this->ReorderParticles(LEV, mfi, cell_list.data());
        // amrex::Gpu::Device::synchronize();
        // BL_PROFILE_VAR_STOP(reorder);
    }
    update_counts();
}

template <typename PType, int NStructReal, int NStructInt>
void SpadesParticleContainer<PType, NStructReal, NStructInt>::
    reposition_particles()
{
    BL_PROFILE("spades::SpadesParticleContainer::reposition_particles()");

    const auto& plo = this->Geom(LEV).ProbLoArray();
    const auto& dx = this->Geom(LEV).CellSizeArray();
    const auto& dxi = this->Geom(LEV).InvCellSizeArray();
    const auto& dom = this->Geom(LEV).Domain();
    const int nbins = 500;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = this->MakeMFIter(LEV); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        const auto& cnt_arr = m_counts.const_array(mfi);
        const auto& offsets_arr = m_offsets.const_array(mfi);
        const auto index = std::make_pair(gid, tid);
        auto& pti = this->GetParticles(LEV)[index];
        auto& particles = pti.GetArrayOfStructs();
        auto* pstruct = particles().dataPtr();

        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(
                     int i, int j, int AMREX_D_PICK(, , k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const auto getter = Get(iv, cnt_arr, offsets_arr, pstruct);

                for (int typ = 0; typ < PType::NTYPES; typ++) {
                    AMREX_ALWAYS_ASSERT(cnt_arr(iv, typ) < nbins);
                    for (int n = 0; n < cnt_arr(iv, typ); n++) {
                        auto& p = getter(n, typ);

                        const amrex::IntVect piv(AMREX_D_DECL(
                            p.idata(CommonIntData::i),
                            p.idata(CommonIntData::j),
                            p.idata(CommonIntData::k)));
                        AMREX_ALWAYS_ASSERT(piv == iv);

                        AMREX_D_TERM(
                            p.pos(0) = plo[0] + iv[0] * dx[0] +
                                       (typ + 1) * dx[0] / (PType::NTYPES + 1);
                            , p.pos(1) = plo[1] + iv[1] * dx[1] +
                                         (n + 1) * dx[1] / nbins;
                            , p.pos(2) =
                                  plo[2] + (iv[2] + constants::HALF) * dx[2];)

                        // ensure the particle didn't change cells
                        AMREX_ALWAYS_ASSERT(
                            piv == getParticleCell(p, plo, dxi, dom));
                    }
                }
            });
    }
}

template <typename PType, int NStructReal, int NStructInt>
void SpadesParticleContainer<PType, NStructReal, NStructInt>::write_plot_file(
    const std::string& plt_filename)
{
    BL_PROFILE("spades::SpadesParticleContainer::write_plot_file()");
    reposition_particles();
    this->WritePlotFile(
        plt_filename, identifier(), m_writeflags_real, m_writeflags_int,
        m_real_data_names, m_int_data_names);
}
