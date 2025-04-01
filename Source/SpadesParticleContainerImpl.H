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

    if (m_sort_type == "encoded") {
        m_min_timestamp.define(
            this->ParticleBoxArray(LEV), this->ParticleDistributionMap(LEV),
            PType::NTYPES, m_ngrow, amrex::MFInfo());
        m_max_timestamp.define(
            this->ParticleBoxArray(LEV), this->ParticleDistributionMap(LEV),
            PType::NTYPES, m_ngrow, amrex::MFInfo());

        m_min_timestamp.setVal(constants::LARGE_NUM);
        m_max_timestamp.setVal(0.0);
    }
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
void SpadesParticleContainer<PType, NStructReal, NStructInt>::check_sort(
    const amrex::MFIter& mfi)
{
    BL_PROFILE("spades::SpadesParticleContainer::check_sort()");
    auto& particle_tile = this->ParticlesAt(LEV, mfi);
    const size_t np = particle_tile.numParticles();
    const auto& particles = particle_tile.GetArrayOfStructs();
    const auto* pstruct = particles().dataPtr();
    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pindex) noexcept {
        if (pindex > 0) {
            const auto& pm = pstruct[pindex - 1];
            const auto& p = pstruct[pindex];
            const amrex::IntVect pivm(AMREX_D_DECL(
                pm.idata(CommonIntData::i), pm.idata(CommonIntData::j),
                pm.idata(CommonIntData::k)));
            const amrex::IntVect piv(AMREX_D_DECL(
                p.idata(CommonIntData::i), p.idata(CommonIntData::j),
                p.idata(CommonIntData::k)));
            if ((pm.idata(CommonIntData::type_id) ==
                 p.idata(CommonIntData::type_id)) &&
                (pivm == piv)) {
                if (pm.rdata(CommonRealData::timestamp) >
                    p.rdata(CommonRealData::timestamp)) {
                    printf(
                        "Unsorted time stamp %.16f > %.16f in ivm: (%d, %d) of "
                        "type "
                        "%d, iv: (%d,%d) of type %d\n",
                        pm.rdata(CommonRealData::timestamp),
                        p.rdata(CommonRealData::timestamp),
                        pm.idata(CommonIntData::i), pm.idata(CommonIntData::j),
                        pm.idata(CommonIntData::type_id),
                        p.idata(CommonIntData::i), p.idata(CommonIntData::j),
                        p.idata(CommonIntData::type_id));
                }
                AMREX_ASSERT(
                    pm.rdata(CommonRealData::timestamp) <=
                    p.rdata(CommonRealData::timestamp));
            }
        }
    });
    amrex::Gpu::streamSynchronize();
}

template <typename PType, int NStructReal, int NStructInt>
template <typename CompareFunc>
void SpadesParticleContainer<PType, NStructReal, NStructInt>::sort_impl(
    const CompareFunc& compare)
{
    BL_PROFILE("spades::SpadesParticleContainer::sort_impl()");
    if (m_sort_type == "nonencoded") {
        nonencoded_sort_impl(compare);
    } else if (m_sort_type == "encoded") {
        encoded_sort_impl();
    } else {
        amrex::Abort("Invalid sort type. Must be nonencoded or encoded");
    }
    update_counts();
}

template <typename PType, int NStructReal, int NStructInt>
template <typename CompareFunc>
void SpadesParticleContainer<PType, NStructReal, NStructInt>::
    nonencoded_sort_impl(const CompareFunc& compare)
{
    // Taking inspiration from AMReX's SortParticlesByBin
    BL_PROFILE("spades::SpadesParticleContainer::nonencoded_sort_impl()");

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = this->MakeMFIter(LEV); mfi.isValid(); ++mfi) {
        auto& particle_tile = this->ParticlesAt(LEV, mfi);
        const size_t np = particle_tile.numParticles();

        if (np == 0) {
            continue;
        }

        BL_PROFILE_VAR(
            "spades::SpadesParticleContainer::nonencoded_sort::sort_prep",
            prep);
        amrex::Gpu::DeviceVector<amrex::Long> cell_list(np);
        auto* p_cell_list = cell_list.data();
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pindex) noexcept {
            p_cell_list[pindex] = pindex;
        });
        amrex::Gpu::streamSynchronize();
        BL_PROFILE_VAR_STOP(prep);

        // Sort particle indices based on the cell index
        BL_PROFILE_VAR(
            "spades::SpadesParticleContainer::nonencoded_sort::sort", sort);
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
                return compare(p1, p2);
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
                return compare(p1, p2);
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
                return compare(p1, p2);
            });
#endif
        amrex::Gpu::streamSynchronize();
        BL_PROFILE_VAR_STOP(sort);

        // Reorder the particles in memory
        BL_PROFILE_VAR(
            "spades::SpadesParticleContainer::nonencoded_sort::"
            "ReorderParticles",
            reorder);
        this->ReorderParticles(LEV, mfi, cell_list.data());
        amrex::Gpu::streamSynchronize();
        BL_PROFILE_VAR_STOP(reorder);

#ifdef AMREX_DEBUG
        check_sort(mfi);
#endif
    }
}

template <typename PType, int NStructReal, int NStructInt>
void SpadesParticleContainer<PType, NStructReal, NStructInt>::
    encoded_sort_impl()
{
    // Taking inspiration from AMReX's SortParticlesByBin
    BL_PROFILE("spades::SpadesParticleContainer::encoded_sort_impl()");

    AMREX_ALWAYS_ASSERT(
        PType::NTYPES <= max_representation(constants::TYPE_NBITS));

    BL_PROFILE_VAR(
        "spades::SpadesParticleContainer::encoded_sort::bounds", bounds);
    m_min_timestamp.setVal(constants::LARGE_NUM);
    m_max_timestamp.setVal(0.0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = this->MakeMFIter(LEV); mfi.isValid(); ++mfi) {

        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        auto& particle_tile = this->GetParticles(LEV)[std::make_pair(gid, tid)];
        const size_t np = particle_tile.numParticles();

        if (np == 0) {
            continue;
        }
        const auto& min_ts_arr = m_min_timestamp.array(mfi);
        const auto& max_ts_arr = m_max_timestamp.array(mfi);
        const auto& particles = particle_tile.GetArrayOfStructs();
        const auto* pstruct = particles().dataPtr();
#ifdef AMREX_DEBUG
        const amrex::Box& box = mfi.tilebox();
#endif
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pindex) noexcept {
            const auto& p = pstruct[pindex];
            const amrex::IntVect piv(AMREX_D_DECL(
                p.idata(CommonIntData::i), p.idata(CommonIntData::j),
                p.idata(CommonIntData::k)));

            AMREX_ASSERT(box.contains(piv));
            amrex::Gpu::Atomic::Min(
                &min_ts_arr(piv, p.idata(CommonIntData::type_id)),
                p.rdata(CommonRealData::timestamp));
            amrex::Gpu::Atomic::Max(
                &max_ts_arr(piv, p.idata(CommonIntData::type_id)),
                p.rdata(CommonRealData::timestamp));
        });
    }
    BL_PROFILE_VAR_STOP(bounds);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = this->MakeMFIter(LEV); mfi.isValid(); ++mfi) {

        const amrex::Box& box = mfi.tilebox();
        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        auto& particle_tile = this->GetParticles(LEV)[std::make_pair(gid, tid)];
        const size_t np = particle_tile.numParticles();

        if (np == 0) {
            continue;
        }

        BL_PROFILE_VAR(
            "spades::SpadesParticleContainer::encoded_sort::sort_prep", prep);
        amrex::Gpu::DeviceVector<amrex::Long> cell_list(np);
        auto* p_cell_list = cell_list.data();
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pindex) noexcept {
            p_cell_list[pindex] = pindex;
        });
        amrex::Gpu::streamSynchronize();
        BL_PROFILE_VAR_STOP(prep);

        BL_PROFILE_VAR(
            "spades::SpadesParticleContainer::encoded_sort::encode", encode);
        const auto& particles = particle_tile.GetArrayOfStructs();
        const auto* pstruct = particles().dataPtr();
        const auto& min_ts_arr = m_min_timestamp.array(mfi);
        const auto& max_ts_arr = m_max_timestamp.array(mfi);
        amrex::Gpu::DeviceVector<std::uint64_t> encode(np);
        auto* p_encode = encode.data();
        const auto box_lo = box.smallEnd();
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pindex) noexcept {
            const auto& p = pstruct[pindex];
            AMREX_D_TERM(
                const auto pi = p.idata(CommonIntData::i) - box_lo[0];
                , const auto pj = p.idata(CommonIntData::j) - box_lo[1];
                , const auto pk = p.idata(CommonIntData::k) - box_lo[2];);

            AMREX_D_EXPR(
                AMREX_ALWAYS_ASSERT(pi >= 0), AMREX_ALWAYS_ASSERT(pj >= 0),
                AMREX_ALWAYS_ASSERT(pk >= 0));
            AMREX_D_EXPR(
                AMREX_ALWAYS_ASSERT(
                    static_cast<std::uint64_t>(pi) <=
                    max_representation(constants::I_NBITS)),
                AMREX_ALWAYS_ASSERT(
                    static_cast<std::uint64_t>(pj) <=
                    max_representation(constants::J_NBITS)),
                AMREX_ALWAYS_ASSERT(
                    static_cast<std::uint64_t>(pk) <=
                    max_representation(constants::K_NBITS)));
            AMREX_ASSERT(p.idata(CommonIntData::type_id) >= 0);
            AMREX_ASSERT(
                p.rdata(CommonRealData::timestamp) >=
                (0.0 - constants::SMALL_NUM));

            int shift = constants::TOTAL_NBITS;
#if AMREX_SPACEDIM == 3
            shift -= constants::K_NBITS;
            const std::uint64_t k =
                static_cast<std::uint64_t>(pk & bitmask(constants::K_NBITS))
                << shift;
#endif
#if AMREX_SPACEDIM >= 2
            shift -= constants::J_NBITS;
            const std::uint64_t j =
                static_cast<std::uint64_t>(pj & bitmask(constants::J_NBITS))
                << shift;
#endif
            shift -= constants::I_NBITS;
            const std::uint64_t i =
                static_cast<std::uint64_t>(pi & bitmask(constants::I_NBITS))
                << shift;
            shift -= constants::TYPE_NBITS;
            const std::uint64_t typ = static_cast<std::uint64_t>(
                                          p.idata(CommonIntData::type_id) &
                                          bitmask(constants::TYPE_NBITS))
                                      << shift;

            const amrex::IntVect piv(AMREX_D_DECL(
                p.idata(CommonIntData::i), p.idata(CommonIntData::j),
                p.idata(CommonIntData::k)));
            const amrex::Real min_t =
                min_ts_arr(piv, p.idata(CommonIntData::type_id));
            const amrex::Real max_t =
                max_ts_arr(piv, p.idata(CommonIntData::type_id));
            AMREX_ASSERT(min_t <= max_t);
            const amrex::Real delta =
                (amrex::Math::abs(max_t - min_t) > constants::EPS)
                    ? max_t - min_t
                    : 1.0;

            const amrex::Real normalized =
                (p.rdata(CommonRealData::timestamp) - min_t) / delta;
            const amrex::Real scaled =
                normalized * static_cast<amrex::Real>((1ULL << shift) - 1);
            const auto ts = static_cast<std::uint64_t>(
                static_cast<std::uint32_t>(scaled) &
                bitmask(constants::TIMESTAMP_NBITS));
            AMREX_ASSERT((shift - constants::TIMESTAMP_NBITS) == 0);

            p_encode[pindex] = AMREX_D_PICK(i, j | i, k | j | i) | typ | ts;
            // printf(
            //     "encode at %ld is %lu with iv: (%d, %d), type: %d, time: "
            //     "%.16f, (normalized: %.16f, scaled: %.16f, encoded: %lu)\n",
            //     pindex, p_encode[pindex], piv[0], piv[1],
            //     p.idata(CommonIntData::type_id),
            //     p.rdata(CommonRealData::timestamp), normalized, scaled, t);
        });
        amrex::Gpu::streamSynchronize();
        BL_PROFILE_VAR_STOP(encode);

        // Sort particle indices based on the cell index
        BL_PROFILE_VAR(
            "spades::SpadesParticleContainer::encoded_sort::sort", sort);
#ifdef AMREX_USE_GPU
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
        thrust::sort(
            thrust::device, cell_list.begin(), cell_list.end(),
            [=] AMREX_GPU_DEVICE(
                const amrex::Long xi, const amrex::Long yi) noexcept {
                return p_encode[xi] < p_encode[yi];
            });
        // thrust::sort_by_key(
        //     thrust::device, encode.begin(), encode.end(), cell_list.begin());
#else
        // Perform sort on CPU, then copy back to device (not good)
        amrex::Vector<amrex::Long> h_cell_list(np, 0);
        amrex::Vector<std::uint64_t> h_encode(np, 0);
        amrex::Gpu::copy(
            amrex::Gpu::deviceToHost, cell_list.begin(), cell_list.end(),
            h_cell_list.begin());
        amrex::Gpu::copy(
            amrex::Gpu::deviceToHost, encode.begin(), encode.end(),
            h_encode.begin());
        std::sort(
            h_cell_list.begin(), h_cell_list.end(),
            [=](const amrex::Long xi, const amrex::Long yi) {
                return h_encode[xi] < h_encode[yi];
            });
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, h_cell_list.begin(), h_cell_list.end(),
            cell_list.begin());
#endif
#else
        std::sort(
            cell_list.begin(), cell_list.end(),
            [=](const amrex::Long xi, const amrex::Long yi) {
                return p_encode[xi] < p_encode[yi];
            });
#endif
        amrex::Gpu::streamSynchronize();
        BL_PROFILE_VAR_STOP(sort);

        // Reorder the particles in memory
        BL_PROFILE_VAR(
            "spades::SpadesParticleContainer::encoded_sort::"
            "ReorderParticles",
            reorder);
        this->ReorderParticles(LEV, mfi, cell_list.data());
        amrex::Gpu::streamSynchronize();
        BL_PROFILE_VAR_STOP(reorder);

#ifdef AMREX_DEBUG
        check_sort(mfi);
#endif
    }
}

template <typename PType, int NStructReal, int NStructInt>
void SpadesParticleContainer<PType, NStructReal, NStructInt>::print_messages(
    const std::string& header)
{
    BL_PROFILE("spades::SpadesParticleContainer::print_messages()");

    amrex::Print() << header << " -----------------------------" << std::endl;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = this->MakeMFIter(LEV); mfi.isValid(); ++mfi) {
        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        const auto index = std::make_pair(gid, tid);
        const auto& particle_tile = this->GetParticles(LEV)[index];
        const size_t np = particle_tile.numParticles();
        const auto& particles = particle_tile.GetArrayOfStructs();
        const auto* pstruct = particles().dataPtr();

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pindex) noexcept {
            const auto& p = pstruct[pindex];
            DevicePrint()(p);
        });
    }
    amrex::Print() << "End " << header << "-----------------------------"
                   << std::endl;
}

template <typename PType, int NStructReal, int NStructInt>
void SpadesParticleContainer<PType, NStructReal, NStructInt>::
    reposition_particles()
{
    BL_PROFILE("spades::SpadesParticleContainer::reposition_particles()");

    const auto& plo = this->Geom(LEV).ProbLoArray();
    const auto& dx = this->Geom(LEV).CellSizeArray();
#ifdef AMREX_DEBUG
    const auto& dxi = this->Geom(LEV).InvCellSizeArray();
    const auto& dom = this->Geom(LEV).Domain();
#endif
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
                    AMREX_ASSERT(cnt_arr(iv, typ) < nbins);
                    for (int n = 0; n < cnt_arr(iv, typ); n++) {
                        auto& p = getter(n, typ);

                        const amrex::IntVect piv(AMREX_D_DECL(
                            p.idata(CommonIntData::i),
                            p.idata(CommonIntData::j),
                            p.idata(CommonIntData::k)));
                        AMREX_ASSERT(piv == iv);

                        AMREX_D_TERM(
                            p.pos(0) = plo[0] + iv[0] * dx[0] +
                                       (typ + 1) * dx[0] / (PType::NTYPES + 1);
                            , p.pos(1) = plo[1] + iv[1] * dx[1] +
                                         (n + 1) * dx[1] / nbins;
                            , p.pos(2) =
                                  plo[2] + (iv[2] + constants::HALF) * dx[2];)

                        // ensure the particle didn't change cells
                        AMREX_ASSERT(piv == getParticleCell(p, plo, dxi, dom));
                    }
                }
            });
    }
}

template <typename PType, int NStructReal, int NStructInt>
void SpadesParticleContainer<PType, NStructReal, NStructInt>::
    write_plot_file_impl(
        const std::string& plt_filename, const std::string& name)
{
    BL_PROFILE("spades::SpadesParticleContainer::write_plot_file()");
    reposition_particles();
    this->WritePlotFile(
        plt_filename, name, m_writeflags_real, m_writeflags_int,
        m_real_data_names, m_int_data_names);
}
