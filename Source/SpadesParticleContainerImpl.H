template <
    int NType,
    int NStructReal,
    int NStructInt,
    int NArrayReal,
    int NArrayInt>
SpadesParticleContainer<NType, NStructReal, NStructInt, NArrayReal, NArrayInt>::
    SpadesParticleContainer(amrex::AmrParGDB* par_gdb, int ngrow)
    : amrex::NeighborParticleContainer<
          NStructReal,
          NStructInt,
          NArrayReal,
          NArrayInt>(par_gdb, ngrow)
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

template <
    int NType,
    int NStructReal,
    int NStructInt,
    int NArrayReal,
    int NArrayInt>
SpadesParticleContainer<NType, NStructReal, NStructInt, NArrayReal, NArrayInt>::
    SpadesParticleContainer(
        const amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<amrex::DistributionMapping>& dmap,
        const amrex::Vector<amrex::BoxArray>& ba,
        int ngrow)
    : amrex::NeighborParticleContainer<
          NStructReal,
          NStructInt,
          NArrayReal,
          NArrayInt>(geom, dmap, ba, {2}, ngrow)
    , m_ngrow(ngrow)
{
    if (geom.size() > 1) {
        amrex::Abort(
            "spades::SPADES::SpadesParticleContainer::"
            "SpadesParticleContainer(): not supporting multilevel right "
            "now");
    }
}

template <
    int NType,
    int NStructReal,
    int NStructInt,
    int NArrayReal,
    int NArrayInt>
void SpadesParticleContainer<
    NType,
    NStructReal,
    NStructInt,
    NArrayReal,
    NArrayInt>::initialize_state()
{
    BL_PROFILE("spades::SpadesParticleContainer::initialize_state()");

    m_counts.define(
        this->ParticleBoxArray(LEV), this->ParticleDistributionMap(LEV), NType,
        m_ngrow, amrex::MFInfo());

    m_offsets.define(
        this->ParticleBoxArray(LEV), this->ParticleDistributionMap(LEV), NType,
        m_ngrow, amrex::MFInfo());

    m_counts.setVal(0);
    m_offsets.setVal(0);

    if (m_sort_type == "encoded") {
        m_min_timestamp.define(
            this->ParticleBoxArray(LEV), this->ParticleDistributionMap(LEV),
            NType, m_ngrow, amrex::MFInfo());
        m_max_timestamp.define(
            this->ParticleBoxArray(LEV), this->ParticleDistributionMap(LEV),
            NType, m_ngrow, amrex::MFInfo());

        m_min_timestamp.setVal(constants::LARGE_NUM);
        m_max_timestamp.setVal(0.0);
    }
}

template <
    int NType,
    int NStructReal,
    int NStructInt,
    int NArrayReal,
    int NArrayInt>
void SpadesParticleContainer<
    NType,
    NStructReal,
    NStructInt,
    NArrayReal,
    NArrayInt>::clear_state()
{
    BL_PROFILE("spades::SpadesParticleContainer::clear_state()");

    m_counts.clear();
    m_offsets.clear();
}
template <
    int NType,
    int NStructReal,
    int NStructInt,
    int NArrayReal,
    int NArrayInt>
void SpadesParticleContainer<
    NType,
    NStructReal,
    NStructInt,
    NArrayReal,
    NArrayInt>::count_particles()
{
    BL_PROFILE("spades::SpadesParticleContainer::count_particles()");

    m_counts.setVal(0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = this->MakeMFIter(LEV); mfi.isValid(); ++mfi) {

        const amrex::Box& box = mfi.tilebox();
        const auto& cnt_arr = m_counts.array(mfi);
        const auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
        auto& pti = this->GetParticles(LEV)[index];
        const auto parrs = particle_arrays(pti);
        const int np = pti.numParticles();

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pidx) noexcept {
            const amrex::IntVect iv(AMREX_D_DECL(
                parrs.m_idata[CommonIntData::i][pidx],
                parrs.m_idata[CommonIntData::j][pidx],
                parrs.m_idata[CommonIntData::k][pidx]));

            if (box.contains(iv)) {
                amrex::Gpu::Atomic::AddNoRet(
                    &cnt_arr(iv, parrs.m_idata[CommonIntData::type_id][pidx]),
                    1);
            }
        });
    }
}

template <
    int NType,
    int NStructReal,
    int NStructInt,
    int NArrayReal,
    int NArrayInt>
void SpadesParticleContainer<
    NType,
    NStructReal,
    NStructInt,
    NArrayReal,
    NArrayInt>::count_offsets()
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
                for (int typ = 0; typ < NType; typ++) {
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
                for (int typ = 1; typ < NType; typ++) {
                    offsets_arr(iv, typ) =
                        offsets_arr(iv, typ - 1) + cnt_arr(iv, typ - 1);
                }
            });
    }
}

template <
    int NType,
    int NStructReal,
    int NStructInt,
    int NArrayReal,
    int NArrayInt>
void SpadesParticleContainer<
    NType,
    NStructReal,
    NStructInt,
    NArrayReal,
    NArrayInt>::update_counts()
{
    BL_PROFILE("spades::SpadesParticleContainer::update_counts()");
    count_particles();
    count_offsets();
}

template <
    int NType,
    int NStructReal,
    int NStructInt,
    int NArrayReal,
    int NArrayInt>
void SpadesParticleContainer<
    NType,
    NStructReal,
    NStructInt,
    NArrayReal,
    NArrayInt>::check_sort(const amrex::MFIter& mfi)
{
    BL_PROFILE("spades::SpadesParticleContainer::check_sort()");
    const auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
    auto& pti = this->GetParticles(LEV)[index];
    const size_t np = pti.numParticles();
    const auto parrs = particle_arrays(pti);

    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pidx) noexcept {
        if (pidx > 0) {
            const auto pm = pidx - 1;
            const amrex::IntVect pivm(AMREX_D_DECL(
                parrs.m_idata[CommonIntData::i][pm],
                parrs.m_idata[CommonIntData::j][pm],
                parrs.m_idata[CommonIntData::k][pm]));
            const amrex::IntVect piv(AMREX_D_DECL(
                parrs.m_idata[CommonIntData::i][pidx],
                parrs.m_idata[CommonIntData::j][pidx],
                parrs.m_idata[CommonIntData::k][pidx]));
            if ((parrs.m_idata[CommonIntData::type_id][pm] ==
                 parrs.m_idata[CommonIntData::type_id][pidx]) &&
                (pivm == piv)) {
                if (parrs.m_rdata[CommonRealData::timestamp][pm] >
                    parrs.m_rdata[CommonRealData::timestamp][pidx]) {
                    printf(
                        "Unsorted time stamp %.16f > %.16f in ivm: (%d, %d) of "
                        "type "
                        "%d, iv: (%d,%d) of type %d\n",
                        parrs.m_rdata[CommonRealData::timestamp][pm],
                        parrs.m_rdata[CommonRealData::timestamp][pidx],
                        parrs.m_idata[CommonIntData::i][pm],
                        parrs.m_idata[CommonIntData::j][pm],
                        parrs.m_idata[CommonIntData::type_id][pm],
                        parrs.m_idata[CommonIntData::i][pidx],
                        parrs.m_idata[CommonIntData::j][pidx],
                        parrs.m_idata[CommonIntData::type_id][pidx]);
                }
                AMREX_ASSERT(
                    parrs.m_rdata[CommonRealData::timestamp][pm] <=
                    parrs.m_rdata[CommonRealData::timestamp][pidx]);
            }
        }
    });
    amrex::Gpu::streamSynchronize();
}

template <
    int NType,
    int NStructReal,
    int NStructInt,
    int NArrayReal,
    int NArrayInt>
template <typename CompareFunc>
void SpadesParticleContainer<
    NType,
    NStructReal,
    NStructInt,
    NArrayReal,
    NArrayInt>::sort_impl(const CompareFunc& compare)
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

template <
    int NType,
    int NStructReal,
    int NStructInt,
    int NArrayReal,
    int NArrayInt>
template <typename CompareFunc>
void SpadesParticleContainer<
    NType,
    NStructReal,
    NStructInt,
    NArrayReal,
    NArrayInt>::nonencoded_sort_impl(const CompareFunc& compare)
{
    // Taking inspiration from AMReX's SortParticlesByBin
    BL_PROFILE("spades::SpadesParticleContainer::nonencoded_sort_impl()");

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = this->MakeMFIter(LEV); mfi.isValid(); ++mfi) {
        const auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
        auto& pti = this->GetParticles(LEV)[index];
        const size_t np = pti.numParticles();

        if (np == 0) {
            continue;
        }

        BL_PROFILE_VAR(
            "spades::SpadesParticleContainer::nonencoded_sort::sort_prep",
            prep);
        amrex::Gpu::DeviceVector<amrex::Long> cell_list(np);
        auto* p_cell_list = cell_list.data();
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pidx) noexcept {
            p_cell_list[pidx] = pidx;
        });
        amrex::Gpu::streamSynchronize();
        BL_PROFILE_VAR_STOP(prep);

        // Sort particle indices based on the cell index
        BL_PROFILE_VAR(
            "spades::SpadesParticleContainer::nonencoded_sort::sort", sort);
        const auto parrs = particle_arrays(pti);
#ifdef AMREX_USE_GPU
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
        thrust::sort(
            thrust::device, cell_list.begin(), cell_list.end(),
            [=] AMREX_GPU_DEVICE(
                const amrex::Long xi, const amrex::Long yi) noexcept {
                return compare(xi, yi, parrs);
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
                return compare(xi, yi, parrs);
            });
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, h_cell_list.begin(), h_cell_list.end(),
            cell_list.begin());
#endif
#else
        std::sort(
            cell_list.begin(), cell_list.end(),
            [=](const amrex::Long xi, const amrex::Long yi) {
                return compare(xi, yi, parrs);
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

template <
    int NType,
    int NStructReal,
    int NStructInt,
    int NArrayReal,
    int NArrayInt>
void SpadesParticleContainer<
    NType,
    NStructReal,
    NStructInt,
    NArrayReal,
    NArrayInt>::encoded_sort_impl()
{
    // Taking inspiration from AMReX's SortParticlesByBin
    BL_PROFILE("spades::SpadesParticleContainer::encoded_sort_impl()");

    AMREX_ALWAYS_ASSERT(NType <= max_representation(constants::TYPE_NBITS));

    BL_PROFILE_VAR(
        "spades::SpadesParticleContainer::encoded_sort::bounds", bounds);
    m_min_timestamp.setVal(constants::LARGE_NUM);
    m_max_timestamp.setVal(0.0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = this->MakeMFIter(LEV); mfi.isValid(); ++mfi) {

        const auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
        auto& pti = this->GetParticles(LEV)[index];
        const size_t np = pti.numParticles();

        if (np == 0) {
            continue;
        }
        const auto& min_ts_arr = m_min_timestamp.array(mfi);
        const auto& max_ts_arr = m_max_timestamp.array(mfi);
        const auto parrs = particle_arrays(pti);
        auto* ts = parrs.m_rdata[CommonRealData::timestamp];
#ifdef AMREX_DEBUG
        const amrex::Box& box = mfi.tilebox();
#endif
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pidx) noexcept {
            const amrex::IntVect piv(AMREX_D_DECL(
                parrs.m_idata[CommonIntData::i][pidx],
                parrs.m_idata[CommonIntData::j][pidx],
                parrs.m_idata[CommonIntData::k][pidx]));

            AMREX_ASSERT(box.contains(piv));
            amrex::Gpu::Atomic::Min(
                &min_ts_arr(piv, parrs.m_idata[CommonIntData::type_id][pidx]),
                ts[pidx]);
            amrex::Gpu::Atomic::Max(
                &max_ts_arr(piv, parrs.m_idata[CommonIntData::type_id][pidx]),
                ts[pidx]);
        });
    }
    BL_PROFILE_VAR_STOP(bounds);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = this->MakeMFIter(LEV); mfi.isValid(); ++mfi) {

        const amrex::Box& box = mfi.tilebox();
        const auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
        auto& pti = this->GetParticles(LEV)[index];
        const size_t np = pti.numParticles();

        if (np == 0) {
            continue;
        }

        BL_PROFILE_VAR(
            "spades::SpadesParticleContainer::encoded_sort::sort_prep", prep);
        amrex::Gpu::DeviceVector<amrex::Long> cell_list(np);
        auto* p_cell_list = cell_list.data();
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pidx) noexcept {
            p_cell_list[pidx] = pidx;
        });
        amrex::Gpu::streamSynchronize();
        BL_PROFILE_VAR_STOP(prep);

        BL_PROFILE_VAR(
            "spades::SpadesParticleContainer::encoded_sort::encode", encode);
        const auto parrs = particle_arrays(pti);
        auto* ts = parrs.m_rdata[CommonRealData::timestamp];
        const auto& min_ts_arr = m_min_timestamp.array(mfi);
        const auto& max_ts_arr = m_max_timestamp.array(mfi);
        amrex::Gpu::DeviceVector<std::uint64_t> encode(np);
        auto* p_encode = encode.data();
        const auto box_lo = box.smallEnd();
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pidx) noexcept {
            AMREX_D_TERM(
                const auto pi =
                    parrs.m_idata[CommonIntData::i][pidx] - box_lo[0];
                , const auto pj =
                      parrs.m_idata[CommonIntData::j][pidx] - box_lo[1];
                , const auto pk =
                      parrs.m_idata[CommonIntData::k][pidx] - box_lo[2];);

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
            AMREX_ASSERT(parrs.m_idata[CommonIntData::type_id][pidx] >= 0);
            AMREX_ASSERT(ts[pidx] >= (0.0 - constants::SMALL_NUM));

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
            const std::uint64_t typ =
                static_cast<std::uint64_t>(
                    parrs.m_idata[CommonIntData::type_id][pidx] &
                    bitmask(constants::TYPE_NBITS))
                << shift;

            const amrex::IntVect piv(AMREX_D_DECL(
                parrs.m_idata[CommonIntData::i][pidx],
                parrs.m_idata[CommonIntData::j][pidx],
                parrs.m_idata[CommonIntData::k][pidx]));
            const amrex::Real min_t =
                min_ts_arr(piv, parrs.m_idata[CommonIntData::type_id][pidx]);
            const amrex::Real max_t =
                max_ts_arr(piv, parrs.m_idata[CommonIntData::type_id][pidx]);
            AMREX_ASSERT(min_t <= max_t);
            const amrex::Real delta =
                (amrex::Math::abs(max_t - min_t) > constants::EPS)
                    ? max_t - min_t
                    : 1.0;

            const amrex::Real normalized = (ts[pidx] - min_t) / delta;
            const amrex::Real scaled =
                normalized * static_cast<amrex::Real>((1ULL << shift) - 1);
            const auto ts = static_cast<std::uint64_t>(
                static_cast<std::uint32_t>(scaled) &
                bitmask(constants::TIMESTAMP_NBITS));
            AMREX_ASSERT((shift - constants::TIMESTAMP_NBITS) == 0);

            p_encode[pidx] = AMREX_D_PICK(i, j | i, k | j | i) | typ | ts;
            // printf(
            //     "encode at %ld is %lu with iv: (%d, %d), type: %d, time: "
            //     "%.16f, (normalized: %.16f, scaled: %.16f, encoded: %lu)\n",
            //     pidx, p_encode[pidx], piv[0], piv[1],
            //     idata(CommonIntData::type_id)[pidx],
            //     ts[pidx], normalized, scaled, t);
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

template <
    int NType,
    int NStructReal,
    int NStructInt,
    int NArrayReal,
    int NArrayInt>
void SpadesParticleContainer<
    NType,
    NStructReal,
    NStructInt,
    NArrayReal,
    NArrayInt>::print_messages(const std::string& header)
{
    BL_PROFILE("spades::SpadesParticleContainer::print_messages()");

    amrex::Print() << header << " -----------------------------" << std::endl;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = this->MakeMFIter(LEV); mfi.isValid(); ++mfi) {
        const auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
        const auto& pti = this->GetParticles(LEV)[index];
        const size_t np = pti.numParticles();
        const auto parrs = particle_arrays(pti);

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pidx) noexcept {
            DevicePrint()(pidx, parrs);
        });
    }
    amrex::Print() << "End " << header << "-----------------------------"
                   << std::endl;
}

template <
    int NType,
    int NStructReal,
    int NStructInt,
    int NArrayReal,
    int NArrayInt>
void SpadesParticleContainer<
    NType,
    NStructReal,
    NStructInt,
    NArrayReal,
    NArrayInt>::reposition_particles()
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
        const auto& cnt_arr = m_counts.const_array(mfi);
        const auto& offsets_arr = m_offsets.const_array(mfi);
        const auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
        auto& pti = this->GetParticles(LEV)[index];
        const auto parrs = particle_arrays(pti);

        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(
                     int i, int j, int AMREX_D_PICK(, , k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const auto getter = Get(iv, cnt_arr, offsets_arr, parrs);

                for (int typ = 0; typ < NType; typ++) {
                    AMREX_ASSERT(cnt_arr(iv, typ) < nbins);
                    for (int n = 0; n < cnt_arr(iv, typ); n++) {
                        const auto pidx = getter(n, typ);
                        auto& p = parrs.m_aos[pidx];

                        const amrex::IntVect piv(AMREX_D_DECL(
                            parrs.m_idata[CommonIntData::i][pidx],
                            parrs.m_idata[CommonIntData::j][pidx],
                            parrs.m_idata[CommonIntData::k][pidx]));
                        AMREX_ASSERT(piv == iv);

                        AMREX_D_TERM(p.pos(0) = plo[0] + iv[0] * dx[0] +
                                                (typ + 1) * dx[0] / (NType + 1);
                                     , p.pos(1) = plo[1] + iv[1] * dx[1] +
                                                  (n + 1) * dx[1] / nbins;
                                     , p.pos(2) =
                                           plo[2] +
                                           (iv[2] + constants::HALF) * dx[2];)

                        // ensure the particle didn't change cells
                        AMREX_ASSERT(piv == getParticleCell(p, plo, dxi, dom));
                    }
                }
            });
    }
}

template <
    int NType,
    int NStructReal,
    int NStructInt,
    int NArrayReal,
    int NArrayInt>
void SpadesParticleContainer<
    NType,
    NStructReal,
    NStructInt,
    NArrayReal,
    NArrayInt>::
    write_plot_file_impl(
        const std::string& plt_filename, const std::string& name)
{
    BL_PROFILE("spades::SpadesParticleContainer::write_plot_file()");
    reposition_particles();
    this->WritePlotFile(
        plt_filename, name, m_writeflags_real, m_writeflags_int,
        m_real_data_names, m_int_data_names);
}
