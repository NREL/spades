#include <utility>
#include "CellSortedParticleContainer.H"

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#endif

namespace spades::particles {

ParticleContainerInfo::ParticleContainerInfo(std::string basename)
    : m_basename(std::move(basename))
{}

ParticleContainerInfo::~ParticleContainerInfo() = default;

CellSortedParticleContainer::CellSortedParticleContainer(
    amrex::AmrParGDB* par_gdb, int ngrow)
    : amrex::NeighborParticleContainer<RealData::ncomps, IntData::ncomps>(
          par_gdb, ngrow)
    , m_info(identifier())
    , m_ngrow(ngrow)
{
    const int nlevs_max = par_gdb->maxLevel() + 1;

    if (nlevs_max > 1) {
        amrex::Abort(
            "spades::SPADES::CellSortedParticleContainer::"
            "CellSortedParticleContainer(): not supporting multilevel right "
            "now");
    }

    initialize_vectors();
}

CellSortedParticleContainer::CellSortedParticleContainer(
    const amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<amrex::DistributionMapping>& dmap,
    const amrex::Vector<amrex::BoxArray>& ba,
    int ngrow)
    : amrex::NeighborParticleContainer<RealData::ncomps, IntData::ncomps>(
          geom, dmap, ba, {2}, ngrow)
    , m_info(identifier())
    , m_ngrow(ngrow)
{
    if (geom.size() > 1) {
        amrex::Abort(
            "spades::SPADES::CellSortedParticleContainer::"
            "CellSortedParticleContainer(): not supporting multilevel right "
            "now");
    }

    initialize_vectors();
}

void CellSortedParticleContainer::initialize_vectors()
{
    BL_PROFILE("spades::CellSortedParticleContainer::initialize_vectors()");

    m_real_data_names.resize(RealData::ncomps, "");
    m_writeflags_real.resize(RealData::ncomps, 0);
    m_int_data_names.resize(IntData::ncomps, "");
    m_writeflags_int.resize(IntData::ncomps, 0);

    m_real_data_names[RealData::timestamp] = "timestamp";
    m_writeflags_real[RealData::timestamp] = 1;
    m_real_data_names[RealData::old_timestamp] = "old_timestamp";
    m_writeflags_real[RealData::old_timestamp] = 1;

    m_int_data_names[IntData::type_id] = "type_id";
    m_writeflags_int[IntData::type_id] = 1;
    m_int_data_names[IntData::sender] = "sender";
    m_writeflags_int[IntData::sender] = 1;
    m_int_data_names[IntData::receiver] = "receiver";
    m_writeflags_int[IntData::receiver] = 1;
    m_int_data_names[IntData::receiver] = "pair";
    m_writeflags_int[IntData::receiver] = 0;

    const int nlevs_max = m_gdb->maxLevel() + 1;
    m_message_counts.resize(nlevs_max);
    m_offsets.resize(nlevs_max);
}

void CellSortedParticleContainer::initialize_state()
{
    BL_PROFILE("spades::CellSortedParticleContainer::initialize_state()");

    const int lev = 0;

    m_message_counts[lev].define(
        ParticleBoxArray(lev), ParticleDistributionMap(lev),
        MessageTypes::NTYPES, m_ngrow, amrex::MFInfo());

    m_offsets[lev].define(
        ParticleBoxArray(lev), ParticleDistributionMap(lev),
        MessageTypes::NTYPES, m_ngrow, amrex::MFInfo());

    m_message_counts[lev].setVal(0);
    m_offsets[lev].setVal(0);
}

void CellSortedParticleContainer::clear_state()
{
    BL_PROFILE("spades::CellSortedParticleContainer::clear_state()");
    const int lev = 0;
    m_message_counts[lev].clear();
    m_offsets[lev].clear();
}

void CellSortedParticleContainer::update_counts()
{
    BL_PROFILE("spades::CellSortedParticleContainer::update_counts()");
    count_messages();
    count_offsets();
}

void CellSortedParticleContainer::count_messages()
{
    BL_PROFILE("spades::CellSortedParticleContainer::count_messages()");

    const int lev = 0;
    m_message_counts[lev].setVal(0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {

        const amrex::Box& box = mfi.tilebox();
        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        const auto& cnt_arr = m_message_counts[lev].array(mfi);
        auto& pti = GetParticles(lev)[std::make_pair(gid, tid)];
        const auto& particles = pti.GetArrayOfStructs();
        const auto* pstruct = particles().dataPtr();
        const int np = pti.numParticles();

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pindex) noexcept {
            const auto& p = pstruct[pindex];
            const amrex::IntVect iv(AMREX_D_DECL(
                p.idata(IntData::i), p.idata(IntData::j), p.idata(IntData::k)));

            if (box.contains(iv)) {
                amrex::Gpu::Atomic::AddNoRet(
                    &cnt_arr(iv, p.idata(IntData::type_id)), 1);
            }
        });
    }
}

void CellSortedParticleContainer::count_offsets()
{
    BL_PROFILE("spades::CellSortedParticleContainer::count_offsets()");

    const int lev = 0;
    m_offsets[lev].setVal(0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        const auto ncell = box.numPts();
        const auto& cnt_arr = m_message_counts[lev].const_array(mfi);
        const auto& offsets_arr = m_offsets[lev].array(mfi);
        int* p_offsets = offsets_arr.dataPtr();
        amrex::Scan::PrefixSum<int>(
            ncell,
            [=] AMREX_GPU_DEVICE(int i) -> int {
                const auto iv = box.atOffset(i);
                int total_messages = 0;
                for (int typ = 0; typ < MessageTypes::NTYPES; typ++) {
                    total_messages += cnt_arr(iv, typ);
                }
                return total_messages;
            },
            [=] AMREX_GPU_DEVICE(int i, const int& x) { p_offsets[i] = x; },
            amrex::Scan::Type::exclusive, amrex::Scan::noRetSum);

        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(
                     int i, int j, int AMREX_D_PICK(, , k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                for (int typ = 1; typ < MessageTypes::NTYPES; typ++) {
                    offsets_arr(iv, typ) =
                        offsets_arr(iv, typ - 1) + cnt_arr(iv, typ - 1);
                }
            });
    }
}

void CellSortedParticleContainer::initialize_particles(
    const amrex::Real lookahead)
{
    BL_PROFILE("spades::CellSortedParticleContainer::initialize_particles()");

    const int lev = 0;

    const auto& plo = Geom(lev).ProbLoArray();
    const auto& dx = Geom(lev).CellSizeArray();
    const auto& dom = Geom(lev).Domain();
    // const auto& dlo = dom.smallEnd();
    // const auto& dhi = dom.bigEnd();

    // // Some test particles
    // #ifdef AMREX_USE_OMP
    // #pragma omp parallel if (Gpu::notInLaunchRegion())
    // #endif
    //     for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
    //         const amrex::Box& box = mfi.tilebox();
    //         const int gid = mfi.index();
    //         const int tid = mfi.LocalTileIndex();
    //         auto& pti = GetParticles(lev)[std::make_pair(gid, tid)];

    //         for (amrex::IntVect iv = box.smallEnd(); iv <= box.bigEnd();
    //              box.next(iv)) {
    //             const auto test = amrex::Random_int(dhi[0] - dlo[0] + 1) +
    //             dlo[0]; amrex::IntVect iv_src(AMREX_D_DECL(0, 0, 0));
    //             amrex::IntVect iv_dest(AMREX_D_DECL(
    //                 amrex::Random_int(dhi[0] - dlo[0] + 1) + dlo[0],
    //                 amrex::Random_int(dhi[1] - dlo[1] + 1) + dlo[1],
    //                 amrex::Random_int(dhi[2] - dlo[2] + 1) + dlo[2]));
    //             amrex::IntVect iv_dest2(AMREX_D_DECL(
    //                 amrex::Random_int(dhi[0] - dlo[0] + 1) + dlo[0],
    //                 amrex::Random_int(dhi[1] - dlo[1] + 1) + dlo[1],
    //                 amrex::Random_int(dhi[2] - dlo[2] + 1) + dlo[2]));
    //             if (iv == iv_src) {
    //                 {
    //                     ParticleType
    //                     p; p.id() =
    //                         ParticleType::NextID();
    //                     p.cpu() = amrex::ParallelDescriptor::MyProc();

    //                     p.idata(IntData::type_id) =
    //                         MessageTypes::anti_message;
    //                     p.idata(IntData::sender) =
    //                     dom.index(iv_src);
    //                     p.idata(IntData::receiver) =
    //                     dom.index(iv_dest);
    //                     p.rdata(RealData::timestamp) =
    //                         random_exponential(1.0, engine) + lookahead + 20;

    //                     AMREX_D_TERM(
    //                         p.pos(0) = plo[0] + (iv_dest[0] + 0.5) * dx[0];
    //                         , p.pos(1) = plo[1] + (iv_dest[1] + 0.5) * dx[1];
    //                         , p.pos(2) = plo[2] + (iv_dest[2] + 0.5) *
    //                         dx[2];)

    //                     AMREX_D_TERM(p.idata(IntData::i) =
    //                     iv_dest[0];
    //                                  , p.idata(IntData::j) =
    //                                  iv_dest[1]; ,
    //                                  p.idata(IntData::k) =
    //                                  iv_dest[2];)

    //                     pti.push_back(p);
    //                 }
    //                 {
    //                     ParticleType
    //                     p; p.id() =
    //                         ParticleType::NextID();
    //                     p.cpu() = amrex::ParallelDescriptor::MyProc();

    //                     p.idata(IntData::type_id) =
    //                         MessageTypes::message;
    //                     p.idata(IntData::sender) =
    //                     dom.index(iv_src);
    //                     p.idata(IntData::receiver) =
    //                     dom.index(iv_dest2);
    //                     p.rdata(RealData::timestamp) =
    //                         random_exponential(1.0, engine) + lookahead;

    //                     AMREX_D_TERM(
    //                         p.pos(0) = plo[0] + (iv_dest2[0] + 0.5) * dx[0];
    //                         , p.pos(1) = plo[1] + (iv_dest2[1] + 0.5) *
    //                         dx[1]; , p.pos(2) = plo[2] + (iv_dest2[2] + 0.5)
    //                         * dx[2];)

    //                     AMREX_D_TERM(p.idata(IntData::i) =
    //                     iv_dest2[0];
    //                                  , p.idata(IntData::j) =
    //                                  iv_dest2[1];
    //                                  ,
    //                                  p.idata(IntData::k) =
    //                                  iv_dest2[2];)

    //                     pti.push_back(p);
    //                 }
    //                 { // creating another particle
    //                     ParticleType
    //                     p; p.id() =
    //                         ParticleType::NextID();
    //                     p.cpu() = amrex::ParallelDescriptor::MyProc();

    //                     p.idata(IntData::type_id) =
    //                         MessageTypes::message;
    //                     p.idata(IntData::sender) =
    //                     dom.index(iv_src);
    //                     p.idata(IntData::receiver) =
    //                     dom.index(iv_dest);
    //                     p.rdata(RealData::timestamp) =
    //                         random_exponential(1.0, engine) + lookahead;

    //                     AMREX_D_TERM(
    //                         p.pos(0) = plo[0] + (iv_dest[0] + 0.5) * dx[0];
    //                         , p.pos(1) = plo[1] + (iv_dest[1] + 0.5) * dx[1];
    //                         , p.pos(2) = plo[2] + (iv_dest[2] + 0.5) *
    //                         dx[2];)

    //                     AMREX_D_TERM(p.idata(IntData::i) =
    //                     iv_dest[0];
    //                                  , p.idata(IntData::j) =
    //                                  iv_dest[1]; ,
    //                                  p.idata(IntData::k) =
    //                                  iv_dest[2];)

    //                     pti.push_back(p);
    //                 }
    //                 { // creating another particle
    //                     ParticleType
    //                     p; p.id() =
    //                         ParticleType::NextID();
    //                     p.cpu() = amrex::ParallelDescriptor::MyProc();

    //                     p.idata(IntData::type_id) =
    //                         MessageTypes::undefined;
    //                     p.idata(IntData::sender) =
    //                     dom.index(iv_src);
    //                     p.idata(IntData::receiver) =
    //                     dom.index(iv_dest2);
    //                     p.rdata(RealData::timestamp) =
    //                         random_exponential(1.0, engine) + lookahead;

    //                     AMREX_D_TERM(
    //                         p.pos(0) = plo[0] + (iv_dest2[0] + 0.5) * dx[0];
    //                         , p.pos(1) = plo[1] + (iv_dest2[1] + 0.5) *
    //                         dx[1]; , p.pos(2) = plo[2] + (iv_dest2[2] + 0.5)
    //                         * dx[2];)

    //                     AMREX_D_TERM(p.idata(IntData::i) =
    //                     iv_dest2[0];
    //                                  , p.idata(IntData::j) =
    //                                  iv_dest2[1];
    //                                  ,
    //                                  p.idata(IntData::k) =
    //                                  iv_dest2[2];)

    //                     pti.push_back(p);
    //                 }
    //             }
    //         }
    //     }

    const int np_per_cell = 100;
    const int msg_per_cell = 10;
    AMREX_ALWAYS_ASSERT(np_per_cell > 2 * msg_per_cell);

    amrex::iMultiFab num_particles(
        ParticleBoxArray(lev), ParticleDistributionMap(lev), 1, 0);
    amrex::iMultiFab offsets(
        ParticleBoxArray(lev), ParticleDistributionMap(lev), 1, 0);
    num_particles.setVal(np_per_cell);
    offsets.setVal(0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();

        const auto ncells = static_cast<int>(box.numPts());
        const int* in = num_particles[mfi].dataPtr();
        int* out = offsets[mfi].dataPtr();
        const auto np = amrex::Scan::PrefixSum<int>(
            ncells, [=] AMREX_GPU_DEVICE(int i) -> int { return in[i]; },
            [=] AMREX_GPU_DEVICE(int i, int const& x) { out[i] = x; },
            amrex::Scan::Type::exclusive, amrex::Scan::retSum);

        const amrex::Long pid = ParticleType::NextID();
        ParticleType::NextID(pid + np);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            static_cast<amrex::Long>(pid + np) < amrex::LastParticleID,
            "Error: overflow on particle id numbers!");

        const auto my_proc = amrex::ParallelDescriptor::MyProc();
        const auto& offset_arr = offsets[mfi].const_array();
        const auto& num_particles_arr = num_particles[mfi].const_array();
        auto& pti = DefineAndReturnParticleTile(lev, mfi);
        pti.resize(np);
        auto* aos = pti.GetArrayOfStructs().dataPtr();
        amrex::ParallelForRNG(
            box, [=] AMREX_GPU_DEVICE(
                     int i, int j, int AMREX_D_PICK(, , k),
                     amrex::RandomEngine const& engine) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const int start = offset_arr(iv);
                for (int n = start; n < start + num_particles_arr(iv); n++) {
                    auto& p = aos[n];
                    p.id() = pid + n;
                    p.cpu() = my_proc;

                    MarkUndefined()(p);
                    p.idata(IntData::sender) = static_cast<int>(dom.index(iv));
                    p.idata(IntData::receiver) =
                        static_cast<int>(dom.index(iv));

                    AMREX_D_TERM(p.pos(0) = plo[0] + (iv[0] + 0.5) * dx[0];
                                 , p.pos(1) = plo[1] + (iv[1] + 0.5) * dx[1];
                                 , p.pos(2) = plo[2] + (iv[2] + 0.5) * dx[2];)

                    AMREX_D_TERM(p.idata(IntData::i) = iv[0];
                                 , p.idata(IntData::j) = iv[1];
                                 , p.idata(IntData::k) = iv[2];)
                }

                for (int n = start; n < start + msg_per_cell; n++) {
                    auto& pmsg = aos[n];
                    const auto pair = static_cast<int>(
                        pairing_function(pmsg.cpu(), pmsg.id()));
                    const amrex::Real ts =
                        random_exponential(1.0, engine) + lookahead;

                    pmsg.rdata(RealData::timestamp) = ts;
                    pmsg.idata(IntData::type_id) = MessageTypes::MESSAGE;
                    pmsg.idata(IntData::pair) = pair;

                    auto& pcnj = aos[n + msg_per_cell];
                    pcnj.rdata(RealData::timestamp) = ts;
                    pcnj.idata(IntData::type_id) = MessageTypes::CONJUGATE;
                    pcnj.idata(IntData::pair) = pair;
                }
            });
    }
    Redistribute();

    // Sanity check all initial particles
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
        auto index = std::make_pair(pti.index(), pti.LocalTileIndex());

        auto& particle_tile = GetParticles(lev)[index];
        const size_t np = particle_tile.numParticles();
        auto& particles = particle_tile.GetArrayOfStructs();
        auto* pstruct = particles().dataPtr();

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pindex) noexcept {
            auto& p = pstruct[pindex];
            bool valid_type = false;
            for (int typ = 0; typ < MessageTypes::NTYPES; typ++) {
                valid_type = p.idata(IntData::type_id) == typ;
                if (valid_type) {
                    break;
                }
            }
            AMREX_ALWAYS_ASSERT(valid_type);
            AMREX_ALWAYS_ASSERT(p.id() > 0);
        });
    }
}

void CellSortedParticleContainer::sort_particles()
{
    // Taking inspiration from AMReX's SortParticlesByBin
    BL_PROFILE("spades::CellSortedParticleContainer::sort_particles()");

    const int lev = 0;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        auto& particle_tile = ParticlesAt(lev, mfi);
        const size_t np = particle_tile.numParticles();

        if (np == 0) {
            continue;
        }

        // BL_PROFILE_VAR(
        //     "spades::CellSortedParticleContainer::sort_particles::sort_prep",
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
        //     "spades::CellSortedParticleContainer::sort_particles::sort",
        //     sort);
        const auto& particles = particle_tile.GetArrayOfStructs();
        const auto* pstruct = particles().dataPtr();
#ifdef AMREX_USE_GPU
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
        thrust::sort(
            thrust::device, cell_list.begin(), cell_list.end(),
            [=] AMREX_GPU_DEVICE(
                const amrex::Long x, const amrex::Long y) noexcept {
                const auto& p1 = pstruct[x];
                const auto& p2 = pstruct[y];
                return Compare()(p1, p2);
            });
#else
        // Perform sort on CPU, then copy back to device (not good)
        amrex::Vector<amrex::Long> h_cell_list(np, 0);
        amrex::Gpu::copy(
            amrex::Gpu::deviceToHost, cell_list.begin(), cell_list.end(),
            h_cell_list.begin());
        std::sort(
            h_cell_list.begin(), h_cell_list.end(),
            [=](const amrex::Long x, const amrex::Long y) {
                const auto& p1 = pstruct[x];
                const auto& p2 = pstruct[y];
                return Compare()(p1, p2);
            });
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, h_cell_list.begin(), h_cell_list.end(),
            cell_list.begin());
#endif
#else
        std::sort(
            cell_list.begin(), cell_list.end(),
            [=](const amrex::Long x, const amrex::Long y) {
                const auto& p1 = pstruct[x];
                const auto& p2 = pstruct[y];
                return Compare()(p1, p2);
            });
#endif
        // amrex::Gpu::Device::synchronize();
        // BL_PROFILE_VAR_STOP(sort);

        // Reorder the particles in memory
        // BL_PROFILE_VAR(
        //     "spades::CellSortedParticleContainer::sort_particles::"
        //     "ReorderParticles",
        //     reorder);
        ReorderParticles(lev, mfi, cell_list.data());
        // amrex::Gpu::Device::synchronize();
        // BL_PROFILE_VAR_STOP(reorder);
    }
    update_counts();
}

void CellSortedParticleContainer::update_undefined()
{
    BL_PROFILE("spades::CellSortedParticleContainer::update_undefined()");

    const int lev = 0;
    const auto& plo = Geom(lev).ProbLoArray();
    const auto& dx = Geom(lev).CellSizeArray();
    const auto& dom = Geom(lev).Domain();
    const int lower_count = m_lower_undefined_count;
    const int upper_count = m_upper_undefined_count;
    const int mid_count =
        (m_upper_undefined_count - m_lower_undefined_count) / 2 +
        m_lower_undefined_count;

    CellSortedParticleContainer pc_adds(
        m_gdb->Geom(), m_gdb->DistributionMap(), m_gdb->boxArray(), ngrow());

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        const auto& cnt_arr = m_message_counts[lev].const_array(mfi);
        const auto& offsets_arr = m_offsets[lev].const_array(mfi);
        auto& particle_tile = GetParticles(lev)[std::make_pair(gid, tid)];
        auto& particles = particle_tile.GetArrayOfStructs();
        auto* pstruct = particles().dataPtr();

        // remove particles
        // BL_PROFILE_VAR(
        //     "spades::CellSortedParticleContainer::update_undefined::remove",
        //     remove);
        const auto ncells = static_cast<int>(box.numPts());
        amrex::Gpu::DeviceVector<int> removals(ncells, 0);
        auto* p_removals = removals.data();
        amrex::ParallelFor(ncells, [=] AMREX_GPU_DEVICE(long icell) noexcept {
            const auto iv = box.atOffset(icell);
            const int current_count = cnt_arr(iv, MessageTypes::UNDEFINED);
            if (current_count > upper_count) {
                p_removals[icell] = current_count - mid_count;
            }
        });
        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(
                     int i, int j, int AMREX_D_PICK(, , k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const auto idx = box.index(iv);
                const auto np = p_removals[idx];
                const auto getter = Get(iv, cnt_arr, offsets_arr, pstruct);
                for (int m = 0; m < np; m++) {
                    auto& p = getter(m, MessageTypes::UNDEFINED);
                    p.id() = -1;
                }
            });
        // amrex::Gpu::Device::synchronize();
        // BL_PROFILE_VAR_STOP(remove);

        // BL_PROFILE_VAR(
        //     "spades::CellSortedParticleContainer::update_undefined::compute_add",
        //     compute_add);
        amrex::Gpu::DeviceVector<int> additions(ncells, 0);
        auto* p_additions = additions.data();
        amrex::ParallelFor(ncells, [=] AMREX_GPU_DEVICE(long icell) noexcept {
            const auto iv = box.atOffset(icell);
            const int current_count = cnt_arr(iv, MessageTypes::UNDEFINED);
            if (lower_count > current_count) {
                p_additions[icell] = mid_count - current_count;
            }
        });

        amrex::Gpu::DeviceVector<int> offsets(ncells, 0);
        auto* p_offsets = offsets.data();
        const auto np = amrex::Scan::PrefixSum<int>(
            ncells,
            [=] AMREX_GPU_DEVICE(int i) -> int { return p_additions[i]; },
            [=] AMREX_GPU_DEVICE(int i, int const& x) { p_offsets[i] = x; },
            amrex::Scan::Type::exclusive, amrex::Scan::retSum);

        const amrex::Long pid = ParticleType::NextID();
        ParticleType::NextID(pid + np);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            static_cast<amrex::Long>(pid + np) < amrex::LastParticleID,
            "Error: overflow on particle id numbers!");

        auto& ptile_adds = pc_adds.DefineAndReturnParticleTile(lev, mfi);
        ptile_adds.resize(np);
        const auto my_proc = amrex::ParallelDescriptor::MyProc();
        auto* aos = ptile_adds.GetArrayOfStructs().dataPtr();
        amrex::ParallelFor(ncells, [=] AMREX_GPU_DEVICE(long icell) noexcept {
            const int start = p_offsets[icell];
            const auto iv = box.atOffset(icell);
            for (int n = start; n < start + p_additions[icell]; n++) {
                auto& p = aos[n];
                p.id() = pid + n;
                p.cpu() = my_proc;

                MarkUndefined()(p);
                p.idata(IntData::sender) = static_cast<int>(dom.index(iv));
                p.idata(IntData::receiver) = static_cast<int>(dom.index(iv));

                AMREX_D_TERM(p.pos(0) = plo[0] + (iv[0] + 0.5) * dx[0];
                             , p.pos(1) = plo[1] + (iv[1] + 0.5) * dx[1];
                             , p.pos(2) = plo[2] + (iv[2] + 0.5) * dx[2];)

                AMREX_D_TERM(p.idata(IntData::i) = iv[0];
                             , p.idata(IntData::j) = iv[1];
                             , p.idata(IntData::k) = iv[2];)
            }
        });
        // amrex::Gpu::Device::synchronize();
        // BL_PROFILE_VAR_STOP(compute_add);
    }

    addParticles(pc_adds, true);
}

void CellSortedParticleContainer::resolve_pairs()
{
    BL_PROFILE("spades::CellSortedParticleContainer::resolve_pairs()");

    const int lev = 0;
    const auto& dom = Geom(lev).Domain();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        const auto& cnt_arr = m_message_counts[lev].const_array(mfi);
        const auto& offsets_arr = m_offsets[lev].const_array(mfi);
        auto& particle_tile = GetParticles(lev)[std::make_pair(gid, tid)];
        auto& particles = particle_tile.GetArrayOfStructs();
        auto* pstruct = particles().dataPtr();

        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(
                     int i, int j, int AMREX_D_PICK(, , k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const auto getter = Get(iv, cnt_arr, offsets_arr, pstruct);

                for (int n = 0; n < cnt_arr(iv, MessageTypes::ANTI_MESSAGE);
                     n++) {

                    auto& pant = getter(n, MessageTypes::ANTI_MESSAGE);
                    AMREX_ALWAYS_ASSERT(pant.idata(IntData::pair) != -1);
                    AMREX_ALWAYS_ASSERT(
                        pant.idata(IntData::receiver) == dom.index(iv));

                    bool found_pair = false;
                    for (int m = 0; m < cnt_arr(iv, MessageTypes::MESSAGE);
                         m++) {
                        // This is a message that was already treated,
                        // expect it to be undefined
                        if (!getter.check(m, MessageTypes::MESSAGE)) {
                            getter.assert_different(
                                m, MessageTypes::MESSAGE,
                                MessageTypes::UNDEFINED);
                            continue;
                        }
                        auto& pmsg = getter(m, MessageTypes::MESSAGE);
                        if ((pmsg.idata(IntData::pair) ==
                             pant.idata(IntData::pair)) &&
                            (std::abs(
                                 pmsg.rdata(RealData::timestamp) -
                                 pant.rdata(RealData::timestamp)) <
                             constants::EPS)) {
                            AMREX_ALWAYS_ASSERT(
                                pmsg.idata(IntData::sender) ==
                                pant.idata(IntData::sender));
                            AMREX_ALWAYS_ASSERT(
                                pmsg.idata(IntData::receiver) ==
                                pant.idata(IntData::receiver));
                            MarkUndefined()(pant);
                            MarkUndefined()(pmsg);
                            found_pair = true;
                            break;
                        }
                    }
                    AMREX_ALWAYS_ASSERT(found_pair);
                }
            });
    }

    sort_particles();
}

void CellSortedParticleContainer::garbage_collect(const amrex::Real gvt)
{
    BL_PROFILE("spades::CellSortedParticleContainer::garbage_collect()");

    const int lev = 0;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
        auto index = std::make_pair(pti.index(), pti.LocalTileIndex());

        auto& particle_tile = GetParticles(lev)[index];
        const size_t np = particle_tile.numParticles();
        auto& particles = particle_tile.GetArrayOfStructs();
        auto* pstruct = particles().dataPtr();

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pindex) noexcept {
            auto& p = pstruct[pindex];
            if ((p.rdata(RealData::timestamp) < gvt) &&
                (p.idata(IntData::type_id) != MessageTypes::UNDEFINED)) {
                MarkUndefined()(p);
            }
        });
    }
}

void CellSortedParticleContainer::reposition_messages()
{
    BL_PROFILE("spades::CellSortedParticleContainer::reposition_messages()");

    const int lev = 0;
    const auto& plo = Geom(lev).ProbLoArray();
    const auto& dx = Geom(lev).CellSizeArray();
    const auto& dxi = Geom(lev).InvCellSizeArray();
    const auto& dom = Geom(lev).Domain();
    const int nbins = 500;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        const auto& cnt_arr = m_message_counts[lev].const_array(mfi);
        const auto& offsets_arr = m_offsets[lev].const_array(mfi);
        const auto index = std::make_pair(gid, tid);
        auto& pti = GetParticles(lev)[index];
        auto& particles = pti.GetArrayOfStructs();
        auto* pstruct = particles().dataPtr();

        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(
                     int i, int j, int AMREX_D_PICK(, , k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const auto getter = Get(iv, cnt_arr, offsets_arr, pstruct);

                for (int typ = 0; typ < MessageTypes::NTYPES; typ++) {
                    AMREX_ALWAYS_ASSERT(cnt_arr(iv, typ) < nbins);
                    for (int n = 0; n < cnt_arr(iv, typ); n++) {
                        auto& p = getter(n, typ);

                        const amrex::IntVect piv(AMREX_D_DECL(
                            p.idata(IntData::i), p.idata(IntData::j),
                            p.idata(IntData::k)));
                        AMREX_ALWAYS_ASSERT(piv == iv);

                        AMREX_D_TERM(p.pos(0) = plo[0] + iv[0] * dx[0] +
                                                (typ + 1) * dx[0] /
                                                    (MessageTypes::NTYPES + 1);
                                     , p.pos(1) = plo[1] + iv[1] * dx[1] +
                                                  (n + 1) * dx[1] / nbins;
                                     ,
                                     p.pos(2) = plo[2] + (iv[2] + 0.5) * dx[2];)

                        // ensure the particle didn't change cells
                        AMREX_ALWAYS_ASSERT(
                            piv == getParticleCell(p, plo, dxi, dom));
                    }
                }
            });
    }
}

void CellSortedParticleContainer::write_plot_file(
    const std::string& plt_filename)
{
    BL_PROFILE("spades::CellSortedParticleContainer::write_plot_file()");
    reposition_messages();
    WritePlotFile(
        plt_filename, identifier(), m_writeflags_real, m_writeflags_int,
        m_real_data_names, m_int_data_names);
}
} // namespace spades::particles
