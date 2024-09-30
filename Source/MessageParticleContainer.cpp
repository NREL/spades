#include <utility>
#include "MessageParticleContainer.H"

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#endif

namespace spades::particles {

MessageParticleContainerInfo::MessageParticleContainerInfo(std::string basename)
    : m_basename(std::move(basename))
{}

MessageParticleContainerInfo::~MessageParticleContainerInfo() = default;

MessageParticleContainer::MessageParticleContainer(
    amrex::AmrParGDB* par_gdb, int ngrow)
    : amrex::NeighborParticleContainer<RealData::ncomps, IntData::ncomps>(
          par_gdb, ngrow)
    , m_info(identifier())
    , m_ngrow(ngrow)
{
    const int nlevs_max = par_gdb->maxLevel() + 1;

    if (nlevs_max > 1) {
        amrex::Abort(
            "spades::SPADES::MessageParticleContainer::"
            "MessageParticleContainer(): not supporting multilevel right "
            "now");
    }

    initialize_vectors();
}

MessageParticleContainer::MessageParticleContainer(
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
            "spades::SPADES::MessageParticleContainer::"
            "MessageParticleContainer(): not supporting multilevel right "
            "now");
    }

    initialize_vectors();
}

void MessageParticleContainer::initialize_vectors()
{
    BL_PROFILE("spades::MessageParticleContainer::initialize_vectors()");

    m_real_data_names.resize(RealData::ncomps, "");
    m_writeflags_real.resize(RealData::ncomps, 0);
    m_int_data_names.resize(IntData::ncomps, "");
    m_writeflags_int.resize(IntData::ncomps, 0);

    m_real_data_names[RealData::timestamp] = "timestamp";
    m_writeflags_real[RealData::timestamp] = 1;
    m_real_data_names[RealData::old_timestamp] = "old_timestamp";
    m_writeflags_real[RealData::old_timestamp] = 1;
    m_real_data_names[RealData::creation_time] = "creation_time";
    m_writeflags_real[RealData::creation_time] = 1;

    m_int_data_names[IntData::type_id] = "type_id";
    m_writeflags_int[IntData::type_id] = 1;
    m_int_data_names[IntData::sender] = "sender";
    m_writeflags_int[IntData::sender] = 1;
    m_int_data_names[IntData::receiver] = "receiver";
    m_writeflags_int[IntData::receiver] = 1;
    m_int_data_names[IntData::receiver] = "pair";
    m_writeflags_int[IntData::receiver] = 0;
}

void MessageParticleContainer::initialize_state()
{
    BL_PROFILE("spades::MessageParticleContainer::initialize_state()");

    m_message_counts.define(
        ParticleBoxArray(LEV), ParticleDistributionMap(LEV),
        MessageTypes::NTYPES, m_ngrow, amrex::MFInfo());

    m_offsets.define(
        ParticleBoxArray(LEV), ParticleDistributionMap(LEV),
        MessageTypes::NTYPES, m_ngrow, amrex::MFInfo());

    m_message_counts.setVal(0);
    m_offsets.setVal(0);
}

void MessageParticleContainer::clear_state()
{
    BL_PROFILE("spades::MessageParticleContainer::clear_state()");

    m_message_counts.clear();
    m_offsets.clear();
}

void MessageParticleContainer::update_counts()
{
    BL_PROFILE("spades::MessageParticleContainer::update_counts()");
    count_messages();
    count_offsets();
}

void MessageParticleContainer::count_messages()
{
    BL_PROFILE("spades::MessageParticleContainer::count_messages()");

    m_message_counts.setVal(0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = MakeMFIter(LEV); mfi.isValid(); ++mfi) {

        const amrex::Box& box = mfi.tilebox();
        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        const auto& cnt_arr = m_message_counts.array(mfi);
        auto& pti = GetParticles(LEV)[std::make_pair(gid, tid)];
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

void MessageParticleContainer::count_offsets()
{
    BL_PROFILE("spades::MessageParticleContainer::count_offsets()");

    m_offsets.setVal(0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = MakeMFIter(LEV); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        const auto ncell = box.numPts();
        const auto& cnt_arr = m_message_counts.const_array(mfi);
        const auto& offsets_arr = m_offsets.array(mfi);
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
            [=] AMREX_GPU_DEVICE(int i, const int& xi) { p_offsets[i] = xi; },
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

void MessageParticleContainer::initialize_messages(const amrex::Real lookahead)
{
    BL_PROFILE("spades::MessageParticleContainer::initialize_messages()");

    const auto& plo = Geom(LEV).ProbLoArray();
    const auto& dx = Geom(LEV).CellSizeArray();
    const auto& dom = Geom(LEV).Domain();
    // const auto& dlo = dom.smallEnd();
    // const auto& dhi = dom.bigEnd();

    // // Some test particles
    // #ifdef AMREX_USE_OMP
    // #pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
    // #endif
    //     for (amrex::MFIter mfi = MakeMFIter(LEV); mfi.isValid(); ++mfi) {
    //         const amrex::Box& box = mfi.tilebox();
    //         const int gid = mfi.index();
    //         const int tid = mfi.LocalTileIndex();
    //         auto& pti = GetParticles(LEV)[std::make_pair(gid, tid)];

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
    //                         p.pos(0) = plo[0] + (iv_dest[0] +
    //                         constants::HALF) * dx[0]; , p.pos(1) = plo[1] +
    //                         (iv_dest[1] + constants::HALF) * dx[1]; ,
    //                         p.pos(2) = plo[2] + (iv_dest[2] +
    //                         constants::HALF) * dx[2];)

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
    //                         p.pos(0) = plo[0] + (iv_dest2[0] +
    //                         constants::HALF) * dx[0]; , p.pos(1) = plo[1] +
    //                         (iv_dest2[1] + constants::HALF) * dx[1]; ,
    //                         p.pos(2) = plo[2] + (iv_dest2[2] +
    //                         constants::HALF)
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
    //                         p.pos(0) = plo[0] + (iv_dest[0] +
    //                         constants::HALF) * dx[0]; , p.pos(1) = plo[1] +
    //                         (iv_dest[1] + constants::HALF) * dx[1]; ,
    //                         p.pos(2) = plo[2] + (iv_dest[2] +
    //                         constants::HALF) * dx[2];)

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
    //                         p.pos(0) = plo[0] + (iv_dest2[0] +
    //                         constants::HALF) * dx[0]; , p.pos(1) = plo[1] +
    //                         (iv_dest2[1] + constants::HALF) * dx[1]; ,
    //                         p.pos(2) = plo[2] + (iv_dest2[2] +
    //                         constants::HALF)
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
    const int msg_per_cell = 1;
    AMREX_ALWAYS_ASSERT(np_per_cell > msg_per_cell);

    amrex::iMultiFab num_particles(
        ParticleBoxArray(LEV), ParticleDistributionMap(LEV), 1, 0,
        amrex::MFInfo());
    amrex::iMultiFab init_offsets(
        ParticleBoxArray(LEV), ParticleDistributionMap(LEV), 1, 0,
        amrex::MFInfo());
    num_particles.setVal(np_per_cell);
    init_offsets.setVal(0);

    for (amrex::MFIter mfi = MakeMFIter(LEV); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();

        const auto ncells = static_cast<int>(box.numPts());
        const int* in = num_particles[mfi].dataPtr();
        int* out = init_offsets[mfi].dataPtr();
        const auto np = amrex::Scan::PrefixSum<int>(
            ncells, [=] AMREX_GPU_DEVICE(int i) -> int { return in[i]; },
            [=] AMREX_GPU_DEVICE(int i, int const& xi) { out[i] = xi; },
            amrex::Scan::Type::exclusive, amrex::Scan::retSum);

        const amrex::Long pid = ParticleType::NextID();
        ParticleType::NextID(pid + np);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            static_cast<amrex::Long>(pid + np) < amrex::LastParticleID,
            "Error: overflow on particle id numbers!");

        const auto my_proc = amrex::ParallelDescriptor::MyProc();
        const auto& offset_arr = init_offsets[mfi].const_array();
        const auto& num_particles_arr = num_particles[mfi].const_array();
        auto& pti = DefineAndReturnParticleTile(LEV, mfi);
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

                    AMREX_D_TERM(
                        p.pos(0) = plo[0] + (iv[0] + constants::HALF) * dx[0];
                        , p.pos(1) = plo[1] + (iv[1] + constants::HALF) * dx[1];
                        ,
                        p.pos(2) = plo[2] + (iv[2] + constants::HALF) * dx[2];)

                    AMREX_D_TERM(p.idata(IntData::i) = iv[0];
                                 , p.idata(IntData::j) = iv[1];
                                 , p.idata(IntData::k) = iv[2];)
                }

                for (int n = start; n < start + msg_per_cell; n++) {
                    auto& pmsg = aos[n];
                    const amrex::Real ts =
                        random_exponential(1.0, engine) + lookahead;

                    pmsg.rdata(RealData::timestamp) = ts;
                    pmsg.rdata(RealData::creation_time) = 0;
                    pmsg.idata(IntData::type_id) = MessageTypes::MESSAGE;
                    pmsg.idata(IntData::pair) = -1;
                }
            });
    }
    Redistribute();

    // Sanity check all initial particles
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MyParIter pti(*this, LEV); pti.isValid(); ++pti) {
        auto index = std::make_pair(pti.index(), pti.LocalTileIndex());

        auto& particle_tile = GetParticles(LEV)[index];
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
            AMREX_ALWAYS_ASSERT(p.id() >= 0);
        });
    }
}

void MessageParticleContainer::sort_messages()
{
    // Taking inspiration from AMReX's SortParticlesByBin
    BL_PROFILE("spades::MessageParticleContainer::sort_messages()");

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = MakeMFIter(LEV); mfi.isValid(); ++mfi) {
        auto& particle_tile = ParticlesAt(LEV, mfi);
        const size_t np = particle_tile.numParticles();

        if (np == 0) {
            continue;
        }

        // BL_PROFILE_VAR(
        //     "spades::MessageParticleContainer::sort_messages::sort_prep",
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
        //     "spades::MessageParticleContainer::sort_messages::sort",
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
            [=](const amrex::Long xi, const amrex::Long yi) {
                const auto& p1 = pstruct[xi];
                const auto& p2 = pstruct[yi];
                return Compare()(p1, p2);
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
                return Compare()(p1, p2);
            });
#endif
        // amrex::Gpu::Device::synchronize();
        // BL_PROFILE_VAR_STOP(sort);

        // Reorder the particles in memory
        // BL_PROFILE_VAR(
        //     "spades::MessageParticleContainer::sort_messages::"
        //     "ReorderParticles",
        //     reorder);
        ReorderParticles(LEV, mfi, cell_list.data());
        // amrex::Gpu::Device::synchronize();
        // BL_PROFILE_VAR_STOP(reorder);
    }
    update_counts();
}

void MessageParticleContainer::update_undefined()
{
    BL_PROFILE("spades::MessageParticleContainer::update_undefined()");

    const auto& plo = Geom(LEV).ProbLoArray();
    const auto& dx = Geom(LEV).CellSizeArray();
    const auto& dom = Geom(LEV).Domain();
    const int lower_count = m_lower_undefined_count;
    const int upper_count = m_upper_undefined_count;
    const int reset_count = m_reset_undefined_count;
    AMREX_ALWAYS_ASSERT(reset_count < upper_count);
    AMREX_ALWAYS_ASSERT(lower_count < reset_count);
    int n_removals = 0;

    MessageParticleContainer pc_adds(
        m_gdb->Geom(), m_gdb->DistributionMap(), m_gdb->boxArray(), ngrow());

    for (amrex::MFIter mfi = MakeMFIter(LEV); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        const auto& cnt_arr = m_message_counts.const_array(mfi);
        const auto& offsets_arr = m_offsets.const_array(mfi);
        auto& particle_tile = GetParticles(LEV)[std::make_pair(gid, tid)];
        auto& particles = particle_tile.GetArrayOfStructs();
        auto* pstruct = particles().dataPtr();

        // remove particles
        // BL_PROFILE_VAR(
        //     "spades::MessageParticleContainer::update_undefined::remove",
        //     remove);
        const auto ncells = static_cast<int>(box.numPts());
        amrex::Gpu::DeviceVector<int> removals(ncells, 0);
        auto* p_removals = removals.data();
        amrex::ParallelFor(ncells, [=] AMREX_GPU_DEVICE(long icell) noexcept {
            const auto iv = box.atOffset(icell);
            const int current_count = cnt_arr(iv, MessageTypes::UNDEFINED);
            if (current_count > upper_count) {
                p_removals[icell] = current_count - reset_count;
            }
        });

        n_removals += amrex::Reduce::Sum(removals.size(), p_removals);

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
        //     "spades::MessageParticleContainer::update_undefined::compute_add",
        //     compute_add);
        amrex::Gpu::DeviceVector<int> additions(ncells, 0);
        auto* p_additions = additions.data();
        amrex::ParallelFor(ncells, [=] AMREX_GPU_DEVICE(long icell) noexcept {
            const auto iv = box.atOffset(icell);
            const int current_count = cnt_arr(iv, MessageTypes::UNDEFINED);
            if (lower_count > current_count) {
                p_additions[icell] = reset_count - current_count;
            }
        });

        amrex::Gpu::DeviceVector<int> update_offsets(ncells, 0);
        auto* p_update_offsets = update_offsets.data();
        const auto np = amrex::Scan::PrefixSum<int>(
            ncells,
            [=] AMREX_GPU_DEVICE(int i) -> int { return p_additions[i]; },
            [=] AMREX_GPU_DEVICE(int i, int const& xi) {
                p_update_offsets[i] = xi;
            },
            amrex::Scan::Type::exclusive, amrex::Scan::retSum);

        const amrex::Long pid = ParticleType::NextID();
        ParticleType::NextID(pid + np);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            static_cast<amrex::Long>(pid + np) < amrex::LastParticleID,
            "Error: overflow on particle id numbers!");

        auto& ptile_adds = pc_adds.DefineAndReturnParticleTile(LEV, mfi);
        ptile_adds.resize(np);
        const auto my_proc = amrex::ParallelDescriptor::MyProc();
        auto* aos = ptile_adds.GetArrayOfStructs().dataPtr();
        amrex::ParallelFor(ncells, [=] AMREX_GPU_DEVICE(long icell) noexcept {
            const int start = p_update_offsets[icell];
            const auto iv = box.atOffset(icell);
            for (int n = start; n < start + p_additions[icell]; n++) {
                auto& p = aos[n];
                p.id() = pid + n;
                p.cpu() = my_proc;

                MarkUndefined()(p);
                p.idata(IntData::sender) = static_cast<int>(dom.index(iv));
                p.idata(IntData::receiver) = static_cast<int>(dom.index(iv));

                AMREX_D_TERM(
                    p.pos(0) = plo[0] + (iv[0] + constants::HALF) * dx[0];
                    , p.pos(1) = plo[1] + (iv[1] + constants::HALF) * dx[1];
                    , p.pos(2) = plo[2] + (iv[2] + constants::HALF) * dx[2];)

                AMREX_D_TERM(p.idata(IntData::i) = iv[0];
                             , p.idata(IntData::j) = iv[1];
                             , p.idata(IntData::k) = iv[2];)
            }
        });
        // amrex::Gpu::Device::synchronize();
        // BL_PROFILE_VAR_STOP(compute_add);
    }

    addParticles(pc_adds, true);

    amrex::ParallelAllReduce::Sum<int>(
        n_removals, amrex::ParallelContext::CommunicatorSub());
    if (n_removals > 0) {
        Redistribute(); // kind of a bummer
    }
}

void MessageParticleContainer::resolve_pairs()
{
    BL_PROFILE("spades::MessageParticleContainer::resolve_pairs()");

    const auto& dom = Geom(LEV).Domain();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = MakeMFIter(LEV); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        const auto& cnt_arr = m_message_counts.const_array(mfi);
        const auto& offsets_arr = m_offsets.const_array(mfi);
        auto& particle_tile = GetParticles(LEV)[std::make_pair(gid, tid)];
        auto& particles = particle_tile.GetArrayOfStructs();
        auto* pstruct = particles().dataPtr();

        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(
                     int i, int j, int AMREX_D_PICK(, , k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const auto getter = Get(iv, cnt_arr, offsets_arr, pstruct);

                for (int n = 0; n < cnt_arr(iv, MessageTypes::ANTI); n++) {

                    auto& pant = getter(n, MessageTypes::ANTI);
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
                            AMREX_ALWAYS_ASSERT(
                                std::abs(
                                    pmsg.rdata(RealData::creation_time) -
                                    pant.rdata(RealData::creation_time)) <
                                constants::EPS);
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

    sort_messages();
}

void MessageParticleContainer::garbage_collect(const amrex::Real gvt)
{
    BL_PROFILE("spades::MessageParticleContainer::garbage_collect()");

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MyParIter pti(*this, LEV); pti.isValid(); ++pti) {
        auto index = std::make_pair(pti.index(), pti.LocalTileIndex());

        auto& particle_tile = GetParticles(LEV)[index];
        const size_t np = particle_tile.numParticles();
        auto& particles = particle_tile.GetArrayOfStructs();
        auto* pstruct = particles().dataPtr();

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pindex) noexcept {
            auto& p = pstruct[pindex];
            if (((p.rdata(RealData::timestamp) < gvt) &&
                 (p.idata(IntData::type_id) != MessageTypes::UNDEFINED)) ||
                ((p.idata(IntData::type_id) == MessageTypes::CONJUGATE) &&
                 (p.rdata(RealData::creation_time) < gvt))) {
                AMREX_ALWAYS_ASSERT(
                    p.idata(IntData::type_id) != MessageTypes::MESSAGE);
                MarkUndefined()(p);
            }
        });
    }
}

amrex::Real MessageParticleContainer::compute_gvt()
{
    BL_PROFILE("spades::MessageParticleContainer::compute_gvt()");
    // If this becomes a performance bottleneck it could be sped up by
    // making a vector of just the message time stamps before the min op

    amrex::Real gvt = constants::LARGE_NUM;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion()) reduction(min : gvt)
#endif
    for (MyParIter pti(*this, LEV); pti.isValid(); ++pti) {
        auto index = std::make_pair(pti.index(), pti.LocalTileIndex());

        const auto& particle_tile = GetParticles(LEV)[index];
        const size_t np = particle_tile.numParticles();
        const auto& particles = particle_tile.GetArrayOfStructs();
        const auto* pstruct = particles().dataPtr();

        amrex::Gpu::DeviceVector<amrex::Real> ts(np, constants::LARGE_NUM);
        auto* p_ts = ts.data();
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pindex) noexcept {
            const auto& p = pstruct[pindex];
            if (p.idata(IntData::type_id) == MessageTypes::MESSAGE) {
                p_ts[pindex] = p.rdata(RealData::timestamp);
            }
        });
        gvt = std::min(amrex::Reduce::Min(np, ts.data()), gvt);
    }
    amrex::ParallelDescriptor::ReduceRealMin(gvt);

    return gvt;
}

void MessageParticleContainer::reposition_messages()
{
    BL_PROFILE("spades::MessageParticleContainer::reposition_messages()");

    const auto& plo = Geom(LEV).ProbLoArray();
    const auto& dx = Geom(LEV).CellSizeArray();
    const auto& dxi = Geom(LEV).InvCellSizeArray();
    const auto& dom = Geom(LEV).Domain();
    const int nbins = 500;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = MakeMFIter(LEV); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        const auto& cnt_arr = m_message_counts.const_array(mfi);
        const auto& offsets_arr = m_offsets.const_array(mfi);
        const auto index = std::make_pair(gid, tid);
        auto& pti = GetParticles(LEV)[index];
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

                        AMREX_D_TERM(
                            p.pos(0) =
                                plo[0] + iv[0] * dx[0] +
                                (typ + 1) * dx[0] / (MessageTypes::NTYPES + 1);
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

void MessageParticleContainer::write_plot_file(const std::string& plt_filename)
{
    BL_PROFILE("spades::MessageParticleContainer::write_plot_file()");
    reposition_messages();
    WritePlotFile(
        plt_filename, identifier(), m_writeflags_real, m_writeflags_int,
        m_real_data_names, m_int_data_names);
}
} // namespace spades::particles
