#include <utility>
#include "MessageParticleContainer.H"

namespace spades::particles {

MessageParticleContainer::MessageParticleContainer(
    amrex::AmrParGDB* par_gdb, int ngrow)
    : SpadesParticleContainer<
          MessageTypes,
          MessageRealData::ncomps,
          MessageIntData::ncomps>(par_gdb, ngrow)
{
    initialize_variable_names();
}

MessageParticleContainer::MessageParticleContainer(
    const amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<amrex::DistributionMapping>& dmap,
    const amrex::Vector<amrex::BoxArray>& ba,
    int ngrow)
    : SpadesParticleContainer<
          MessageTypes,
          MessageRealData::ncomps,
          MessageIntData::ncomps>(geom, dmap, ba, ngrow)
{
    initialize_variable_names();
}

void MessageParticleContainer::initialize_variable_names()
{
    BL_PROFILE("spades::MessageParticleContainer::initialize_variable_names()");

    m_real_data_names.resize(MessageRealData::ncomps, "");
    m_writeflags_real.resize(MessageRealData::ncomps, 0);
    m_int_data_names.resize(MessageIntData::ncomps, "");
    m_writeflags_int.resize(MessageIntData::ncomps, 0);

    m_real_data_names[MessageRealData::timestamp] = "timestamp";
    m_writeflags_real[MessageRealData::timestamp] = 1;
    m_real_data_names[MessageRealData::old_timestamp] = "old_timestamp";
    m_writeflags_real[MessageRealData::old_timestamp] = 1;
    m_real_data_names[MessageRealData::creation_time] = "creation_time";
    m_writeflags_real[MessageRealData::creation_time] = 1;

    m_int_data_names[MessageIntData::type_id] = "type_id";
    m_writeflags_int[MessageIntData::type_id] = 1;
    m_int_data_names[MessageIntData::sender] = "sender";
    m_writeflags_int[MessageIntData::sender] = 1;
    m_int_data_names[MessageIntData::receiver] = "receiver";
    m_writeflags_int[MessageIntData::receiver] = 1;
    m_int_data_names[MessageIntData::receiver] = "pair";
    m_writeflags_int[MessageIntData::receiver] = 0;
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

    //                     p.idata(MessageIntData::type_id) =
    //                         MessageTypes::anti_message;
    //                     p.idata(MessageIntData::sender) =
    //                     dom.index(iv_src);
    //                     p.idata(MessageIntData::receiver) =
    //                     dom.index(iv_dest);
    //                     p.rdata(MessageRealData::timestamp) =
    //                         random_exponential(1.0, engine) + lookahead + 20;

    //                     AMREX_D_TERM(
    //                         p.pos(0) = plo[0] + (iv_dest[0] +
    //                         constants::HALF) * dx[0]; , p.pos(1) = plo[1] +
    //                         (iv_dest[1] + constants::HALF) * dx[1]; ,
    //                         p.pos(2) = plo[2] + (iv_dest[2] +
    //                         constants::HALF) * dx[2];)

    //                     AMREX_D_TERM(p.idata(MessageIntData::i) =
    //                     iv_dest[0];
    //                                  , p.idata(MessageIntData::j) =
    //                                  iv_dest[1]; ,
    //                                  p.idata(MessageIntData::k) =
    //                                  iv_dest[2];)

    //                     pti.push_back(p);
    //                 }
    //                 {
    //                     ParticleType
    //                     p; p.id() =
    //                         ParticleType::NextID();
    //                     p.cpu() = amrex::ParallelDescriptor::MyProc();

    //                     p.idata(MessageIntData::type_id) =
    //                         MessageTypes::message;
    //                     p.idata(MessageIntData::sender) =
    //                     dom.index(iv_src);
    //                     p.idata(MessageIntData::receiver) =
    //                     dom.index(iv_dest2);
    //                     p.rdata(MessageRealData::timestamp) =
    //                         random_exponential(1.0, engine) + lookahead;

    //                     AMREX_D_TERM(
    //                         p.pos(0) = plo[0] + (iv_dest2[0] +
    //                         constants::HALF) * dx[0]; , p.pos(1) = plo[1] +
    //                         (iv_dest2[1] + constants::HALF) * dx[1]; ,
    //                         p.pos(2) = plo[2] + (iv_dest2[2] +
    //                         constants::HALF)
    //                         * dx[2];)

    //                     AMREX_D_TERM(p.idata(MessageIntData::i) =
    //                     iv_dest2[0];
    //                                  , p.idata(MessageIntData::j) =
    //                                  iv_dest2[1];
    //                                  ,
    //                                  p.idata(MessageIntData::k) =
    //                                  iv_dest2[2];)

    //                     pti.push_back(p);
    //                 }
    //                 { // creating another particle
    //                     ParticleType
    //                     p; p.id() =
    //                         ParticleType::NextID();
    //                     p.cpu() = amrex::ParallelDescriptor::MyProc();

    //                     p.idata(MessageIntData::type_id) =
    //                         MessageTypes::message;
    //                     p.idata(MessageIntData::sender) =
    //                     dom.index(iv_src);
    //                     p.idata(MessageIntData::receiver) =
    //                     dom.index(iv_dest);
    //                     p.rdata(MessageRealData::timestamp) =
    //                         random_exponential(1.0, engine) + lookahead;

    //                     AMREX_D_TERM(
    //                         p.pos(0) = plo[0] + (iv_dest[0] +
    //                         constants::HALF) * dx[0]; , p.pos(1) = plo[1] +
    //                         (iv_dest[1] + constants::HALF) * dx[1]; ,
    //                         p.pos(2) = plo[2] + (iv_dest[2] +
    //                         constants::HALF) * dx[2];)

    //                     AMREX_D_TERM(p.idata(MessageIntData::i) =
    //                     iv_dest[0];
    //                                  , p.idata(MessageIntData::j) =
    //                                  iv_dest[1]; ,
    //                                  p.idata(MessageIntData::k) =
    //                                  iv_dest[2];)

    //                     pti.push_back(p);
    //                 }
    //                 { // creating another particle
    //                     ParticleType
    //                     p; p.id() =
    //                         ParticleType::NextID();
    //                     p.cpu() = amrex::ParallelDescriptor::MyProc();

    //                     p.idata(MessageIntData::type_id) =
    //                         MessageTypes::undefined;
    //                     p.idata(MessageIntData::sender) =
    //                     dom.index(iv_src);
    //                     p.idata(MessageIntData::receiver) =
    //                     dom.index(iv_dest2);
    //                     p.rdata(MessageRealData::timestamp) =
    //                         random_exponential(1.0, engine) + lookahead;

    //                     AMREX_D_TERM(
    //                         p.pos(0) = plo[0] + (iv_dest2[0] +
    //                         constants::HALF) * dx[0]; , p.pos(1) = plo[1] +
    //                         (iv_dest2[1] + constants::HALF) * dx[1]; ,
    //                         p.pos(2) = plo[2] + (iv_dest2[2] +
    //                         constants::HALF)
    //                         * dx[2];)

    //                     AMREX_D_TERM(p.idata(MessageIntData::i) =
    //                     iv_dest2[0];
    //                                  , p.idata(MessageIntData::j) =
    //                                  iv_dest2[1];
    //                                  ,
    //                                  p.idata(MessageIntData::k) =
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

                    MarkMessageUndefined()(p);
                    p.idata(MessageIntData::sender) =
                        static_cast<int>(dom.index(iv));
                    p.idata(MessageIntData::receiver) =
                        static_cast<int>(dom.index(iv));

                    AMREX_D_TERM(
                        p.pos(0) = plo[0] + (iv[0] + constants::HALF) * dx[0];
                        , p.pos(1) = plo[1] + (iv[1] + constants::HALF) * dx[1];
                        ,
                        p.pos(2) = plo[2] + (iv[2] + constants::HALF) * dx[2];)

                    AMREX_D_TERM(p.idata(MessageIntData::i) = iv[0];
                                 , p.idata(MessageIntData::j) = iv[1];
                                 , p.idata(MessageIntData::k) = iv[2];)
                }

                for (int n = start; n < start + msg_per_cell; n++) {
                    auto& pmsg = aos[n];
                    const amrex::Real ts =
                        random_exponential(1.0, engine) + lookahead;

                    pmsg.rdata(MessageRealData::timestamp) = ts;
                    pmsg.rdata(MessageRealData::creation_time) = 0;
                    pmsg.idata(MessageIntData::type_id) = MessageTypes::MESSAGE;
                    pmsg.idata(MessageIntData::pair) = -1;
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
                valid_type = p.idata(MessageIntData::type_id) == typ;
                if (valid_type) {
                    break;
                }
            }
            AMREX_ALWAYS_ASSERT(valid_type);
            AMREX_ALWAYS_ASSERT(p.id() >= 0);
        });
    }
}
void MessageParticleContainer::sort()
{
    BL_PROFILE("spades::MessageParticleContainer::sort()");

    sort_impl(CompareMessage());
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
        const auto& cnt_arr = m_counts.const_array(mfi);
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

                MarkMessageUndefined()(p);
                p.idata(MessageIntData::sender) =
                    static_cast<int>(dom.index(iv));
                p.idata(MessageIntData::receiver) =
                    static_cast<int>(dom.index(iv));

                AMREX_D_TERM(
                    p.pos(0) = plo[0] + (iv[0] + constants::HALF) * dx[0];
                    , p.pos(1) = plo[1] + (iv[1] + constants::HALF) * dx[1];
                    , p.pos(2) = plo[2] + (iv[2] + constants::HALF) * dx[2];)

                AMREX_D_TERM(p.idata(MessageIntData::i) = iv[0];
                             , p.idata(MessageIntData::j) = iv[1];
                             , p.idata(MessageIntData::k) = iv[2];)
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
        const auto& cnt_arr = m_counts.const_array(mfi);
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
                    AMREX_ALWAYS_ASSERT(pant.idata(MessageIntData::pair) != -1);
                    AMREX_ALWAYS_ASSERT(
                        pant.idata(MessageIntData::receiver) == dom.index(iv));

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
                        if ((pmsg.idata(MessageIntData::pair) ==
                             pant.idata(MessageIntData::pair)) &&
                            (std::abs(
                                 pmsg.rdata(MessageRealData::timestamp) -
                                 pant.rdata(MessageRealData::timestamp)) <
                             constants::EPS)) {
                            AMREX_ALWAYS_ASSERT(
                                pmsg.idata(MessageIntData::sender) ==
                                pant.idata(MessageIntData::sender));
                            AMREX_ALWAYS_ASSERT(
                                pmsg.idata(MessageIntData::receiver) ==
                                pant.idata(MessageIntData::receiver));
                            AMREX_ALWAYS_ASSERT(
                                std::abs(
                                    pmsg.rdata(MessageRealData::creation_time) -
                                    pant.rdata(
                                        MessageRealData::creation_time)) <
                                constants::EPS);
                            MarkMessageUndefined()(pant);
                            MarkMessageUndefined()(pmsg);
                            found_pair = true;
                            break;
                        }
                    }
                    AMREX_ALWAYS_ASSERT(found_pair);
                }
            });
    }

    sort();
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
            if (((p.rdata(MessageRealData::timestamp) < gvt) &&
                 (p.idata(MessageIntData::type_id) !=
                  MessageTypes::UNDEFINED)) ||
                ((p.idata(MessageIntData::type_id) ==
                  MessageTypes::CONJUGATE) &&
                 (p.rdata(MessageRealData::creation_time) < gvt))) {
                AMREX_ALWAYS_ASSERT(
                    p.idata(MessageIntData::type_id) != MessageTypes::MESSAGE);
                MarkMessageUndefined()(p);
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
            if (p.idata(MessageIntData::type_id) == MessageTypes::MESSAGE) {
                p_ts[pindex] = p.rdata(MessageRealData::timestamp);
            }
        });
        gvt = std::min(amrex::Reduce::Min(np, ts.data()), gvt);
    }
    amrex::ParallelDescriptor::ReduceRealMin(gvt);

    return gvt;
}

} // namespace spades::particles
