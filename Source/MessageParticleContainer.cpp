#include <utility>
#include "MessageParticleContainer.H"

namespace spades::particles {

MessageParticleContainer::MessageParticleContainer(
    amrex::AmrParGDB* par_gdb, int ngrow)
    : SpadesParticleContainer<
          MessageTypes,
          MessageRealData::ncomps,
          MessageIntData::ncomps>(par_gdb, ngrow)
{}

MessageParticleContainer::MessageParticleContainer(
    const amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<amrex::DistributionMapping>& dmap,
    const amrex::Vector<amrex::BoxArray>& ba,
    int ngrow)
    : SpadesParticleContainer<
          MessageTypes,
          MessageRealData::ncomps,
          MessageIntData::ncomps>(geom, dmap, ba, ngrow)
{}

void MessageParticleContainer::read_parameters()
{
    SpadesParticleContainer::read_parameters();
    {
        amrex::ParmParse pp("spades");
        pp.query("messages_per_lp", m_messages_per_lp);
    }
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

    AMREX_D_TERM(m_int_data_names[MessageIntData::i] = "i";
                 m_writeflags_int[MessageIntData::i] = 1;
                 , m_int_data_names[MessageIntData::j] = "j";
                 m_writeflags_int[MessageIntData::j] = 1;
                 , m_int_data_names[MessageIntData::k] = "k";
                 m_writeflags_int[MessageIntData::k] = 1;)
    m_int_data_names[MessageIntData::type_id] = "type_id";
    m_writeflags_int[MessageIntData::type_id] = 1;
    m_int_data_names[MessageIntData::sender_lp] = "sender_lp";
    m_writeflags_int[MessageIntData::sender_lp] = 1;
    m_int_data_names[MessageIntData::sender_entity] = "sender_entity";
    m_writeflags_int[MessageIntData::sender_entity] = 1;
    m_int_data_names[MessageIntData::receiver_lp] = "receiver_lp";
    m_writeflags_int[MessageIntData::receiver_lp] = 1;
    m_int_data_names[MessageIntData::receiver_entity] = "receiver_entity";
    m_writeflags_int[MessageIntData::receiver_entity] = 1;
    m_int_data_names[MessageIntData::pair] = "pair";
    m_writeflags_int[MessageIntData::pair] = 0;
}

void MessageParticleContainer::initialize_messages(const amrex::Real lookahead)
{
    BL_PROFILE("spades::MessageParticleContainer::initialize_messages()");

    const auto& plo = Geom(LEV).ProbLoArray();
    const auto& dx = Geom(LEV).CellSizeArray();
    const auto& dom = Geom(LEV).Domain();
    const auto messages_per_lp = m_messages_per_lp;
    const int total_messages_per_lp = 3 * messages_per_lp;
    AMREX_ALWAYS_ASSERT(total_messages_per_lp > messages_per_lp);

    amrex::iMultiFab num_particles(
        ParticleBoxArray(LEV), ParticleDistributionMap(LEV), 1, 0,
        amrex::MFInfo());
    amrex::iMultiFab init_offsets(
        ParticleBoxArray(LEV), ParticleDistributionMap(LEV), 1, 0,
        amrex::MFInfo());
    num_particles.setVal(total_messages_per_lp);
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
                    p.idata(MessageIntData::sender_lp) =
                        static_cast<int>(dom.index(iv));
                    p.idata(MessageIntData::sender_entity) = 0;
                    p.idata(MessageIntData::receiver_lp) =
                        static_cast<int>(dom.index(iv));
                    p.idata(MessageIntData::receiver_entity) = 0;

                    AMREX_D_TERM(
                        p.pos(0) = plo[0] + (iv[0] + constants::HALF) * dx[0];
                        , p.pos(1) = plo[1] + (iv[1] + constants::HALF) * dx[1];
                        ,
                        p.pos(2) = plo[2] + (iv[2] + constants::HALF) * dx[2];)

                    AMREX_D_TERM(p.idata(MessageIntData::i) = iv[0];
                                 , p.idata(MessageIntData::j) = iv[1];
                                 , p.idata(MessageIntData::k) = iv[2];)
                }

                for (int n = start; n < start + messages_per_lp; n++) {
                    auto& pmsg = aos[n];
                    const amrex::Real ts =
                        random_exponential(1.0, engine) + lookahead;

                    pmsg.rdata(MessageRealData::timestamp) = ts;
                    pmsg.rdata(MessageRealData::creation_time) = 0;
                    pmsg.idata(MessageIntData::type_id) = MessageTypes::MESSAGE;
                    pmsg.idata(MessageIntData::pair) = -1;
                }
            });

        // This is necessary
        amrex::Gpu::streamSynchronize();
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
            AMREX_ASSERT(valid_type);
            AMREX_ASSERT(p.id() >= 0);
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
            const int msg_count = cnt_arr(iv, MessageTypes::MESSAGE);
            const int target_count = std::max(2 * msg_count, 2);
            const int current_count = cnt_arr(iv, MessageTypes::UNDEFINED);
            if (current_count > target_count) {
                p_removals[icell] = current_count - target_count;
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
        // amrex::Gpu::streamSynchronize();
        // BL_PROFILE_VAR_STOP(remove);

        // BL_PROFILE_VAR(
        //     "spades::MessageParticleContainer::update_undefined::compute_add",
        //     compute_add);
        amrex::Gpu::DeviceVector<int> additions(ncells, 0);
        auto* p_additions = additions.data();
        amrex::ParallelFor(ncells, [=] AMREX_GPU_DEVICE(long icell) noexcept {
            const auto iv = box.atOffset(icell);
            const int msg_count = cnt_arr(iv, MessageTypes::MESSAGE);
            const int target_count = std::max(2 * msg_count, 2);
            const int current_count = cnt_arr(iv, MessageTypes::UNDEFINED);
            if (target_count > current_count) {
                p_additions[icell] = target_count - current_count;
            }
#ifdef AMREX_DEBUG
            if (p_removals[icell] > 0) {
                AMREX_ASSERT(p_additions[icell] == 0);
            }
            if (p_additions[icell] > 0) {
                AMREX_ASSERT(p_removals[icell] == 0);
            }
#endif
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
                p.idata(MessageIntData::sender_lp) =
                    static_cast<int>(dom.index(iv));
                p.idata(MessageIntData::sender_entity) = 0;
                p.idata(MessageIntData::receiver_lp) =
                    static_cast<int>(dom.index(iv));
                p.idata(MessageIntData::receiver_entity) = 0;

                AMREX_D_TERM(
                    p.pos(0) = plo[0] + (iv[0] + constants::HALF) * dx[0];
                    , p.pos(1) = plo[1] + (iv[1] + constants::HALF) * dx[1];
                    , p.pos(2) = plo[2] + (iv[2] + constants::HALF) * dx[2];)

                AMREX_D_TERM(p.idata(MessageIntData::i) = iv[0];
                             , p.idata(MessageIntData::j) = iv[1];
                             , p.idata(MessageIntData::k) = iv[2];)
            }
        });
        // amrex::Gpu::streamSynchronize();
        // BL_PROFILE_VAR_STOP(compute_add);

        // This is necessary
        amrex::Gpu::streamSynchronize();
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

#ifdef AMREX_DEBUG
    const auto& dom = Geom(LEV).Domain();
#endif

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
                    AMREX_ASSERT(pant.idata(MessageIntData::pair) != -1);
                    AMREX_ASSERT(
                        pant.idata(MessageIntData::receiver_lp) ==
                        dom.index(iv));

#ifdef AMREX_DEBUG
                    bool found_pair = false;
#endif
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
                            AMREX_ASSERT(
                                pmsg.idata(MessageIntData::sender_lp) ==
                                pant.idata(MessageIntData::sender_lp));
                            AMREX_ASSERT(
                                pmsg.idata(MessageIntData::sender_entity) ==
                                pant.idata(MessageIntData::sender_entity));
                            AMREX_ASSERT(
                                pmsg.idata(MessageIntData::receiver_lp) ==
                                pant.idata(MessageIntData::receiver_lp));
                            AMREX_ASSERT(
                                pmsg.idata(MessageIntData::receiver_entity) ==
                                pant.idata(MessageIntData::receiver_entity));
                            AMREX_ASSERT(
                                std::abs(
                                    pmsg.rdata(MessageRealData::creation_time) -
                                    pant.rdata(
                                        MessageRealData::creation_time)) <
                                constants::EPS);
                            MarkMessageUndefined()(pant);
                            MarkMessageUndefined()(pmsg);
#ifdef AMREX_DEBUG
                            found_pair = true;
#endif
                            break;
                        }
                    }
                    AMREX_ASSERT(found_pair);
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
                AMREX_ASSERT(
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
