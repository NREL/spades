#include <memory>

#include "SPADES.H"

namespace spades {
SPADES::SPADES()
{
    BL_PROFILE("spades::SPADES::SPADES()");
    read_parameters();
    int nlevs_max = max_level + 1;

    if (max_level > 0) {
        amrex::Abort(
            "spades::SPADES::SPADES(): not supporting multilevel right now");
    }

    m_state_varnames.push_back("lvt");
    m_state_varnames.push_back("rollback");
    AMREX_ALWAYS_ASSERT(m_state_varnames.size() == constants::N_STATES);

    m_message_counts_varnames.resize(particles::MessageTypes::NTYPES);
    m_message_counts_varnames[particles::MessageTypes::MESSAGE] = "message";
    m_message_counts_varnames[particles::MessageTypes::PROCESSED] = "processed";
    m_message_counts_varnames[particles::MessageTypes::ANTI_MESSAGE] =
        "anti_message";
    m_message_counts_varnames[particles::MessageTypes::CONJUGATE] = "conjugate";
    m_message_counts_varnames[particles::MessageTypes::UNDEFINED] = "undefined";
    for (const auto& mm : m_message_counts_varnames) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            !mm.empty(), "Variable name needs to be specified");
    }

    for (const auto& vname : m_state_varnames) {
        m_spades_varnames.push_back(vname);
    }
    for (const auto& vname : m_message_counts_varnames) {
        m_spades_varnames.push_back(vname);
    }
    AMREX_ALWAYS_ASSERT(
        m_spades_varnames.size() ==
        (constants::N_STATES + particles::MessageTypes::NTYPES));

    m_isteps.resize(nlevs_max, 0);
    m_nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= max_level; ++lev) {
        m_nsubsteps[lev] = MaxRefRatio(lev - 1);
    }

    m_gvts.resize(nlevs_max, constants::LOW_NUM);
    m_ts_new.resize(nlevs_max, 0.0);
    m_ts_old.resize(nlevs_max, constants::LOW_NUM);
    m_dts.resize(nlevs_max, constants::LARGE_NUM);
    m_lbts.resize(nlevs_max, constants::LOW_NUM);

    m_state.resize(nlevs_max);
    m_plt_mf.resize(nlevs_max);
}

SPADES::~SPADES() = default;

void SPADES::init_data()
{
    BL_PROFILE("spades::SPADES::init_data()");

    if (m_restart_chkfile.empty()) {

        init_rng();

        // start simulation from the beginning
        const amrex::Real time = 0.0;
        set_ics();

        init_particle_container();

        InitFromScratch(time);

        const int lev = 0;
        update_gvt(lev);
        update_lbts(lev);

        compute_dt();

        if (m_chk_int > 0) {
            write_checkpoint_file();
        }
    } else {
        // restart from a checkpoint
        read_checkpoint_file();
    }

    if (m_plot_int > 0) {
        write_plot_file();
    }

    amrex::Print() << "Grid summary: " << std::endl;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        printGridSummary(amrex::OutStream(), 0, finest_level);
    }
    amrex::Print() << "Particle summary: " << std::endl;
    for (int lev = 0; lev <= finest_level; ++lev) {
        const auto np = m_pc->NumberOfParticlesAtLevel(lev, true);
        amrex::Print() << "  Level " << lev << "   " << np << " particles"
                       << std::endl;
    }
}

void SPADES::init_particle_container()
{
    BL_PROFILE("spades::SPADES::init_particle_container()");
    m_pc =
        std::make_unique<particles::CellSortedParticleContainer>(GetParGDB());
}

void SPADES::read_parameters()
{
    BL_PROFILE("spades::SPADES::read_parameters()");

    {
        amrex::ParmParse pp;
        pp.query("max_step", m_max_step);
        pp.query("stop_time", m_stop_time);
    }

    {
        amrex::ParmParse pp("amr");

        pp.query("regrid_int", m_regrid_int);
        pp.query("plot_file", m_plot_file);
        pp.query("plot_int", m_plot_int);
        pp.query("chk_file", m_chk_file);
        pp.query("chk_int", m_chk_int);
        pp.query("restart", m_restart_chkfile);
        pp.query("file_name_digits", m_file_name_digits);
        pp.query("rng_file_name_digits", m_rng_file_name_digits);
        AMREX_ALWAYS_ASSERT(
            m_rng_file_name_digits >= amrex::ParallelDescriptor::NProcs());
    }

    {
        amrex::ParmParse pp("spades");
        pp.get("ic_type", m_ic_type);
        pp.query("write_particles", m_write_particles);
        pp.query("seed", m_seed);
        pp.query("lookahead", m_lookahead);
        pp.query("window_size", m_window_size);
        pp.query("messages_per_step", m_messages_per_step);
    }

    // force periodic bcs
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
        if (!amrex::DefaultGeometry().isPeriodic(dir)) {
            amrex::Abort(
                "spades::SPADES::read_parameters(): geometry needs to be "
                "periodic");
        }
    }
}

void SPADES::evolve()
{
    BL_PROFILE("spades::SPADES::evolve()");

    amrex::Real cur_time = m_ts_new[0];
    int last_plot_file_step = 0;
    int last_chk_file_step = 0;

    for (int step = m_isteps[0]; step < m_max_step && cur_time < m_stop_time;
         ++step) {
        compute_dt();

        amrex::Print() << "\n==============================================="
                          "================================================="
                       << std::endl;
        amrex::Print() << "Step: " << step << " dt : " << m_dts[0]
                       << " time: " << cur_time << " to " << cur_time + m_dts[0]
                       << std::endl;

        // m_fillpatch_op->fillpatch(0, cur_time, m_f[0]);
        time_step(0, cur_time, 1);

        post_time_step();

        cur_time += m_dts[0];

        // sync up time
        for (int lev = 0; lev <= finest_level; ++lev) {
            m_ts_new[lev] = cur_time;
        }

        if (m_plot_int > 0 && (step + 1) % m_plot_int == 0) {
            last_plot_file_step = step + 1;
            write_plot_file();
        }

        if (m_chk_int > 0 && (step + 1) % m_chk_int == 0) {
            last_chk_file_step = step + 1;
            write_checkpoint_file();
        }

        if (cur_time >= m_stop_time - 1.e-6 * m_dts[0]) {
            break;
        }
    }
    if (m_plot_int > 0 && m_isteps[0] > last_plot_file_step) {
        write_plot_file();
    }
    if (m_chk_int > 0 && m_isteps[0] > last_chk_file_step) {
        write_checkpoint_file();
    }
}

// advance a level by dt
// includes a recursive call for finer levels
void SPADES::time_step(
    const int lev, const amrex::Real time, const int iteration)
{
    BL_PROFILE("spades::SPADES::time_step()");
    if (m_regrid_int > 0) // We may need to regrid
    {

        // help keep track of whether a level was already regridded
        // from a coarser level call to regrid
        static amrex::Vector<int> last_regrid_step(max_level + 1, 0);

        // regrid changes level "lev+1" so we don't regrid on max_level
        // also make sure we don't regrid fine levels again if
        // it was taken care of during a coarser regrid
        if (lev < max_level && m_isteps[lev] > last_regrid_step[lev]) {
            if (m_isteps[lev] % m_regrid_int == 0) {
                // regrid could add newly refine levels (if finest_level <
                // max_level) so we save the previous finest level index
                int old_finest = finest_level;
                regrid(lev, time);

                // mark that we have regridded this level already
                for (int k = lev; k <= finest_level; ++k) {
                    last_regrid_step[k] = m_isteps[k];
                }

                // if there are newly created levels, set the time step
                // dt gets halved here
                for (int k = old_finest + 1; k <= finest_level; ++k) {
                    m_dts[k] = m_dts[k - 1] / MaxRefRatio(k - 1);
                }
                if (amrex::ParallelDescriptor::IOProcessor()) {
                    amrex::Print()
                        << "Grid summary after regrid: " << std::endl;
                    printGridSummary(amrex::OutStream(), 0, finest_level);
                }
                m_pc->sort_particles();
            }
        }
    }

    if (Verbose() != 0) {
        amrex::Print() << "[Level " << lev << " step " << m_isteps[lev] + 1
                       << "] ";
        amrex::Print() << "Advance with time = " << m_ts_new[lev]
                       << " dt = " << m_dts[lev] << " gvt = " << m_gvts[lev]
                       << " lbts = " << m_lbts[lev] << std::endl;
    }

    if (lev < finest_level) {
        // m_fillpatch_op->fillpatch(lev + 1, m_ts_new[lev + 1], m_f[lev + 1]);
        for (int i = 1; i <= m_nsubsteps[lev + 1]; ++i) {
            time_step(lev + 1, time + (i - 1) * m_dts[lev + 1], i);
        }
    }

    advance(lev, time, m_dts[lev], iteration, m_nsubsteps[lev]);

    ++m_isteps[lev];

    if (Verbose() != 0) {
        amrex::Print() << "[Level " << lev << " step " << m_isteps[lev]
                       << "] statistics:" << std::endl;
        const auto n_cells = CountCells(lev);
        amrex::Print() << "  Advanced " << n_cells << " nodes and "
                       << m_pc->TotalNumberOfParticles(lev != 0) << " particles"
                       << std::endl;
        const auto n_messages =
            m_pc->total_count(lev, particles::MessageTypes::MESSAGE);
        amrex::Print() << "  " << n_messages << " messages" << std::endl;
        AMREX_ALWAYS_ASSERT(n_messages == n_cells);
    }
}

void SPADES::advance(
    const int lev,
    const amrex::Real /*time*/,
    const amrex::Real dt_lev,
    const int /*iteration*/,
    const int /*ncycle*/)
{
    BL_PROFILE("spades::SPADES::advance()");

    m_pc->update_undefined();

    update_gvt(lev);
    m_pc->garbage_collect(m_gvts[lev]);
    m_pc->sort_particles();

    update_lbts(lev);

    process_messages(lev);
    m_pc->Redistribute();
    m_pc->sort_particles();

    rollback(lev);

    m_ts_old[lev] = m_ts_new[lev]; // old time is now current time (time)
    m_ts_new[lev] += dt_lev;       // new time is ahead
}

void SPADES::post_time_step()
{
    BL_PROFILE("spades::SPADES::post_time_step()");
}

void SPADES::process_messages(const int lev)
{
    BL_PROFILE("spades::SPADES::process_messages()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_state_ngrow == m_pc->ngrow(),
        "Particle and state cells must be equal for now");

    const auto& plo = Geom(lev).ProbLoArray();
    const auto& dx = Geom(lev).CellSizeArray();
    const auto& dom = Geom(lev).Domain();
    const auto& dlo = dom.smallEnd();
    const auto& dhi = dom.bigEnd();
    const auto lbts = m_lbts[lev];
    const auto lookahead = m_lookahead;
    const auto window_size = m_window_size;
    const auto messages_per_step = m_messages_per_step;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = m_pc->MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        const auto& sarr = m_state[lev].array(mfi);
        const auto& cnt_arr = m_pc->message_counts(lev).const_array(mfi);
        const auto& offsets_arr = m_pc->offsets(lev).const_array(mfi);
        const auto index = std::make_pair(gid, tid);
        auto& pti = m_pc->GetParticles(lev)[index];
        auto& particles = pti.GetArrayOfStructs();
        auto* pstruct = particles().dataPtr();

        amrex::ParallelForRNG(
            box, [=] AMREX_GPU_DEVICE(
                     int i, int j, int AMREX_D_PICK(, , k),
                     amrex::RandomEngine const& engine) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const auto getter =
                    particles::Get(iv, cnt_arr, offsets_arr, pstruct);

                for (int n = 0;
                     n < amrex::min<int>(
                             cnt_arr(iv, particles::MessageTypes::MESSAGE),
                             messages_per_step);
                     n++) {

                    auto& prcv = getter(n, particles::MessageTypes::MESSAGE);
                    if (prcv.rdata(particles::RealData::timestamp) >=
                        lbts + lookahead + window_size) {
                        return;
                    }
                    AMREX_ALWAYS_ASSERT(
                        sarr(iv, constants::LVT_IDX) <
                        prcv.rdata(particles::RealData::timestamp));

                    // process the event
                    AMREX_ALWAYS_ASSERT(
                        dom.atOffset(
                            prcv.idata(particles::IntData::receiver)) == iv);
                    prcv.rdata(particles::RealData::old_timestamp) =
                        sarr(iv, constants::LVT_IDX);
                    sarr(iv, constants::LVT_IDX) =
                        prcv.rdata(particles::RealData::timestamp);
                    prcv.idata(particles::IntData::type_id) =
                        particles::MessageTypes::PROCESSED;

                    // Create a new message to send
                    auto& psnd =
                        getter(2 * n, particles::MessageTypes::UNDEFINED);
                    amrex::IntVect iv_dest(AMREX_D_DECL(
                        amrex::Random_int(dhi[0] - dlo[0] + 1, engine) + dlo[0],
                        amrex::Random_int(dhi[1] - dlo[1] + 1, engine) + dlo[1],
                        amrex::Random_int(dhi[2] - dlo[2] + 1, engine) +
                            dlo[2]));
                    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> pos = {
                        AMREX_D_DECL(
                            plo[0] + (iv_dest[0] + 0.5) * dx[0],
                            plo[1] + (iv_dest[1] + 0.5) * dx[1],
                            plo[2] + (iv_dest[2] + 0.5) * dx[2])};
                    const amrex::Real ts = sarr(iv, constants::LVT_IDX) +
                                           random_exponential(1.0, engine) +
                                           lookahead;
                    // FIXME, in general, Create is clunky. Better way?
                    particles::Create()(
                        psnd, ts, pos, iv_dest, static_cast<int>(dom.index(iv)),
                        static_cast<int>(dom.index(iv_dest)));
                    const auto pair = static_cast<int>(
                        pairing_function(prcv.cpu(), prcv.id()));
                    psnd.idata(particles::IntData::pair) = pair;
                    psnd.rdata(particles::RealData::creation_time) =
                        sarr(iv, constants::LVT_IDX);
                    // Create the conjugate message
                    auto& pcnj =
                        getter(2 * n + 1, particles::MessageTypes::UNDEFINED);

                    // FIXME, could do a copy. Or just pass p.pos to Create
                    // This is weird. The conjugate
                    // position is iv but the receiver is
                    // still updated (we need to know who to
                    // send this to)
                    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> conj_pos = {
                        AMREX_D_DECL(
                            plo[0] + (iv[0] + 0.5) * dx[0],
                            plo[1] + (iv[1] + 0.5) * dx[1],
                            plo[2] + (iv[2] + 0.5) * dx[2])};

                    particles::Create()(
                        pcnj, ts, conj_pos, iv, static_cast<int>(dom.index(iv)),
                        static_cast<int>(dom.index(iv_dest)));
                    pcnj.idata(particles::IntData::pair) = pair;
                    pcnj.rdata(particles::RealData::creation_time) =
                        sarr(iv, constants::LVT_IDX);
                    pcnj.idata(particles::IntData::type_id) =
                        particles::MessageTypes::CONJUGATE;
                }
            });
    }
}

void SPADES::rollback(const int lev)
{
    BL_PROFILE("spades::SPADES::rollback()");

    if (m_window_size <= 0) {
        return;
    }

    const auto& plo = Geom(lev).ProbLoArray();
    const auto& dx = Geom(lev).CellSizeArray();
    const auto& dom = Geom(lev).Domain();

    amrex::iMultiFab rollback(boxArray(lev), DistributionMap(lev), 1, 0);
    rollback.setVal(1.0);
    amrex::Long require_rollback = rollback.sum(0);

    const int max_iter = 100;
    int iter = 0;
    m_state[lev].setVal(0.0, constants::RLB_IDX, 1);
    while ((require_rollback > 0) && (iter < max_iter)) {

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi = m_pc->MakeMFIter(lev); mfi.isValid(); ++mfi) {
            const amrex::Box& box = mfi.tilebox();
            const int gid = mfi.index();
            const int tid = mfi.LocalTileIndex();
            const auto& sarr = m_state[lev].array(mfi);
            const auto& rarr = rollback.array(mfi);
            const auto& cnt_arr = m_pc->message_counts(lev).const_array(mfi);
            const auto& offsets_arr = m_pc->offsets(lev).const_array(mfi);
            const auto index = std::make_pair(gid, tid);
            auto& pti = m_pc->GetParticles(lev)[index];
            auto& particles = pti.GetArrayOfStructs();
            auto* pstruct = particles().dataPtr();

            amrex::ParallelFor(
                box, [=] AMREX_GPU_DEVICE(
                         int i, int j, int AMREX_D_PICK(, , k)) noexcept {
                    const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                    const auto getter =
                        particles::Get(iv, cnt_arr, offsets_arr, pstruct);

                    const amrex::Real msg_lvt =
                        cnt_arr(iv, particles::MessageTypes::MESSAGE) > 0
                            ? getter(0, particles::MessageTypes::MESSAGE)
                                  .rdata(particles::RealData::timestamp)
                            : constants::LARGE_NUM;
                    const amrex::Real ant_lvt =
                        cnt_arr(iv, particles::MessageTypes::ANTI_MESSAGE) > 0
                            ? getter(0, particles::MessageTypes::ANTI_MESSAGE)
                                  .rdata(particles::RealData::timestamp)
                            : constants::LARGE_NUM;
                    const amrex::Real rollback_timestamp =
                        amrex::min<amrex::Real>(msg_lvt, ant_lvt);

                    rarr(iv) = 0;
                    if (rollback_timestamp <= sarr(iv, constants::LVT_IDX)) {
                        AMREX_ALWAYS_ASSERT(
                            cnt_arr(iv, particles::MessageTypes::PROCESSED) >
                            0);
                        rarr(iv) = 1;
                        sarr(iv, constants::RLB_IDX) += 1;

                        for (int n = 0;
                             n <
                             cnt_arr(iv, particles::MessageTypes::PROCESSED);
                             n++) {
                            auto& pprd =
                                getter(n, particles::MessageTypes::PROCESSED);
                            if (pprd.rdata(particles::RealData::timestamp) >=
                                rollback_timestamp) {

                                const int pair = static_cast<int>(
                                    pairing_function(pprd.cpu(), pprd.id()));
                                for (int m = 0;
                                     m <
                                     cnt_arr(
                                         iv,
                                         particles::MessageTypes::CONJUGATE);
                                     m++) {

                                    // This is a conjugate message that was
                                    // already treated, expect it to be an anti
                                    // message
                                    if (!getter.check(
                                            m, particles::MessageTypes::
                                                   CONJUGATE)) {
                                        getter.assert_different(
                                            m,
                                            particles::MessageTypes::CONJUGATE,
                                            particles::MessageTypes::
                                                ANTI_MESSAGE);
                                        continue;
                                    }

                                    auto& pcnj = getter(
                                        m, particles::MessageTypes::CONJUGATE);

                                    if (pair ==
                                        pcnj.idata(particles::IntData::pair)) {
                                        AMREX_ALWAYS_ASSERT(
                                            std::abs(
                                                pprd.rdata(particles::RealData::
                                                               timestamp) -
                                                pcnj.rdata(particles::RealData::
                                                               creation_time)) <
                                            constants::EPS);
                                        pcnj.idata(
                                            particles::IntData::type_id) =
                                            particles::MessageTypes::
                                                ANTI_MESSAGE;
                                        const auto piv =
                                            dom.atOffset(pcnj.idata(
                                                particles::IntData::receiver));
                                        AMREX_D_TERM(
                                            pcnj.pos(0) =
                                                plo[0] + (piv[0] + 0.5) * dx[0];
                                            ,
                                            pcnj.pos(1) =
                                                plo[1] + (piv[1] + 0.5) * dx[1];
                                            , pcnj.pos(2) =
                                                  plo[2] +
                                                  (piv[2] + 0.5) * dx[2];)
                                        AMREX_D_TERM(
                                            pcnj.idata(particles::IntData::i) =
                                                piv[0];
                                            ,
                                            pcnj.idata(particles::IntData::j) =
                                                piv[1];
                                            ,
                                            pcnj.idata(particles::IntData::k) =
                                                piv[2];)
                                    }
                                }

                                pprd.idata(particles::IntData::type_id) =
                                    particles::MessageTypes::MESSAGE;

                                // restore the state
                                sarr(iv, constants::LVT_IDX) =
                                    sarr(iv, constants::LVT_IDX) >
                                            pprd.rdata(particles::RealData::
                                                           old_timestamp)
                                        ? pprd.rdata(particles::RealData::
                                                         old_timestamp)
                                        : sarr(iv, constants::LVT_IDX);
                            }
                        }
                    }
                });
        }

        m_pc->Redistribute();
        m_pc->sort_particles();
        require_rollback = rollback.sum(0);
        iter++;
    }
    if (iter == max_iter) {
        amrex::Abort(
            "Maximum number of rollbacks reached (max_iter = " +
            std::to_string(max_iter) + ")");
    } else {
        rollback_statistics(lev);
    }

    m_pc->resolve_pairs();

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_pc->message_counts(lev).sum(particles::MessageTypes::ANTI_MESSAGE) ==
            0,
        "There should be no anti-messages left after rollback");
}

void SPADES::rollback_statistics(const int lev)
{
    BL_PROFILE("spades::SPADES::rollback_statistics()");

    auto max_rlb = static_cast<int>(m_state[lev].max(constants::RLB_IDX));
    amrex::ParallelAllReduce::Max<int>(
        max_rlb, amrex::ParallelContext::CommunicatorSub());

    amrex::Vector<amrex::Real> nrlbks(max_rlb + 1, 0);
    for (int n = 0; n < nrlbks.size(); n++) {
        nrlbks[n] = amrex::ReduceSum(
            m_state[lev], 0,
            [=] AMREX_GPU_HOST_DEVICE(
                const amrex::Box bx,
                const amrex::Array4<amrex::Real const>& sarr) -> amrex::Real {
                amrex::Real nrlbk = 0;
                amrex::Loop(
                    bx, [=, &nrlbk](
                            int i, int j, int AMREX_D_PICK(, , k)) noexcept {
                        const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                        nrlbk += (sarr(iv, constants::RLB_IDX) == n) ? 1 : 0;
                    });
                return nrlbk;
            });
    }
    amrex::ParallelDescriptor::ReduceRealSum(
        nrlbks.data(), static_cast<int>(nrlbks.size()),
        amrex::ParallelDescriptor::IOProcessorNumber());

    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "Rollback statistics at level " << lev << std::endl;
        for (int n = 0; n < nrlbks.size(); n++) {
            amrex::Print() << "  number of nodes doing " << n
                           << " rollbacks: " << nrlbks[n] << std::endl;
        }
        const auto nt = static_cast<amrex::Long>(
            std::accumulate(nrlbks.begin(), nrlbks.end(), 0.0));
        AMREX_ALWAYS_ASSERT(nt == boxArray(lev).numPts());
    }
}

void SPADES::update_gvt(const int lev)
{
    BL_PROFILE("spades::SPADES::update_gvt()");
    m_gvts[lev] = m_state[lev].min(constants::LVT_IDX, 0);
}

void SPADES::update_lbts(const int lev)
{
    BL_PROFILE("spades::SPADES::update_lbts()");
    amrex::MultiFab lbts(boxArray(lev), DistributionMap(lev), 1, 0);
    lbts.setVal(constants::LARGE_NUM);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = m_pc->MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        const auto& cnt_arr = m_pc->message_counts(lev).const_array(mfi);
        const auto& offsets_arr = m_pc->offsets(lev).const_array(mfi);
        const auto& lbts_arr = lbts.array(mfi);
        const auto index = std::make_pair(gid, tid);
        auto& pti = m_pc->GetParticles(lev)[index];
        auto& particles = pti.GetArrayOfStructs();
        auto* pstruct = particles().dataPtr();

        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(
                     int i, int j, int AMREX_D_PICK(, , k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const auto getter =
                    particles::Get(iv, cnt_arr, offsets_arr, pstruct);

                if (cnt_arr(iv, particles::MessageTypes::MESSAGE) > 0) {
                    auto& prcv = getter(0, particles::MessageTypes::MESSAGE);
                    lbts_arr(iv) = prcv.rdata(particles::RealData::timestamp);
                }
            });
    }

    m_lbts[lev] = lbts.min(0, 0);
    AMREX_ALWAYS_ASSERT(m_lbts[lev] >= m_gvts[lev]);
}

// a wrapper for EstTimeStep
void SPADES::compute_dt()
{
    BL_PROFILE("spades::SPADES::compute_dt()");
    amrex::Vector<amrex::Real> dt_tmp(finest_level + 1);

    for (int lev = 0; lev <= finest_level; ++lev) {
        dt_tmp[lev] = est_time_step(lev);
    }
    amrex::ParallelDescriptor::ReduceRealMin(
        dt_tmp.data(), static_cast<int>(dt_tmp.size()));

    constexpr amrex::Real change_max = 1.1;
    amrex::Real dt_0 = dt_tmp[0];
    int n_factor = 1;
    for (int lev = 0; lev <= finest_level; ++lev) {
        dt_tmp[lev] = std::min(dt_tmp[lev], change_max * m_dts[lev]);
        n_factor *= m_nsubsteps[lev];
        dt_0 = std::min(dt_0, n_factor * dt_tmp[lev]);
    }

    // Limit dt's by the value of stop_time.
    const amrex::Real eps = 1.e-3 * dt_0;
    if (m_ts_new[0] + dt_0 > m_stop_time - eps) {
        dt_0 = m_stop_time - m_ts_new[0];
    }

    m_dts[0] = dt_0;
    for (int lev = 1; lev <= finest_level; ++lev) {
        m_dts[lev] = m_dts[lev - 1] / m_nsubsteps[lev];
    }
}

// compute dt
amrex::Real SPADES::est_time_step(const int /*lev*/)
{
    BL_PROFILE("spades::SPADES::est_time_step()");
    return 1.0;
}

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
void SPADES::MakeNewLevelFromCoarse(
    int /*lev*/,
    amrex::Real /*time*/,
    const amrex::BoxArray& /*ba*/,
    const amrex::DistributionMapping& /*dm*/)
{
    BL_PROFILE("spades::SPADES::MakeNewLevelFromCoarse()");
    amrex::Abort("spades::SPADES::MakeNewLevelFromCoarse(): not implemented");
}

// Make a new level from scratch using provided BoxArray and
// DistributionMapping. Only used during initialization. overrides the pure
// virtual function in AmrCore
void SPADES::MakeNewLevelFromScratch(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    BL_PROFILE("spades::SPADES::MakeNewLevelFromScratch()");

    m_state[lev].define(
        ba, dm, constants::N_STATES, m_state_ngrow, amrex::MFInfo());

    m_ts_new[lev] = time;
    m_ts_old[lev] = constants::LOW_NUM;

    // Initialize the data
    initialize_state(lev);

    // Update particle container
    m_pc->Define(Geom(lev), dm, ba);
    m_pc->initialize_particles(m_lookahead);
    m_pc->initialize_state();
    m_pc->sort_particles();
}

void SPADES::initialize_state(const int lev)
{
    BL_PROFILE("spades::SPADES::initialize_state()");

    m_ic_op->initialize(lev, geom[lev].data());

    m_state[lev].FillBoundary(Geom(lev).periodicity());
}

// Remake an existing level using provided BoxArray and DistributionMapping
// and fill with existing fine and coarse data.
void SPADES::RemakeLevel(
    int /*lev*/,
    amrex::Real /* time*/,
    const amrex::BoxArray& /*ba*/,
    const amrex::DistributionMapping& /*dm*/)
{
    BL_PROFILE("spades::SPADES::RemakeLevel()");
    amrex::Abort("spades::SPADES::RemakeLevel(): not implemented");
}

// Delete level data
void SPADES::ClearLevel(int lev)
{
    BL_PROFILE("spades::SPADES::ClearLevel()");
    m_pc->clear_state();
    m_state[lev].clear();
    m_plt_mf[lev].clear();
}

void SPADES::set_ics()
{
    BL_PROFILE("spades::SPADES::set_ics()");
    if (m_ic_type == "constant") {
        m_ic_op = std::make_unique<ic::Initializer<ic::Constant>>(
            ic::Constant(ic::Constant()), m_state);
    } else {
        amrex::Abort(
            "spades::SPADES::set_ics(): User must specify a valid initial "
            "condition");
    }
}

// Check if a field exists
bool SPADES::check_field_existence(const std::string& name)
{
    BL_PROFILE("spades::SPADES::check_field_existence()");
    const auto vnames = {m_state_varnames, m_message_counts_varnames};
    return std::any_of(vnames.begin(), vnames.end(), [=](const auto& vn) {
        return get_field_component(name, vn) != -1;
    });
}

// Get field component
int SPADES::get_field_component(
    const std::string& name, const amrex::Vector<std::string>& varnames)
{
    BL_PROFILE("spades::SPADES::get_field_component()");
    const auto itr = std::find(varnames.begin(), varnames.end(), name);
    if (itr != varnames.cend()) {
        return static_cast<int>(std::distance(varnames.begin(), itr));
    }
    return -1;
}

// get a field based on a variable name
std::unique_ptr<amrex::MultiFab>
SPADES::get_field(const std::string& name, const int lev, const int ngrow)
{
    BL_PROFILE("spades::SPADES::get_field()");

    if (!check_field_existence(name)) {
        amrex::Abort(
            "spades::SPADES::get_field(): this field was not found: " + name);
    }

    const int nc = 1;
    std::unique_ptr<amrex::MultiFab> mf = std::make_unique<amrex::MultiFab>(
        boxArray(lev), DistributionMap(lev), nc, ngrow);

    const int srccomp_sd = get_field_component(name, m_state_varnames);
    if (srccomp_sd != -1) {
        amrex::MultiFab::Copy(*mf, m_state[lev], srccomp_sd, 0, nc, ngrow);
    }
    const int srccomp_id = get_field_component(name, m_message_counts_varnames);
    if (srccomp_id != -1) {
        auto const& cnt_arrs = m_pc->message_counts(lev).const_arrays();
        auto const& mf_arrs = mf->arrays();
        amrex::ParallelFor(
            *mf, mf->nGrowVect(), m_pc->message_counts(lev).nComp(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
                mf_arrs[nbx](i, j, k, n) = cnt_arrs[nbx](i, j, k, n);
            });
        amrex::Gpu::synchronize();
    }
    return mf;
}

amrex::Vector<std::string> SPADES::plot_file_var_names() const
{
    return m_spades_varnames;
}

std::string SPADES::plot_file_name(const int step) const
{
    return amrex::Concatenate(m_plot_file, step, m_file_name_digits);
}

std::string SPADES::chk_file_name(const int step) const
{
    return amrex::Concatenate(m_chk_file, step, m_file_name_digits);
}

// put together an array of multifabs for writing
amrex::Vector<const amrex::MultiFab*> SPADES::plot_file_mf()
{
    amrex::Vector<const amrex::MultiFab*> r;
    for (int lev = 0; lev <= finest_level; ++lev) {

        m_plt_mf[lev].define(
            boxArray(lev), DistributionMap(lev),
            static_cast<int>(plot_file_var_names().size()), 0);
        int cnt = 0;
        amrex::MultiFab::Copy(
            m_plt_mf[lev], m_state[lev], 0, cnt, m_state[lev].nComp(), 0);
        cnt += m_state[lev].nComp();
        auto const& cnt_arrs = m_pc->message_counts(lev).const_arrays();
        auto const& plt_mf_arrs = m_plt_mf[lev].arrays();
        amrex::ParallelFor(
            m_plt_mf[lev], m_plt_mf[lev].nGrowVect(),
            m_pc->message_counts(lev).nComp(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
                plt_mf_arrs[nbx](i, j, k, n + cnt) = cnt_arrs[nbx](i, j, k, n);
            });
        amrex::Gpu::synchronize();

        r.push_back(&m_plt_mf[lev]);
    }
    return r;
}

void SPADES::write_plot_file()
{
    BL_PROFILE("spades::SPADES::write_plot_file()");
    const std::string& plotfilename = plot_file_name(m_isteps[0]);
    const auto& mf = plot_file_mf();
    const auto& varnames = plot_file_var_names();

    amrex::Print() << "Writing plot file " << plotfilename << " at time "
                   << m_ts_new[0] << std::endl;

    amrex::WriteMultiLevelPlotfile(
        plotfilename, finest_level + 1, mf, varnames, Geom(), m_ts_new[0],
        m_isteps, refRatio());
    if (m_write_particles) {
        m_pc->write_plot_file(plotfilename);
    }

    write_info_file(plotfilename);
}

void SPADES::write_checkpoint_file() const
{
    BL_PROFILE("spades::SPADES::write_checkpoint_file()");
    const auto& varnames = m_state_varnames;

    // chk00010            write a checkpoint file with this root directory
    // chk00010/Header     this contains information you need to save (e.g.,
    // finest_level, t_new, etc.) and also
    //                     the BoxArrays at each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                these subdirectories will hold the MultiFab data
    // at each level of refinement

    const std::string& checkpointname = chk_file_name(m_isteps[0]);

    amrex::Print() << "Writing checkpoint file " << checkpointname
                   << " at time " << m_ts_new[0] << std::endl;

    const int nlevels = finest_level + 1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then
    // build
    // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    // write Header file
    if (amrex::ParallelDescriptor::IOProcessor()) {

        const std::string header_file_name(checkpointname + "/Header");
        amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);
        std::ofstream header_file;
        header_file.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        header_file.open(
            header_file_name.c_str(),
            std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);

        if (!header_file.good()) {
            amrex::FileOpenFailed(header_file_name);
        }

        header_file.precision(17);

        // write out title line
        header_file << "Checkpoint file for SPADES\n";

        // write out finest_level
        header_file << finest_level << "\n";

        // write out array of istep
        for (int m_istep : m_isteps) {
            header_file << m_istep << " ";
        }
        header_file << "\n";

        // write out array of dt
        for (double m_dt : m_dts) {
            header_file << m_dt << " ";
        }
        header_file << "\n";

        // write out array of t_new
        for (double i : m_ts_new) {
            header_file << i << " ";
        }
        header_file << "\n";

        // write the BoxArray at each level
        for (int lev = 0; lev <= finest_level; ++lev) {
            boxArray(lev).writeOn(header_file);
            header_file << '\n';
        }
    }

    // write the MultiFab data to, e.g., chk00010/Level_0/
    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::VisMF::Write(
            m_state[lev], amrex::MultiFabFileFullPrefix(
                              lev, checkpointname, "Level_", varnames[0]));
    }

    m_pc->Checkpoint(checkpointname, m_pc->identifier());

    write_info_file(checkpointname);

    write_rng_file(checkpointname);
}

void SPADES::read_checkpoint_file()
{
    BL_PROFILE("spades::SPADES::read_checkpoint_file()");
    const auto& varnames = m_state_varnames;

    amrex::Print() << "Restarting from checkpoint file " << m_restart_chkfile
                   << std::endl;

    // Header
    const std::string file(m_restart_chkfile + "/Header");

    amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::GetIOBufferSize());

    amrex::Vector<char> file_char_ptr;
    amrex::ParallelDescriptor::ReadAndBcastFile(file, file_char_ptr);
    std::string file_char_ptr_string(file_char_ptr.dataPtr());
    std::istringstream is(file_char_ptr_string, std::istringstream::in);

    std::string line, word;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level;
    goto_next_line(is);

    // read in array of istep
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            m_isteps[i++] = std::stoi(word);
        }
    }

    // read in array of dt
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            m_dts[i++] = std::stod(word);
        }
    }

    // read in array of t_new
    std::getline(is, line);
    {
        std::istringstream lis(line);
        int i = 0;
        while (lis >> word) {
            m_ts_new[i++] = std::stod(word);
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {

        // read in level 'lev' BoxArray from Header
        amrex::BoxArray ba;
        ba.readFrom(is);
        goto_next_line(is);

        // create a distribution mapping
        amrex::DistributionMapping dm{ba, amrex::ParallelDescriptor::NProcs()};

        // set BoxArray grids and DistributionMapping dmap in
        // AMReX_AmrMesh.H class
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);

        // build MultiFabs
        const int ncomp = static_cast<int>(varnames.size());
        AMREX_ALWAYS_ASSERT(ncomp == constants::N_STATES);
        m_state[lev].define(
            ba, dm, constants::N_STATES, m_state_ngrow, amrex::MFInfo());
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::VisMF::Read(
            m_state[lev], amrex::MultiFabFileFullPrefix(
                              lev, m_restart_chkfile, "Level_", varnames[0]));
        update_gvt(lev);
    }

    read_rng_file(m_restart_chkfile);

    init_particle_container();
    m_pc->Restart(m_restart_chkfile, m_pc->identifier());
    m_pc->initialize_state();
    m_pc->sort_particles();
}

void SPADES::write_info_file(const std::string& path) const
{
    BL_PROFILE("spades::SPADES::write_info_file()");
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    const std::string dash_line = "\n" + std::string(78, '-') + "\n";
    const std::string fname(path + "/spades_info");
    std::ofstream fh(fname.c_str(), std::ios::out);
    if (!fh.good()) {
        amrex::FileOpenFailed(fname);
    }

    fh << dash_line << "Grid information: " << std::endl;
    for (int lev = 0; lev < finestLevel() + 1; ++lev) {
        fh << "  Level: " << lev << "\n"
           << "    num. boxes = " << boxArray().size() << "\n"
           << "    maximum zones = ";

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            fh << Geom(lev).Domain().length(dir) << " ";
        }
        fh << "\n";
    }

    fh << dash_line << "Input file parameters: " << std::endl;
    amrex::ParmParse::dumpTable(fh, true);
    fh.close();
}

void SPADES::init_rng() const
{
    BL_PROFILE("spades::SPADES::init_rng()");

    if (m_seed > 0) {
        amrex::InitRandom(
            m_seed + amrex::ParallelDescriptor::MyProc(),
            amrex::ParallelDescriptor::NProcs(),
            m_seed + amrex::ParallelDescriptor::MyProc());
    } else if (m_seed == 0) {
        auto now = std::chrono::system_clock::now();
        auto now_ns =
            std::chrono::time_point_cast<std::chrono::nanoseconds>(now);
        auto rand_seed = now_ns.time_since_epoch().count();
        amrex::ParallelDescriptor::Bcast(
            &rand_seed, 1, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::InitRandom(
            rand_seed + amrex::ParallelDescriptor::MyProc(),
            amrex::ParallelDescriptor::NProcs(),
            rand_seed + amrex::ParallelDescriptor::MyProc());
    } else {
        amrex::Abort("RNG seed must be non-negative");
    }
}

void SPADES::write_rng_file(const std::string& path) const
{
    BL_PROFILE("spades::SPADES::write_rng_file()");

    const std::string base(path + "/rng");
    const std::string& fname = amrex::Concatenate(
        base, amrex::ParallelDescriptor::MyProc(), m_rng_file_name_digits);
    std::ofstream fh(fname.c_str(), std::ios::out);
    if (!fh.good()) {
        amrex::FileOpenFailed(fname);
    }
    amrex::SaveRandomState(fh);
}

void SPADES::read_rng_file(const std::string& path) const
{
    BL_PROFILE("spades::SPADES::read_rng_file()");

    if (m_seed < 0) {
        const std::string base(path + "/rng");
        const std::string& fname = amrex::Concatenate(
            base, amrex::ParallelDescriptor::MyProc(), m_rng_file_name_digits);
        amrex::Vector<char> file_char_ptr;
        read_file(fname, file_char_ptr, true);
        std::string file_str(file_char_ptr.dataPtr());
        std::istringstream is(file_str, std::istringstream::in);
        amrex::RestoreRandomState(is, 1, 0);
    } else {
        amrex::Warning(
            "Warning: Restarting without using the saved RNG seeds. Run with "
            "spades.seed=-1 to use those");
        AMREX_ASSERT_WITH_MESSAGE(
            m_seed < 0,
            "Use spades.seed=-1 when debugging and using restart files");
    }
}

} // namespace spades
