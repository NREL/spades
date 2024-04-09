#include <memory>

#include "SPADES.H"

namespace spades {
SPADES::SPADES()
{
    BL_PROFILE("SPADES::SPADES()");
    read_parameters();
    int nlevs_max = max_level + 1;

    if (max_level > 0) {
        amrex::Abort("Not supporting multilevel right now");
    }

    m_state_varnames.push_back("lvt");
    m_events_varnames.push_back("anti_message");
    m_events_varnames.push_back("undefined");
    m_events_varnames.push_back("message");
    m_events_varnames.push_back("processed");
    AMREX_ALWAYS_ASSERT(
        m_events_varnames.size() == particles::MessageTypes::ntypes);
    for (const auto& vname : m_state_varnames) {
        m_spades_varnames.push_back(vname);
    }
    for (const auto& vname : m_events_varnames) {
        m_spades_varnames.push_back(vname);
    }

    m_isteps.resize(nlevs_max, 0);
    m_nsubsteps.resize(nlevs_max, 1);
    for (int lev = 1; lev <= max_level; ++lev) {
        m_nsubsteps[lev] = MaxRefRatio(lev - 1);
    }

    m_gvts.resize(nlevs_max, constants::LOW_NUM);
    m_ts_new.resize(nlevs_max, 0.0);
    m_ts_old.resize(nlevs_max, constants::LOW_NUM);
    m_dts.resize(nlevs_max, constants::LARGE_NUM);

    m_state.resize(nlevs_max);
    m_events.resize(nlevs_max);
    m_plt_mf.resize(nlevs_max);
}

SPADES::~SPADES() = default;

void SPADES::init_data()
{
    BL_PROFILE("SPADES::init_data()");

    if (m_restart_chkfile.empty()) {
        // start simulation from the beginning
        const amrex::Real time = 0.0;
        set_ics();

        init_particle_containter();

        InitFromScratch(time);

        const int lev = 0;
        update_gvt(lev);

        compute_dt();

        count_events();

        if (m_chk_int > 0) {
            write_checkpoint_file();
        }
    } else {
        // restart from a checkpoint
        read_checkpoint_file();
        amrex::Abort("need to create then grab the particle chkp");
    }

    if (m_plot_int > 0) {
        write_plot_file();
    }

    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "Grid summary: " << std::endl;
        printGridSummary(amrex::OutStream(), 0, finest_level);
    }
}

void SPADES::init_particle_containter()
{
    BL_PROFILE("SPADES::init_particle_container()");
    m_pc =
        std::make_unique<particles::CellSortedParticleContainer>(GetParGDB());
}

void SPADES::read_parameters()
{
    BL_PROFILE("SPADES::read_parameters()");

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
    }

    {
        amrex::ParmParse pp("spades");
        pp.get("ic_type", m_ic_type);
        pp.query("write_particles", m_write_particles);
    }

    // force periodic bcs
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
        if (!amrex::DefaultGeometry().isPeriodic(dir)) {
            amrex::Abort("Geometry needs to be periodic");
        }
    }
}

void SPADES::evolve()
{
    BL_PROFILE("SPADES::evolve()");

    amrex::Real cur_time = m_ts_new[0];
    int last_plot_file_step = 0;

    for (int step = m_isteps[0]; step < m_max_step && cur_time < m_stop_time;
         ++step) {
        compute_dt();

        amrex::Print() << "\n==============================================="
                          "==============================="
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
            write_checkpoint_file();
        }

        if (cur_time >= m_stop_time - 1.e-6 * m_dts[0]) {
            break;
        }
    }
    if (m_plot_int > 0 && m_isteps[0] > last_plot_file_step) {
        write_plot_file();
    }
}

// advance a level by dt
// includes a recursive call for finer levels
void SPADES::time_step(
    const int lev, const amrex::Real time, const int iteration)
{
    BL_PROFILE("SPADES::time_step()");
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
                m_pc->update_cell_lists();
            }
        }
    }

    if (Verbose() != 0) {
        amrex::Print() << "[Level " << lev << " step " << m_isteps[lev] + 1
                       << "] ";
        amrex::Print() << "Advance with time = " << m_ts_new[lev]
                       << " dt = " << m_dts[lev] << std::endl;
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
        amrex::Print() << "[Level " << lev << " step " << m_isteps[lev] << "] ";
        amrex::Print() << "Advanced " << CountCells(lev) << " cells and "
                       << m_pc->TotalNumberOfParticles(lev) << " particles"
                       << std::endl;
    }
}

void SPADES::advance(
    const int lev,
    const amrex::Real /*time*/,
    const amrex::Real dt_lev,
    const int /*iteration*/,
    const int /*ncycle*/)
{
    BL_PROFILE("SPADES::advance()");

    process_events(lev);

    update_gvt(lev);

    m_pc->garbage_collect(m_gvts[lev]);

    m_pc->Redistribute();
    m_pc->sort_particles();
    m_pc->assign_cell_lists();

    m_ts_old[lev] = m_ts_new[lev]; // old time is now current time (time)
    m_ts_new[lev] += dt_lev;       // new time is ahead
}

void SPADES::post_time_step()
{
    BL_PROFILE("SPADES::post_time_step()");

    count_events();
}

void SPADES::count_events()
{
    BL_PROFILE("SPADES::count_events()");

    const int lev = 0;
    m_events[lev].setVal(0);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (amrex::MFIter mfi = m_pc->MakeMFIter(lev); mfi.isValid(); ++mfi) {

        const amrex::Box& box = mfi.tilebox();
        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        const auto& events_arr = m_events[lev].array(mfi);
        auto& pti = m_pc->GetParticles(lev)[std::make_pair(gid, tid)];
        const auto& particles = pti.GetArrayOfStructs();
        const auto* pstruct = particles().dataPtr();
        const int np = pti.numParticles();

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pindex) noexcept {
            const auto& p = pstruct[pindex];
            const amrex::IntVect iv(AMREX_D_DECL(
                p.idata(particles::IntData::i), p.idata(particles::IntData::j),
                p.idata(particles::IntData::k)));

            if (box.contains(iv)) {
                amrex::Gpu::Atomic::AddNoRet(
                    &events_arr(iv, p.idata(particles::IntData::type_id)), 1);
            }
        });
    }
}

void SPADES::process_events(const int lev)
{
    BL_PROFILE("SPADES::process_events()");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_state_ngrow == m_pc->ngrow(),
        "Particle and state cells must be equal for now");

    const auto& plo = Geom(lev).ProbLoArray();
    const auto& dx = Geom(lev).CellSizeArray();
    const auto& dlo = Geom(lev).Domain().smallEnd();
    const auto& dhi = Geom(lev).Domain().bigEnd();
    const auto lookahead = m_lookahead;

    // #ifdef _OPENMP
    // #pragma omp parallel
    // #endif
    //     for (amrex::MFIter mfi = m_pc->MakeMFIter(lev); mfi.isValid(); ++mfi)
    //     {
    //         const amrex::Box& box = mfi.tilebox();
    //         const int gid = mfi.index();
    //         const int tid = mfi.LocalTileIndex();
    //         const auto& sarr = m_state[lev].array(mfi);
    //         const auto index = std::make_pair(gid, tid);
    //         const auto& cell_lists = m_pc->cell_lists(lev, index);
    //         auto& pti = m_pc->GetParticles(lev)[index];
    //         auto& particles = pti.GetArrayOfStructs();
    //         auto* pstruct = particles().dataPtr();

    //         const auto* p_cell_list = cell_lists.list().data();
    //         const auto* p_cell_counts = cell_lists.counts().data();
    //         const auto* p_cell_offsets = cell_lists.offsets().data();

    //         amrex::ParallelForRNG(
    //             box, [=] AMREX_GPU_DEVICE(
    //                      int i, int j, int k,
    //                      amrex::RandomEngine const& engine) noexcept {
    //                 const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
    //                 const auto ni = box.index(iv);
    //                 const auto* l_cell_list =
    //                 &p_cell_list[p_cell_offsets[ni]]; const auto np =
    //                 p_cell_counts[ni];

    //                 for (int pidx = 0; pidx < np; pidx++) {
    //                     particles::CellSortedParticleContainer::ParticleType&
    //                     p =
    //                         pstruct[l_cell_list[pidx]];

    //                     amrex::Print()
    //                         << "particle: " << pidx << " in cell: " << iv
    //                         << " has timestamp: "
    //                         << p.rdata(particles::RealData::timestamp)
    //                         << " and type " <<
    //                         p.idata(particles::IntData::type_id)
    //                         << std::endl;

    //                     // // only process messages, fixme: get a pointer to
    //                     the
    //                     // // messages and skip the test
    //                     // if (p.idata(particles::IntData::type_id) ==
    //                     //     particles::MessageTypes::message) {

    //                     //     if (sarr(iv, constants::LVT_IDX) >
    //                     //         p.rdata(particles::RealData::timestamp)) {
    //                     //         amrex::Print()
    //                     //             << "LVT: " << sarr(iv,
    //                     constants::LVT_IDX)
    //                     //             << " and message time stamp is "
    //                     //             <<
    //                     p.rdata(particles::RealData::timestamp)
    //                     //             << std::endl;
    //                     //         // amrex::Abort("I got a message from the
    //                     past");
    //                     //         amrex::Print()
    //                     //             << "I got a message from the past" <<
    //                     std::endl;
    //                     //     }

    //                     //     // process the event
    //                     //     sarr(iv, constants::LVT_IDX) =
    //                     //         p.rdata(particles::RealData::timestamp);
    //                     //     p.idata(particles::IntData::type_id) =
    //                     //         particles::MessageTypes::processed;

    //                     //     // find undefined event to add to the queue
    //                     //     // fixme just get a pointer to this thing
    //                     //     // fixme add abort if you can't find an
    //                     undefined
    //                     //     // particles
    //                     //     int pidx_undef1 = -1;
    //                     //     for (int pidx2 = 0; pidx2 < np; pidx2++) {
    //                     //         particles::CellSortedParticleContainer::
    //                     //             ParticleType& pl =
    //                     pstruct[l_cell_list[pidx2]];
    //                     //         if (pl.idata(particles::IntData::type_id)
    //                     ==
    //                     //             particles::MessageTypes::undefined) {
    //                     //             pidx_undef1 = pidx2;
    //                     //             break;
    //                     //         }
    //                     //     }
    //                     //     AMREX_ALWAYS_ASSERT(pidx_undef1 > -1);

    //                     //
    //                     particles::CellSortedParticleContainer::ParticleType&
    //                     //         p2 = pstruct[l_cell_list[pidx_undef1]];
    //                     //     amrex::IntVect iv_dest(AMREX_D_DECL(
    //                     //         amrex::Random_int(dhi[0] - dlo[0] + 1) +
    //                     dlo[0],
    //                     //         amrex::Random_int(dhi[1] - dlo[1] + 1) +
    //                     dlo[1],
    //                     //         amrex::Random_int(dhi[2] - dlo[2] + 1) +
    //                     dlo[2]));
    //                     //     amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>
    //                     pos = {
    //                     //         AMREX_D_DECL(
    //                     //             plo[0] + (iv_dest[0] + 0.5) * dx[0],
    //                     //             plo[1] + (iv_dest[1] + 0.5) * dx[1],
    //                     //             plo[2] + (iv_dest[2] + 0.5) * dx[2])};
    //                     //     const amrex::Real ts = sarr(iv,
    //                     constants::LVT_IDX) +
    //                     //                            random_exponential(1.0)
    //                     +
    //                     //                            lookahead;

    //                     //     particles::Create()(
    //                     //         p2, ts, pos, iv_dest, box.index(iv),
    //                     //         box.index(iv_dest));

    //                     //     // find _another_ undefined particle for the
    //                     //     // antimessage fixme: oooff this is ugly
    //                     //     int pidx_undef2 = -1;
    //                     //     for (int pidx2 = 0; pidx2 < np; pidx2++) {
    //                     //         particles::CellSortedParticleContainer::
    //                     //             ParticleType& pl =
    //                     pstruct[l_cell_list[pidx2]];
    //                     //         if (pl.idata(particles::IntData::type_id)
    //                     ==
    //                     //             particles::MessageTypes::undefined) {
    //                     //             pidx_undef2 = pidx2;
    //                     //             break;
    //                     //         }
    //                     //     }
    //                     //     AMREX_ALWAYS_ASSERT(pidx_undef2 > -1);
    //                     //
    //                     particles::CellSortedParticleContainer::ParticleType&
    //                     //         p3 = pstruct[l_cell_list[pidx_undef2]];

    //                     //     // fixme, could do a copy. Or just pass p.pos
    //                     to Create
    //                     //     amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>
    //                     anti_pos =
    //                     //         {AMREX_D_DECL(
    //                     //             plo[0] + (iv[0] + 0.5) * dx[0],
    //                     //             plo[1] + (iv[1] + 0.5) * dx[1],
    //                     //             plo[2] + (iv[2] + 0.5) * dx[2])};

    //                     //     // This is weird. The antimessage
    //                     //     // position is iv but the receiver is
    //                     //     // still updated (we need to know who to
    //                     //     // send this to)
    //                     //     particles::Create()(
    //                     //         p3, ts, anti_pos, iv_dest, box.index(iv),
    //                     //         box.index(iv_dest));
    //                     //     p3.idata(particles::IntData::type_id) =
    //                     //         particles::MessageTypes::anti_message;

    //                     //     break; // only do 1 event
    //                     // }
    //                 }
    //             });
    //     }
}

void SPADES::update_gvt(const int lev)
{
    BL_PROFILE("SPADES::update_gvt()");
    m_gvts[lev] = m_state[lev].min(0, 0);
}

// a wrapper for EstTimeStep
void SPADES::compute_dt()
{
    BL_PROFILE("SPADES::compute_dt()");
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
    BL_PROFILE("SPADES::est_time_step()");
    return 1.0;
}

// Make a new level using provided BoxArray and DistributionMapping and
// fill with interpolated coarse level data.
void SPADES::MakeNewLevelFromCoarse(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    BL_PROFILE("SPADES::MakeNewLevelFromCoarse()");
    amrex::Abort("MakeNewLevelFromCoarse not implemented");
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
    BL_PROFILE("SPADES::MakeNewLevelFromScratch()");

    m_state[lev].define(
        ba, dm, constants::N_STATES, m_state_ngrow, amrex::MFInfo());

    m_events[lev].define(
        ba, dm, particles::MessageTypes::ntypes, m_state[lev].nGrow(),
        amrex::MFInfo());

    m_ts_new[lev] = time;
    m_ts_old[lev] = constants::LOW_NUM;

    // Initialize the data
    initialize_state(lev);
    m_events[lev].setVal(0);

    // Update particle container
    m_pc->Define(Geom(lev), dm, ba);
    m_pc->initialize(m_lookahead);
    m_pc->sort_particles();
    m_pc->update_cell_lists();
}

void SPADES::initialize_state(const int lev)
{
    BL_PROFILE("SPADES::initialize_f()");

    m_ic_op->initialize(lev, geom[lev].data());

    m_state[lev].FillBoundary(Geom(lev).periodicity());
}

// Remake an existing level using provided BoxArray and DistributionMapping
// and fill with existing fine and coarse data.
void SPADES::RemakeLevel(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm)
{
    BL_PROFILE("SPADES::RemakeLevel()");
    amrex::Abort("RemakeLevel not implemented");
}

// Delete level data
void SPADES::ClearLevel(int lev)
{
    BL_PROFILE("SPADES::ClearLevel()");
    m_state[lev].clear();
    m_events[lev].clear();
    m_plt_mf[lev].clear();
}

void SPADES::set_ics()
{
    BL_PROFILE("SPADES::set_ics()");
    if (m_ic_type == "constant") {
        m_ic_op = std::make_unique<ic::Initializer<ic::Constant>>(
            ic::Constant(ic::Constant()), m_state);
    } else {
        amrex::Abort(
            "SPADES::set_ics(): User must specify a valid initial condition");
    }
}

// Check if a field exists
bool SPADES::check_field_existence(const std::string& name)
{
    BL_PROFILE("SPADES::check_field_existence()");
    const auto vnames = {m_state_varnames, m_events_varnames};
    return std::any_of(vnames.begin(), vnames.end(), [=](const auto& vn) {
        return get_field_component(name, vn) != -1;
    });
}

// Get field component
int SPADES::get_field_component(
    const std::string& name, const amrex::Vector<std::string>& varnames)
{
    BL_PROFILE("SPADES::get_field_component()");
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
    BL_PROFILE("SPADES::get_field()");

    if (!check_field_existence(name)) {
        amrex::Abort("SPADES::get_field(): this field was not found: " + name);
    }

    const int nc = 1;
    std::unique_ptr<amrex::MultiFab> mf = std::make_unique<amrex::MultiFab>(
        boxArray(lev), DistributionMap(lev), nc, ngrow);

    const int srccomp_sd = get_field_component(name, m_state_varnames);
    if (srccomp_sd != -1) {
        amrex::MultiFab::Copy(*mf, m_state[lev], srccomp_sd, 0, nc, ngrow);
    }
    const int srccomp_id = get_field_component(name, m_events_varnames);
    if (srccomp_id != -1) {
        auto const& events_arrs = m_events[lev].const_arrays();
        auto const& mf_arrs = mf->arrays();
        amrex::ParallelFor(
            *mf, mf->nGrowVect(), m_events[lev].nComp(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
                mf_arrs[nbx](i, j, k, n) = events_arrs[nbx](i, j, k, n);
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
        auto const& events_arrs = m_events[lev].const_arrays();
        auto const& plt_mf_arrs = m_plt_mf[lev].arrays();
        amrex::ParallelFor(
            m_plt_mf[lev], m_plt_mf[lev].nGrowVect(), m_events[lev].nComp(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
                plt_mf_arrs[nbx](i, j, k, n + cnt) =
                    events_arrs[nbx](i, j, k, n);
            });
        amrex::Gpu::synchronize();

        r.push_back(&m_plt_mf[lev]);
    }
    return r;
}

void SPADES::write_plot_file()
{
    BL_PROFILE("SPADES::write_plot_file()");
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
    BL_PROFILE("SPADES::write_checkpoint_file()");
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

    write_info_file(checkpointname);
}

void SPADES::read_checkpoint_file()
{
    BL_PROFILE("SPADES::read_checkpoint_file()");
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
        AMREX_ASSERT(ncomp == constants::N_STATES);
        m_state[lev].define(
            ba, dm, constants::N_STATES, m_state_ngrow, amrex::MFInfo());
    }

    // read in the MultiFab data
    for (int lev = 0; lev <= finest_level; ++lev) {
        amrex::VisMF::Read(
            m_state[lev], amrex::MultiFabFileFullPrefix(
                              lev, m_restart_chkfile, "Level_", varnames[0]));
    }

    // Populate the other data
    for (int lev = 0; lev <= finest_level; ++lev) {
        m_state[lev].FillBoundary(Geom(lev).periodicity());
        m_state[lev].setVal(0.0);
    }
}

void SPADES::write_info_file(const std::string& path) const
{
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    const std::string dash_line = "\n" + std::string(78, '-') + "\n";
    const std::string fname(path + "/poudre_info");
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

// utility to skip to next line in Header
void SPADES::goto_next_line(std::istream& is)
{
    constexpr std::streamsize bl_ignore_max{100000};
    is.ignore(bl_ignore_max, '\n');
}
} // namespace spades
