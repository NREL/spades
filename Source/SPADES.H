#ifndef SPADES_H
#define SPADES_H
#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ErrorList.H>
#include <AMReX_Random.H>
#include "Constants.H"
#include "IC.H"
#include "Utilities.H"
#include "CellSortedParticleContainer.H"

/**
   @brief Scalable PArallel Discrete Events Simulation
 **/
namespace spades {

class SPADES : public amrex::AmrCore
{
public:
    SPADES();
    ~SPADES() override;

    // initializes multilevel data
    void init_data();

    // advance solution to final time
    void evolve();

    // Make a new level using provided BoxArray and DistributionMapping and
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
    void MakeNewLevelFromCoarse(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping
    // and fill with existing fine and coarse data. overrides the pure virtual
    // function in AmrCore
    void RemakeLevel(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;

    // Delete level data
    // overrides the pure virtual function in AmrCore
    void ClearLevel(int lev) override;

    // Make a new level from scratch using provided BoxArray and
    // DistributionMapping. Only used during initialization. overrides the pure
    // virtual function in AmrCore
    void MakeNewLevelFromScratch(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;

    // tag all cells for refinement
    // overrides the pure virtual function in AmrCore
    void ErrorEst(
        int /*lev*/,
        amrex::TagBoxArray& /*tags*/,
        amrex::Real /*time*/,
        int /*ngrow*/) override
    {}

    // compute dt
    static amrex::Real est_time_step(const int lev);

    // get a field based on a variable name
    std::unique_ptr<amrex::MultiFab>
    get_field(const std::string& name, const int lev, const int ngrow);

    // put together an array of multifabs for writing
    amrex::Vector<const amrex::MultiFab*> plot_file_mf();

    // Initialize state
    void initialize_state(const int lev);

    // Initialize the particle container
    void init_particle_container();

    // process messages
    void process_messages(const int lev);

    // rollback
    void rollback(const int lev);

    // print rollback statistics
    void rollback_statistics(const int lev);

    // update gvt
    void update_gvt(const int lev);

    // update lbts
    void update_lbts(const int lev);

private:
    // read parameters
    void read_parameters();

    // a wrapper for EstTimeStep(0
    void compute_dt();

    // advance a level by dt
    void time_step(const int lev, const amrex::Real time, const int iteration);

    // advance a single level for a single time step, updates flux registers
    void advance(
        const int lev,
        const amrex::Real time,
        const amrex::Real dt_lev,
        const int iteration,
        const int ncycle);

    // perform work after a time step
    static void post_time_step();

    // Check if a field exists
    bool check_field_existence(const std::string& name);

    // Get field component
    static int get_field_component(
        const std::string& name, const amrex::Vector<std::string>& varnames);

    // Set the user defined IC functions
    void set_ics();

    // get plotfile name
    [[nodiscard]] std::string plot_file_name(const int step) const;

    // get chkfile name
    [[nodiscard]] std::string chk_file_name(const int step) const;

    // set plotfile variables names
    [[nodiscard]] amrex::Vector<std::string> plot_file_var_names() const;

    // write plotfile to disk
    void write_plot_file();

    // write checkpoint file to disk
    void write_checkpoint_file() const;

    // read checkpoint file from disk
    void read_checkpoint_file();

    // write job info to disk
    void write_info_file(const std::string& /*path*/) const;

    // write rng info
    void write_rng_file(const std::string& /*path*/) const;

    // read rng info
    void read_rng_file(const std::string& /*path*/) const;

    // initialize the rng
    void init_rng() const;

    // write simulation information
    void write_data_file(const bool is_init) const;

    // compute summary data
    void level_summary(const int lev);

    amrex::Vector<int> m_isteps;    // which step?
    amrex::Vector<int> m_nsubsteps; // how many substeps on each level?

    // keep track of old time, new time, and time step at each level
    amrex::Vector<amrex::Real> m_gvts;
    amrex::Vector<amrex::Real> m_ts_new;
    amrex::Vector<amrex::Real> m_ts_old;
    amrex::Vector<amrex::Real> m_dts;

    // lower bound on incoming time stamp
    amrex::Vector<amrex::Real> m_lbts;

    // cell count at each level
    amrex::Vector<amrex::Long> m_ncells;

    // total message count at each level
    amrex::Vector<amrex::Long> m_ntotal_messages;

    //  message count at each level
    amrex::Vector<amrex::Long> m_nmessages;

    // count of processed messages
    amrex::Vector<int> m_nprocessed_messages;

    // number of rollbacks at each level
    amrex::Vector<amrex::Vector<int>> m_nrollbacks;

    // timings for each step
    amrex::Vector<amrex::Real> m_min_timings;
    amrex::Vector<amrex::Real> m_max_timings;
    amrex::Vector<amrex::Real> m_avg_timings;

    // maximum number of steps and stop time
    int m_max_step = std::numeric_limits<int>::max();
    amrex::Real m_stop_time = std::numeric_limits<amrex::Real>::max();

    // variable names for output
    amrex::Vector<std::string> m_spades_varnames;
    amrex::Vector<std::string> m_state_varnames;
    amrex::Vector<std::string> m_message_counts_varnames;

    // array of multifabs to store the solution at each level of refinement
    amrex::Vector<amrex::MultiFab> m_state;
    amrex::Vector<amrex::MultiFab> m_plt_mf;

    const int m_state_ngrow = 0;

    // how often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int m_regrid_int = 2;

    // if >= 0 we restart from a checkpoint
    std::string m_restart_chkfile;

    // plotfile prefix and frequency
    std::string m_plot_file{"plt"};

    int m_plot_int = -1;

    // checkpoint prefix and frequency
    std::string m_chk_file{"chk"};
    int m_chk_int = -1;

    // digits used in chk and plt file names
    int m_file_name_digits{constants::FILE_NAME_DIGITS};
    int m_rng_file_name_digits{constants::RNG_FILE_NAME_DIGITS};

    // output particles
    bool m_write_particles{false};

    // filename for simulation data
    const std::string m_data_fname{"data.csv"};

    // data precision for data output
    const int m_data_precision{6};

    amrex::Real m_lookahead{1.0};
    amrex::Real m_window_size{constants::LARGE_NUM};
    int m_messages_per_step{1};

    int m_seed{0};

    std::unique_ptr<ic::InitializerBase> m_ic_op;
    std::string m_ic_type;

    std::unique_ptr<particles::CellSortedParticleContainer> m_pc;
};
} // namespace spades
#endif
