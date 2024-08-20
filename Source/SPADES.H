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
   @brief Scalable PArallel Discrete Event Simulation
 **/
namespace spades {

//! Main SPADES class
class SPADES : public amrex::AmrCore
{
public:
    SPADES();
    ~SPADES() override;

    //! Initializes multilevel data
    void init_data();

    //! Advance solution to final time
    void evolve();

    /**
       @brief Make a new level

       Make a new level using provided BoxArray and DistributionMapping
       and fill with interpolated coarse level data. Overrides the pure
       virtual function in AmrCore.

       @param lev [in] level
       @param time [in] time
       @param ba [in] box array
       @param dm [in] distribution map
     **/
    void MakeNewLevelFromCoarse(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;

    /**
       @brief Remake an existing level

       Remake an existing level using provided BoxArray and
       DistributionMapping and fill with existing fine and coarse
       data. Overrides the pure virtual function in AmrCore.

       @param lev [in] level
       @param time [in] time
       @param ba [in] box array
       @param dm [in] distribution map
     **/
    void RemakeLevel(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;

    /**
       @brief Delete level data

       Overrides the pure virtual function in AmrCore.

       @param lev [in] level
     **/
    void ClearLevel(int lev) override;

    /**
       @brief Make a level from scratch

       Make a new level from scratch using provided BoxArray and
       DistributionMapping. Only used during initialization. Overrides
       the pure virtual function in AmrCore

       @param lev [in] level
       @param time [in] time
       @param ba [in] box array
       @param dm [in] distribution map
     **/
    void MakeNewLevelFromScratch(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;

    /**
       @brief Tag cells for refinement

       Overrides the pure virtual function in AmrCore
     **/
    void ErrorEst(
        int /*lev*/,
        amrex::TagBoxArray& /*tags*/,
        amrex::Real /*time*/,
        int /*ngrow*/) override
    {}

    /**
       @brief Compute the time step
       @param lev [in] level
       @return the time step
     **/
    static amrex::Real est_time_step(const int lev);

    /**
       @brief Get a field based on a variable name
       @param name [in] field name
       @param lev [in] level
       @param ngrow [in] number of grow cells
       @return the requested field
     **/
    std::unique_ptr<amrex::MultiFab>
    get_field(const std::string& name, const int lev, const int ngrow);

    /**
       @brief Put together an array of MultiFab for writing
       @return array of MultiFab
     **/
    //!
    void plot_file_mf();

    /**
       @brief Initialize state
       @param lev [in] level
     **/
    void initialize_state(const int lev);

    //! Initialize the particle container
    void init_particle_container();

    /**
       @brief Process messages
       @param lev [in] level
     **/
    void process_messages(const int lev);

    /**
       @brief Perform rollback
       @param lev [in] level
     **/
    void rollback(const int lev);

    /**
       @brief Print rollback statistics
       @param lev [in] level
     **/
    void rollback_statistics(const int lev);

    /**
       @brief Update the global virtual time
       @param lev [in] level
     **/
    void update_gvt(const int lev);

    /**
       @brief Update the Lower Bound on Incoming Time Stamp
       @param lev [in] level
     **/
    void update_lbts(const int lev);

private:
    //! Read parameters
    void read_parameters();

    //! Wrapper for EstTimeStep
    void compute_dt();

    /**
       @brief Advance a level by the time step

       Includes a recursive call for finer levels

       @param lev [in] level
       @param time [in] time
       @param iteration [in] iteration
     **/
    void time_step(const int lev, const amrex::Real time, const int iteration);

    /**
       @brief Advance a single level for a single time step
       @param lev [in] level
       @param time [in] time
       @param dt_lev [in] time step for this level
       @param iteration [in] iteration
       @param ncycle [in] number of cycles for this level
     **/
    void advance(
        const int lev,
        const amrex::Real time,
        const amrex::Real dt_lev,
        const int iteration,
        const int ncycle);

    //! Perform work after a time step
    static void post_time_step();

    /**
       @brief Check if a field exists
       @param name [in] field name
       @return True if field of \p name exists
     **/
    bool check_field_existence(const std::string& name);

    /**
       @brief Get field component
       @param name [in] field name
       @param varnames [in] vector of the field names
       @return field component
     **/
    static int get_field_component(
        const std::string& name, const amrex::Vector<std::string>& varnames);

    //! Set the user defined IC functions
    void set_ics();

    /**
       @brief Get plotfile name
       @param step [in] current time step
       @return plotfile name
     **/
    [[nodiscard]] std::string plot_file_name(const int step) const;

    /**
       @brief Get checkpoint file name
       @param step [in] current time step
       @return chkfile name
     **/
    [[nodiscard]] std::string chk_file_name(const int step) const;

    /**
       @brief Set plotfile variables names
       @return vector of variable names
     **/
    [[nodiscard]] amrex::Vector<std::string> plot_file_var_names() const;

    //! Write plotfile to disk
    void write_plot_file();

    //! Write checkpoint file to disk
    void write_checkpoint_file() const;

    //! Read checkpoint file from disk
    void read_checkpoint_file();

    /**
       @brief Write job info to disk
       @param path [in] path for file
     **/
    void write_info_file(const std::string& path) const;

    /**
       @brief Write random number generator seed info
       @param path [in] path for file
     **/
    void write_rng_file(const std::string& path) const;

    /**
       @brief Read random number generator seed info
       @param path [in] path for file
     **/
    void read_rng_file(const std::string& path) const;

    //! Initialize the random number generator
    void init_rng() const;

    /**
       @brief Write simulation information
       @param is_init [in] boolean indicating if this is the initializing step
     **/
    void write_data_file(const bool is_init) const;

    /**
       @brief Compute summary data for this level
       @param lev [in] level
     **/
    void level_summary(const int lev);

    //! Current steps at each level
    amrex::Vector<int> m_isteps;

    //! Number of substeps at each level
    amrex::Vector<int> m_nsubsteps;

    //! Global virtual time at each level
    amrex::Vector<amrex::Real> m_gvts;

    //! New time at each level
    amrex::Vector<amrex::Real> m_ts_new;

    //! Old time at each level
    amrex::Vector<amrex::Real> m_ts_old;

    //! Time step at each level
    amrex::Vector<amrex::Real> m_dts;

    //! Lower bound on incoming time stamp at each level
    amrex::Vector<amrex::Real> m_lbts;

    //! Cell count at each level
    amrex::Vector<amrex::Long> m_ncells;

    //! Total message count at each level
    amrex::Vector<amrex::Long> m_ntotal_messages;

    //! Message count at each level
    amrex::Vector<amrex::Long> m_nmessages;

    //! Count of processed messages at each level
    amrex::Vector<int> m_nprocessed_messages;

    //! Number of rollbacks at each level
    amrex::Vector<amrex::Vector<int>> m_nrollbacks;

    //! Min timings for each step at each level
    amrex::Vector<amrex::Real> m_min_timings;

    //! Max timings for each step at each level
    amrex::Vector<amrex::Real> m_max_timings;

    //! Average timings for each step at each level
    amrex::Vector<amrex::Real> m_avg_timings;

    //! Maximum number of steps
    int m_max_step = std::numeric_limits<int>::max();

    //! Stop time
    amrex::Real m_stop_time = std::numeric_limits<amrex::Real>::max();

    //! All variable names for output
    amrex::Vector<std::string> m_spades_varnames;

    //! State variable names for output
    amrex::Vector<std::string> m_state_varnames;

    //! Message count names for output
    amrex::Vector<std::string> m_message_counts_varnames;

    //! Multifabs to store the solution at each level of refinement
    amrex::Vector<amrex::MultiFab> m_state;

    //! Multifabs for output at each level of refinement
    amrex::Vector<amrex::MultiFab> m_plt_mf;

    //! Number of ghost cells
    const int m_state_ngrow = 0;

    //! Regrid interval (optional user input)
    int m_regrid_int = 2;

    //! Restart file name, restart from this checkpoint if it is not empty
    std::string m_restart_chkfile;

    //! Plotfile prefix (optional user input)
    std::string m_plot_file{"plt"};

    //! Plotfile frequency (optional user input)
    int m_plot_int = -1;

    //! Checkpoint prefix (optional user input)
    std::string m_chk_file{"chk"};

    //! Checkpoint frequency (optional user input)
    int m_chk_int = -1;

    //! Digits used in the plotfile and checkpoint file names
    int m_file_name_digits{constants::FILE_NAME_DIGITS};

    //! Digits used in the rng seed file names
    int m_rng_file_name_digits{constants::RNG_FILE_NAME_DIGITS};

    //! Boolean to output particles (optional user input)
    bool m_write_particles{false};

    //! Filename for simulation data
    std::string m_data_fname{"data.csv"};

    //! Data precision for data output
    const int m_data_precision{6};

    //! Lookahead value (optional user input)
    amrex::Real m_lookahead{1.0};

    //! Window size for processing messages (optional user input)
    amrex::Real m_window_size{constants::LARGE_NUM};

    //! Number of messages to process in each step (optional user input)
    int m_messages_per_step{1};

    //! Random number generator seed (optional user input)
    int m_seed{0};

    //! Initial condition operator
    std::unique_ptr<ic::InitializerBase> m_ic_op;

    //! Initial condition type (optional user input)
    std::string m_ic_type;

    //! Particle container
    std::unique_ptr<particles::CellSortedParticleContainer> m_pc;
};
} // namespace spades
#endif
