#ifndef MESSAGEPARTICLECONTAINER_H
#define MESSAGEPARTICLECONTAINER_H
#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_Random.H>
#include "SpadesParticleContainer.H"
#include "MessageData.H"
#include "MessageOps.H"
#include "Utilities.H"

/**
   @brief SPADES particles
 **/
namespace spades::particles {

//! Main SPADES message container
class MessageParticleContainer
    : public SpadesParticleContainer<
          MessageTypes::NTYPES,
          0,
          0,
          MessageRealData::ncomps,
          MessageIntData::ncomps>
{
public:
    /**
       @brief Class identifier name
       @return class identifier
     **/
    static std::string identifier() { return "messages"; }

    /**
       @brief Constructor
       @param par_gdb [in] particle database
       @param ngrow [in] number of grow cells
     **/
    explicit MessageParticleContainer(amrex::AmrParGDB* par_gdb, int ngrow = 0);

    /**
       @brief Constructor
       @param geom [in] geometry
       @param dmap [in] distribution map
       @param ba [in] box array
       @param ngrow [in] number of grow cells
     **/
    explicit MessageParticleContainer(
        const amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<amrex::DistributionMapping>& dmap,
        const amrex::Vector<amrex::BoxArray>& ba,
        int ngrow = 0);

    /**
       @brief Initialize the messages
       @param lookahead [in] lookahead
     **/
    void initialize_messages(const amrex::Real lookahead);

    //! Sort the messages
    void sort() override;

    //! Update the undefined messages
    void update_undefined();

    //! Resolve message pairs (remove message/anti-message pairs)
    void resolve_pairs();

    //! Compute the minimum time stamp of the messages
    amrex::Real compute_gvt();

    /**
       @brief Perform garbage collection
       @param gvt [in] global virtual time
     **/
    void garbage_collect(const amrex::Real gvt);

    void write_plot_file(const std::string& plt_filename) override
    {
        write_plot_file_impl(plt_filename, identifier());
    };

    //! Read user parameters
    void read_parameters() override;

    //! Initialize variable names
    void initialize_variable_names() override;

private:
    //! Initial messages per LP
    int m_messages_per_lp{1};

    //! Width of exponential distribution (optional user input)
    //! smaller values of lambda = larger variance in random values
    //! larger values of lambda = smaller variance in random values
    amrex::Real m_lambda{1.0};
};
} // namespace spades::particles
#endif
