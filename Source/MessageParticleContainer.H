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
          MessageTypes,
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

protected:
    //! Initialize variable names
    void initialize_variable_names() override;

private:
    //! Lower bound for undefined message counts
    const int m_lower_undefined_count{40};

    //! Reset value for undefined message counts
    const int m_reset_undefined_count{50};

    //! Upper bound for undefined message counts
    const int m_upper_undefined_count{100};
};
} // namespace spades::particles
#endif
