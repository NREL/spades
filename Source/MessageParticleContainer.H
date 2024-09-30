#ifndef MESSAGEPARTICLECONTAINER_H
#define MESSAGEPARTICLECONTAINER_H
#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_NeighborParticles.H>
#include <AMReX_Random.H>
#include "MessageData.H"
#include "MessageOps.H"
#include "Utilities.H"

/**
   @brief SPADES particles
 **/
namespace spades::particles {

//! Information for a message container
struct MessageParticleContainerInfo
{
    explicit MessageParticleContainerInfo(std::string basename);

    ~MessageParticleContainerInfo();

    //! Field name without state information
    std::string m_basename;
};

//! Main SPADES message container
class MessageParticleContainer
    : public amrex::NeighborParticleContainer<RealData::ncomps, IntData::ncomps>
{
public:
    /**
       @brief Class identifier name
       @return class identifier
     **/
    static std::string identifier() { return "particles"; }

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

    //! Initialize particle states (counts and offsets)
    void initialize_state();

    //! Delete particle states (counts and offsets)
    void clear_state();

    //! Update the message counts and offsets
    void update_counts();

    //! Update the message counts
    void count_messages();

    //! Update the message offsets
    void count_offsets();

    /**
       @brief Get the message counts
       @return message counts
     **/
    const amrex::iMultiFab& message_counts() const { return m_message_counts; };

    /**
       @brief Get the message offsets
       @return message offsets
     **/
    const amrex::iMultiFab& offsets() const { return m_offsets; };

    /**
       @brief Get the total number of messages of \p typ
       @param typ [in] message type
       @return message counts
     **/
    amrex::Long total_count(const int typ) const
    {
        BL_PROFILE("spades::MessageParticleContainer::total_count()");
        return m_message_counts.sum(typ);
    }

    /**
       @brief Initialize the messages
       @param lookahead [in] lookahead
     **/
    void initialize_messages(const amrex::Real lookahead);

    //! Sort the messages
    void sort_messages();

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

    //! Reposition the messages inside a cell for visualization
    void reposition_messages();

    /**
       @brief Write the particles to file
       @param plt_filename [in] file name for the plot file
     **/
    void write_plot_file(const std::string& plt_filename);

    //! Number of grow cells
    int ngrow() const { return m_ngrow; }

    //! Level index
    static constexpr int LEV{0};

private:
    //! Initialize vectors (write flags, message counts, offsets, etc)
    void initialize_vectors();

    //! Data for particle container
    MessageParticleContainerInfo m_info;

    //! Flags for real data to write to file
    amrex::Vector<int> m_writeflags_real;

    //! Flags for int data to write to file
    amrex::Vector<int> m_writeflags_int;

    //! Names for real data to write to file
    amrex::Vector<std::string> m_real_data_names;

    //! Names for int data to write to file
    amrex::Vector<std::string> m_int_data_names;

    //! Number of grow cells
    int m_ngrow;

    //! Count of message types in each cell
    amrex::iMultiFab m_message_counts;

    //! Offsets of message types in each cell
    amrex::iMultiFab m_offsets;

    //! Lower bound for undefined message counts
    const int m_lower_undefined_count{40};

    //! Reset value for undefined message counts
    const int m_reset_undefined_count{50};

    //! Upper bound for undefined message counts
    const int m_upper_undefined_count{100};
};
} // namespace spades::particles
#endif
