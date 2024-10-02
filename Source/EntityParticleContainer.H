#ifndef ENTITYPARTICLECONTAINER_H
#define ENTITYPARTICLECONTAINER_H
#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_NeighborParticles.H>
#include <AMReX_Random.H>
#include "EntityData.H"
#include "EntityOps.H"
#include "Utilities.H"

/**
   @brief SPADES particles
 **/
namespace spades::particles {

//! Main SPADES entity container
class EntityParticleContainer
    : public amrex::NeighborParticleContainer<
          EntityRealData::ncomps,
          EntityIntData::ncomps>
{
public:
    /**
       @brief Class identifier name
       @return class identifier
     **/
    static std::string identifier() { return "entities"; }

    /**
       @brief Constructor
       @param par_gdb [in] particle database
       @param ngrow [in] number of grow cells
     **/
    explicit EntityParticleContainer(amrex::AmrParGDB* par_gdb, int ngrow = 0);

    /**
       @brief Constructor
       @param geom [in] geometry
       @param dmap [in] distribution map
       @param ba [in] box array
       @param ngrow [in] number of grow cells
     **/
    explicit EntityParticleContainer(
        const amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<amrex::DistributionMapping>& dmap,
        const amrex::Vector<amrex::BoxArray>& ba,
        int ngrow = 0);

    //! Initialize particle states (counts and offsets)
    void initialize_state();

    //! Delete particle states (counts and offsets)
    void clear_state();

    //! Update the entity counts and offsets
    void update_counts();

    //! Update the entity counts
    void count_entities();

    //! Update the entity offsets
    void count_offsets();

    /**
       @brief Get the entity counts
       @return entity counts
     **/
    const amrex::iMultiFab& entity_counts() const { return m_entity_counts; };

    /**
       @brief Get the entity offsets
       @return entity offsets
     **/
    const amrex::iMultiFab& offsets() const { return m_offsets; };

    /**
       @brief Get the total number of entities of \p typ
       @param typ [in] entity type
       @return entity counts
     **/
    amrex::Long total_count(const int typ) const
    {
        BL_PROFILE("spades::EntityParticleContainer::total_count()");
        return m_entity_counts.sum(typ);
    }

    //! Initialize the entities
    void initialize_entities();

    //! Sort the entities
    void sort_entities();

    //! Update the undefined entities
    void update_undefined();

    /**
       @brief Perform garbage collection
       @param gvt [in] global virtual time
     **/
    void garbage_collect(const amrex::Real gvt);

    //! Reposition the entities inside a cell for visualization
    void reposition_entities();

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
    //! Initialize vectors (write flags, entity counts, offsets, etc)
    void initialize_vectors();

    //! Data for particle container
    ParticleContainerInfo m_info;

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

    //! Count of entity types in each cell
    amrex::iMultiFab m_entity_counts;

    //! Offsets of entity types in each cell
    amrex::iMultiFab m_offsets;

    //! Lower bound for undefined entity counts
    const int m_lower_undefined_count{40};

    //! Reset value for undefined entity counts
    const int m_reset_undefined_count{50};

    //! Upper bound for undefined entity counts
    const int m_upper_undefined_count{100};
};
} // namespace spades::particles
#endif
