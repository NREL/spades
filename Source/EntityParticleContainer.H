#ifndef ENTITYPARTICLECONTAINER_H
#define ENTITYPARTICLECONTAINER_H
#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_Random.H>
#include "SpadesParticleContainer.H"
#include "EntityData.H"
#include "EntityOps.H"
#include "Utilities.H"

/**
   @brief SPADES particles
 **/
namespace spades::particles {

//! Main SPADES entity container
class EntityParticleContainer
    : public SpadesParticleContainer<
          EntityTypes,
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

    //! Initialize the entities
    void initialize_entities();

    //! Sort the entities
    void sort() override;

    /**
       @brief Perform garbage collection
       @param gvt [in] global virtual time
     **/
    void garbage_collect(const amrex::Real gvt);

protected:
    //! Initialize variable names
    void initialize_variable_names() override;

private:
    //! Lower bound for undefined entity counts
    const int m_lower_undefined_count{40};

    //! Reset value for undefined entity counts
    const int m_reset_undefined_count{50};

    //! Upper bound for undefined entity counts
    const int m_upper_undefined_count{100};
};
} // namespace spades::particles
#endif
