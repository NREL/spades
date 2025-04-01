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
          0,
          0,
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

    void write_plot_file(const std::string& plt_filename) override
    {
        write_plot_file_impl(plt_filename, identifier());
    };

    //! Read user parameters
    void read_parameters() override;

    //! Initialize variable names
    void initialize_variable_names() override;

private:
    //! Initial entities per LP
    int m_entities_per_lp{1};
};
} // namespace spades::particles
#endif
