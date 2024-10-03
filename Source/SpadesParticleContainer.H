#ifndef SPADESPARTICLECONTAINER_H
#define SPADESPARTICLECONTAINER_H
#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_NeighborParticles.H>
#include <AMReX_Random.H>
#include "ParticleData.H"
#include "ParticleOps.H"
#include "Utilities.H"

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#endif

/**
   @brief SPADES particles
 **/
namespace spades::particles {

//! Main SPADES particle container
template <typename PType, int NStructReal, int NStructInt>
class SpadesParticleContainer
    : public amrex::NeighborParticleContainer<NStructReal, NStructInt>
{
public:
    /**
       @brief Class identifier name
       @return class identifier
     **/
    static std::string identifier() { return "spades_particles"; }

    /**
       @brief Constructor
       @param par_gdb [in] particle database
       @param ngrow [in] number of grow cells
     **/
    explicit SpadesParticleContainer(amrex::AmrParGDB* par_gdb, int ngrow = 0);

    /**
       @brief Constructor
       @param geom [in] geometry
       @param dmap [in] distribution map
       @param ba [in] box array
       @param ngrow [in] number of grow cells
     **/
    explicit SpadesParticleContainer(
        const amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<amrex::DistributionMapping>& dmap,
        const amrex::Vector<amrex::BoxArray>& ba,
        int ngrow = 0);

    //! Initialize particle states (counts and offsets)
    void initialize_state();

    //! Delete particle states (counts and offsets)
    void clear_state();

    //! Update the particle counts and offsets
    void update_counts();

    //! Update the particle counts
    void count_particles();

    //! Update the particle offsets
    void count_offsets();

    /**
       @brief Get the particle counts
       @return particle counts
     **/
    const amrex::iMultiFab& counts() const { return m_counts; };

    /**
       @brief Get the particle offsets
       @return particle offsets
     **/
    const amrex::iMultiFab& offsets() const { return m_offsets; };

    /**
       @brief Get the total number of particles of \p typ
       @param typ [in] particle type
       @return particle counts
     **/
    amrex::Long total_count(const int typ) const
    {
        BL_PROFILE("spades::SpadesParticleContainer::total_count()");
        return m_counts.sum(typ);
    }

    //! Sort the particles
    virtual void sort() = 0;

    //! Reposition the particles inside a cell for visualization
    void reposition_particles();

    /**
      @brief Write the particles to file
      @param plt_filename [in] file name for the plot file
    **/
    virtual void write_plot_file(const std::string& plt_filename) = 0;

    //! Number of grow cells
    int ngrow() const { return m_ngrow; }

    //! Level index
    static constexpr int LEV{0};

protected:
    //! Initialize variable names
    virtual void initialize_variable_names() = 0;

    /**
      @brief Sort the particles implementation
      @param compare [in] comparison functor
    **/
    //!
    template <typename CompareFunctor>
    void sort_impl(const CompareFunctor& compare);

    /**
      @brief Write the particles to file (implementation)
      @param plt_filename [in] file name for the plot file
      @param name [in] file name for the particle file
    **/
    void write_plot_file_impl(
        const std::string& plt_filename, const std::string& name);

    //! Number of grow cells
    int m_ngrow;

    //! Flags for real data to write to file
    amrex::Vector<int> m_writeflags_real;

    //! Flags for int data to write to file
    amrex::Vector<int> m_writeflags_int;

    //! Names for real data to write to file
    amrex::Vector<std::string> m_real_data_names;

    //! Names for int data to write to file
    amrex::Vector<std::string> m_int_data_names;

    //! Count of types in each cell
    amrex::iMultiFab m_counts;

    //! Offsets of types in each cell
    amrex::iMultiFab m_offsets;
};

#include "SpadesParticleContainerImpl.H"

} // namespace spades::particles
#endif
