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
template <
    int NType,
    int NStructReal,
    int NStructInt,
    int NArrayReal,
    int NArrayInt>
class SpadesParticleContainer
    : public amrex::NeighborParticleContainer<
          NStructReal,
          NStructInt,
          NArrayReal,
          NArrayInt>
{
public:
    using ParticleType = typename amrex::NeighborParticleContainer<
        NStructReal,
        NStructInt,
        NArrayReal,
        NArrayInt>::ParticleType;
    using ParticleTileType = typename amrex::NeighborParticleContainer<
        NStructReal,
        NStructInt,
        NArrayReal,
        NArrayInt>::ParticleTileType;
    using IntVector = typename amrex::NeighborParticleContainer<
        NStructReal,
        NStructInt,
        NArrayReal,
        NArrayInt>::IntVector;
    using RealVector = typename amrex::NeighborParticleContainer<
        NStructReal,
        NStructInt,
        NArrayReal,
        NArrayInt>::RealVector;

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

    /**
       @brief Get the minimum number of particles of \p typ
       @param typ [in] particle type
       @return particle counts min
     **/
    amrex::Long min_count(const int typ) const
    {
        BL_PROFILE("spades::SpadesParticleContainer::min_count()");
        return m_counts.min(typ);
    }

    /**
       @brief Get the maximum number of particles of \p typ
       @param typ [in] particle type
       @return particle counts min
     **/
    amrex::Long max_count(const int typ) const
    {
        BL_PROFILE("spades::SpadesParticleContainer::max_count()");
        return m_counts.max(typ);
    }

    //! Check the result of the sort operation
    void check_sort(const amrex::MFIter& mfi);

    //! Sort the particles
    virtual void sort() = 0;

    /**
      @brief Sort the particles implementation
      @param compare [in] comparison functor
    **/
    template <typename CompareFunctor>
    void sort_impl(const CompareFunctor& compare);

    /**
      @brief Non-encoded sort the particles implementation
      @param compare [in] comparison functor
    **/
    template <typename CompareFunctor>
    void nonencoded_sort_impl(const CompareFunctor& compare);

    //! Encoded sort the particles implementation
    void encoded_sort_impl();

    //! Print all the particles to screen
    void print_messages(const std::string& header);

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

    //! Initialize variable names
    virtual void initialize_variable_names() = 0;

    /**
      @brief Write the particles to file (implementation)
      @param plt_filename [in] file name for the plot file
      @param name [in] file name for the particle file
    **/
    void write_plot_file_impl(
        const std::string& plt_filename, const std::string& name);

    //! Read user parameters
    virtual void read_parameters()
    {
        {
            amrex::ParmParse pp("spades");
            pp.query("sort_type", m_sort_type);
            check_sort_type(m_sort_type);
        }
    }

    //! Check valid sort type
    void check_sort_type(const std::string& sort_type)
    {
        const amrex::Vector<std::string> valid_types = {
            "nonencoded", "encoded"};
        if (std::find(valid_types.cbegin(), valid_types.cend(), sort_type) ==
            valid_types.cend()) {
            amrex::Abort("Invalid sort type. Must be nonencoded or encoded");
        }
    }

    ParticleOps<NArrayReal, NArrayInt, ParticleType, RealVector, IntVector>
    particle_arrays(ParticleTileType& pti) const
    {
        auto& aos = pti.GetArrayOfStructs();
        auto& soa = pti.GetStructOfArrays();
        return ParticleOps<
            NArrayReal, NArrayInt, ParticleType, RealVector, IntVector>(
            aos.dataPtr(), soa.GetRealData(), soa.GetIntData());
    }

protected:
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

    //! Sort type
    std::string m_sort_type{"nonencoded"};

    //! Count of particles in each cell
    amrex::iMultiFab m_counts;

    //! Offsets of particles in each cell
    amrex::iMultiFab m_offsets;

    //! Minimum timestamp in each cell for each type
    amrex::MultiFab m_min_timestamp;

    //! Maximum timestamp in each cell for each type
    amrex::MultiFab m_max_timestamp;
};

#include "SpadesParticleContainerImpl.H"

} // namespace spades::particles
#endif
