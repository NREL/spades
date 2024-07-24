#ifndef CELLSORTEDPARTICLECONTAINER_H
#define CELLSORTEDPARTICLECONTAINER_H
#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_NeighborParticles.H>
#include <AMReX_Random.H>
#include "ParticleData.H"
#include "ParticleOps.H"
#include "Utilities.H"

/**
   @brief SPADES particles
 **/
namespace spades::particles {

struct CheckPair
{
    template <class P>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
    operator()(const P& /*p1*/, const P& /*p2*/) const
    {
        return false;
    }
};

//! Information for a particle container
struct ParticleContainerInfo
{
    explicit ParticleContainerInfo(std::string basename);

    ~ParticleContainerInfo();

    //! Field name without state information
    std::string m_basename;
};

class CellSortedParticleContainer
    : public amrex::NeighborParticleContainer<RealData::ncomps, IntData::ncomps>
{
public:
    static std::string identifier() { return "particles"; }

    explicit CellSortedParticleContainer(
        amrex::AmrParGDB* par_gdb, int ngrow = 0);

    explicit CellSortedParticleContainer(
        const amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<amrex::DistributionMapping>& dmap,
        const amrex::Vector<amrex::BoxArray>& ba,
        int ngrow = 0);

    void initialize_state();

    void clear_state();

    // count messages and offsets
    void update_counts();

    // count messages
    void count_messages();

    // count offsets
    void count_offsets();

    const amrex::Vector<amrex::iMultiFab>& message_counts() const
    {
        return m_message_counts;
    };

    const amrex::Vector<amrex::iMultiFab>& offsets() const
    {
        return m_offsets;
    };

    const amrex::iMultiFab& message_counts(const int lev) const
    {
        return m_message_counts[lev];
    };

    const amrex::iMultiFab& offsets(const int lev) const
    {
        return m_offsets[lev];
    };

    amrex::Long total_count(const int lev, const int typ) const
    {
        BL_PROFILE("spades::CellSortedParticleContainer::total_count()");
        return m_message_counts[lev].sum(typ);
    }

    void initialize_messages(const amrex::Real lookahead);

    void sort_messages();

    // update the undefined messages
    void update_undefined();

    void resolve_pairs();

    amrex::Real gvt();

    void garbage_collect(const amrex::Real gvt);

    // reposition the messages inside a cell
    void reposition_messages();

    void write_plot_file(const std::string& plt_filename);

    int ngrow() const { return m_ngrow; }

private:
    void initialize_vectors();

    //! Data for particle container
    ParticleContainerInfo m_info;

    amrex::Vector<int> m_writeflags_real;
    amrex::Vector<int> m_writeflags_int;

    amrex::Vector<std::string> m_real_data_names;
    amrex::Vector<std::string> m_int_data_names;

    int m_ngrow;

    amrex::Vector<amrex::iMultiFab> m_message_counts;
    amrex::Vector<amrex::iMultiFab> m_offsets;

    const int m_lower_undefined_count{40};
    const int m_reset_undefined_count{50};
    const int m_upper_undefined_count{100};
};
} // namespace spades::particles
#endif
