#ifndef CELLSORTEDPARTICLECONTAINER_H
#define CELLSORTEDPARTICLECONTAINER_H
#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_NeighborParticles.H>
#include "ParticleData.H"
#include "ParticleOps.H"
#include "Utilities.H"

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

/** Information for a particle container
 *  \ingroup particles
 */
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
    explicit CellSortedParticleContainer(
        amrex::AmrParGDB* par_gdb, int ngrow = 0);

    void initialize_state();

    void clear_state();

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

    void initialize_particles(const amrex::Real lookahead);

    void add_messages();

    void sort_particles();

    void garbage_collect(const amrex::Real gvt);

    void write_plot_file(const std::string& plt_filename);

    int ngrow() const { return m_ngrow; }

protected:
    amrex::BoxArray::RefID m_BARef;

    amrex::DistributionMapping::RefID m_DMRef;

private:
    //! Data for particle container
    ParticleContainerInfo m_info;

    amrex::Vector<int> m_writeflags_real;
    amrex::Vector<int> m_writeflags_int;

    amrex::Vector<std::string> m_real_data_names;
    amrex::Vector<std::string> m_int_data_names;

    int m_ngrow;

    amrex::Vector<amrex::iMultiFab> m_message_counts;
    amrex::Vector<amrex::iMultiFab> m_offsets;
};
} // namespace spades::particles
#endif
