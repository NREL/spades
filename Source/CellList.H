#ifndef CELLLIST_H
#define CELLLIST_H

#include "AMReX_Gpu.H"
#include "ParticleOps.H"
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#endif

namespace spades::particles {

template <class ParticleType>
class CellList
{
public:
    template <class PTile>
    void build(PTile& pti, const amrex::Box& box)
    {
        BL_PROFILE("spades::CellList::build");

        // Note: I think another way to do all this would be to use
        // `ReorderParticles` and `SortParticles` in
        // AMReX_ParticleContainer.H. Then they would be sorted in
        // memory and things would be more efficient

        const auto& particles = pti.GetArrayOfStructs();
        const size_t np = pti.numParticles();
        auto* pstruct = particles().dataPtr();

        // Get the particle indices and their cell index
        // BL_PROFILE_VAR("spades::CellList::prep", prep);
        m_cell_list.resize(np);
        amrex::Gpu::DeviceVector<amrex::Long> cell_index(np, 0);
        auto* p_cell_index = cell_index.data();
        auto* p_cell_list = m_cell_list.data();
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pindex) noexcept {
            const auto& p = pstruct[pindex];
            const amrex::IntVect piv(AMREX_D_DECL(
                p.idata(particles::IntData::i), p.idata(particles::IntData::j),
                p.idata(particles::IntData::k)));
            p_cell_list[pindex] = pindex;
            p_cell_index[pindex] = box.index(piv);
        });
        // amrex::Gpu::Device::synchronize();
        // BL_PROFILE_VAR_STOP(prep);

        // Sort particle indices based on the cell index
        // BL_PROFILE_VAR("spades::CellList::sort", sort);
#ifdef AMREX_USE_GPU
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
        thrust::sort(
            thrust::device, m_cell_list.begin(), m_cell_list.end(),
            [=] AMREX_GPU_DEVICE(
                const amrex::Long x, const amrex::Long y) noexcept {
                const auto& p1 = pstruct[x];
                const auto& p2 = pstruct[y];
                return Compare()(p1, p2);
            });
#else
        // Perform sort on CPU, then copy back to device (not good)
        amrex::Vector<amrex::Long> h_cell_list(np, 0);
        amrex::Gpu::copy(
            amrex::Gpu::deviceToHost, m_cell_list.begin(), m_cell_list.end(),
            h_cell_list.begin());
        std::sort(
            h_cell_list.begin(), h_cell_list.end(),
            [=](const amrex::Long x, const amrex::Long y) {
                const auto& p1 = pstruct[x];
                const auto& p2 = pstruct[y];
                return Compare()(p1, p2);
            });
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, h_cell_list.begin(), h_cell_list.end(),
            m_cell_list.begin());
#endif
#else
        std::sort(
            m_cell_list.begin(), m_cell_list.end(),
            [=](const amrex::Long x, const amrex::Long y) {
                const auto& p1 = pstruct[x];
                const auto& p2 = pstruct[y];
                return Compare()(p1, p2);
            });
#endif
        // amrex::Gpu::Device::synchronize();
        // BL_PROFILE_VAR_STOP(sort);

        // Count the particles
        // BL_PROFILE_VAR("spades::CellList::count", count);
        const auto npts = box.numPts();
        m_cell_counts.resize(npts);
        m_cell_counts.assign(npts, 0);
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
        amrex::Gpu::DeviceVector<amrex::Long> ones(np, 1);
        amrex::Gpu::DeviceVector<amrex::Long> okeys(np, 0);
        thrust::reduce_by_key(
            thrust::device, cell_index.begin(), cell_index.end(), ones.begin(),
            okeys.begin(), m_cell_counts.begin());
#else
        auto* p_cell_counts = m_cell_counts.data();
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pindex) noexcept {
            amrex::Gpu::Atomic::AddNoRet(
                &p_cell_counts[p_cell_index[pindex]],
                static_cast<amrex::Long>(1));
        });
#endif
        // amrex::Gpu::Device::synchronize();
        // BL_PROFILE_VAR_STOP(count);

        // Compute the offsets
        // BL_PROFILE_VAR("spades::CellList::offset", offset);
        m_cell_offsets.resize(npts);
        amrex::Gpu::exclusive_scan(
            m_cell_counts.begin(), m_cell_counts.end(), m_cell_offsets.begin());
        // amrex::Gpu::Device::synchronize();
        // BL_PROFILE_VAR_STOP(offset);
    }

    const amrex::Gpu::DeviceVector<amrex::Long>& list() const
    {
        return m_cell_list;
    }
    const amrex::Gpu::DeviceVector<amrex::Long>& counts() const
    {
        return m_cell_counts;
    }
    const amrex::Gpu::DeviceVector<amrex::Long>& offsets() const
    {
        return m_cell_offsets;
    }

protected:
    amrex::Gpu::DeviceVector<amrex::Long> m_cell_offsets;
    amrex::Gpu::DeviceVector<amrex::Long> m_cell_list;
    amrex::Gpu::DeviceVector<amrex::Long> m_cell_counts;
};
} // namespace spades::particles
#endif
