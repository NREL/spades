#ifndef PARTICLEOPS_H
#define PARTICLEOPS_H

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Print.H>
#include "ParticleData.H"

namespace spades::particles {

//! Particle operations
template <size_t NReal, size_t NInt, class P, class RT, class IT>
struct ParticleArrays
{
    /**
       @brief Constructor
       @param aos [in] particle AOS data
       @param rdata [in] particle real data
       @param idata [in] particle int data
     **/
    ParticleArrays(
        P* aos, std::array<RT, NReal>& rdata, std::array<IT, NInt>& idata)
        : m_aos(aos)
    {
        for (size_t i = 0; i < NReal; i++) {
            m_rdata[i] = rdata[i].data();
        }
        for (size_t i = 0; i < NInt; i++) {
            m_idata[i] = idata[i].data();
        }
    }

    //! AOS data
    P* m_aos = nullptr;

    //! Pointers to real data
    amrex::GpuArray<amrex::Real*, NReal> m_rdata;

    //! Pointers to int data
    amrex::GpuArray<int*, NInt> m_idata;
};

//! Functor for printing particle data on host
struct Print
{
    /**
       @brief Print particle data to screen
       @param n [in] particle index whose data is printed
       @param parrs [in] particle arrays
     **/
    template <class PArrs>
    void operator()(const amrex::Long n, PArrs& parrs) const
    {
        auto& p = parrs.m_aos[n];
        amrex::Print() << "Particle data: \n"
                       << "  id = " << p.id() << "\n"
                       << "  cpu = " << p.cpu() << "\n"
                       << "  type = "
                       << parrs.m_idata[CommonIntData::type_id][n] << "\n";
        for (int i = 0; i < p.NInt; i++) {
            amrex::Print() << "  int comp(" << i
                           << ") = " << parrs.m_idata[i][n] << "\n";
        }
        for (int i = 0; i < p.NReal; i++) {
            amrex::Print() << "  real comp(" << i
                           << ") = " << parrs.m_rdata[i][n] << "\n";
        }
    }
};

//! Functor for printing particle data on device
struct DevicePrint
{
    /**
         @brief Print particle data to screen
         @param n [in] particle index whose data is printed
         @param parrs [in] particle arrays
       **/
    template <class PArrs>
    AMREX_GPU_DEVICE void operator()(const amrex::Long n, PArrs& parrs) const
    {
        auto& p = parrs.m_aos[n];
        printf(
            "Particle data: id = %ld, cpu = %ld, type = %d, timestamp = %.8e\n",
            static_cast<amrex::Long>(p.id()), static_cast<amrex::Long>(p.cpu()),
            parrs.m_idata[CommonIntData::type_id][n],
            parrs.m_rdata[CommonRealData::timestamp][n]);
    }
};

//! Functor for copying a particle's data to another
struct Copy
{
    /**
         @brief Copy particle data
         @param p1 [in] first particle index
         @param p2 [in] second particle index
         @param parrs [in, out] particle arrays
       **/
    template <class PArrs>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(const amrex::Long p1, const amrex::Long p2, PArrs& parrs) const
    {
        for (int i = 0; i < PArrs::NReal; i++) {
            parrs.m_rdata[i][p2] = parrs.m_rdata[i][p1];
        }
        for (int i = 0; i < PArrs::NInt; i++) {
            parrs.m_idata[i][p2] = parrs.m_idata[i][p1];
        }
    }
};

//! Functor for accessing a particle in a cell
template <class PArrs>
struct Get
{
    /**
       @brief Constructor
       @param iv [in] cell index of the particle
       @param counts [in] cell counts of particle types
       @param offsets [in] cell offsets of particle types
       @param parrs [in] particle arrays
     **/
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    Get(const amrex::IntVect& iv,
        const amrex::Array4<const int>& counts,
        const amrex::Array4<const int>& offsets,
        PArrs& parrs)
        : m_iv(iv), m_counts(counts), m_offsets(offsets), m_parrs(parrs)
    {}

    /**
       @brief Get a particle in a cell
       @param n [in] nth particle to return
       @param typ [in] particle type that is requested
       @return particle index of type \p typ and at position \p n
     **/
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
    operator()(const int n, const int typ) const
    {

        AMREX_ASSERT(m_counts(m_iv, typ) > n);
        const int idx = m_offsets(m_iv, typ) + n;

        AMREX_ASSERT(m_parrs.m_aos[idx].id() >= 0);
        AMREX_ASSERT(m_parrs.m_idata[CommonIntData::type_id][idx] == typ);

        return idx;
    }

    /**
       @brief Check validity of the particle type
       @param n [in] nth particle to return
       @param typ [in] particle type that is requested
       @return boolean if the particle is of the requested type
     **/
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
    check(const int n, const int typ) const
    {
        AMREX_ASSERT(m_counts(m_iv, typ) > n);
        const int idx = m_offsets(m_iv, typ) + n;
        return m_parrs.m_idata[CommonIntData::type_id][idx] == typ;
    }

    /**
       @brief Assert that the particle is of a different but expected
       type

       A particle sometimes changes type (e.g., marked as undefined)
       but the counts and offsets haven't been updated yet. This
       asserts that the requested particle of a certain type \p typ
       changed to another type \p expected_type.

       @param n [in] nth particle to return
       @param typ [in] particle type that is requested
       @param expected_type [in] expected particle type
     **/
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    assert_different(const int n, const int typ, const int expected_type) const
    {
#ifdef AMREX_DEBUG
        const int idx = m_offsets(m_iv, typ) + n;
        AMREX_ASSERT(
            m_parrs.m_idata[CommonIntData::type_id][idx] == expected_type);
#else
        amrex::ignore_unused(n, typ, expected_type);
#endif
    }

    //! Cell index
    const amrex::IntVect& m_iv;

    //! Cell counts of particle types
    const amrex::Array4<const int>& m_counts;

    //! Cell offsets of particle types
    const amrex::Array4<const int>& m_offsets;

    //! Particle array
    PArrs& m_parrs;
};

//! Functor for pairing particles
struct CheckPair
{
    /**
       @brief Compare particles for pairing
       @return boolean of p1 and p2 are pairs
     **/
    template <class P>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
    operator()(const P& /*p1*/, const P& /*p2*/) const
    {
        return false;
    }
};

/**
   @brief Functor for comparing particles

   This is used to sort particles based on cell index, entity type,
   and timestamp (in that order).
 **/
struct CompareParticle
{
    /**
     @brief Compare particles
     @param p1 [in] first particle
     @param p2 [in] second particle
     @param parrs [in] particle arrays
     @return boolean of p1 < p2 based on criteria
   **/
    template <class PArrs>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
    operator()(const amrex::Long p1, const amrex::Long p2, PArrs& parrs) const
    {
        // sort by iv, then particle type then by timestamp
        const amrex::IntVect piv1(AMREX_D_DECL(
            parrs.m_idata[CommonIntData::i][p1],
            parrs.m_idata[CommonIntData::j][p1],
            parrs.m_idata[CommonIntData::k][p1]));
        const amrex::IntVect piv2(AMREX_D_DECL(
            parrs.m_idata[CommonIntData::i][p2],
            parrs.m_idata[CommonIntData::j][p2],
            parrs.m_idata[CommonIntData::k][p2]));

        const auto m1 = parrs.m_idata[CommonIntData::type_id][p1];
        const auto m2 = parrs.m_idata[CommonIntData::type_id][p2];
        const auto t1 = parrs.m_rdata[CommonRealData::timestamp][p1];
        const auto t2 = parrs.m_rdata[CommonRealData::timestamp][p2];
        return (piv1 < piv2) ||
               (piv1 == piv2 && ((m1 < m2) || (m1 == m2 && t1 < t2)));
    }
};

} // namespace spades::particles
#endif
