#ifndef PARTICLEOPS_H
#define PARTICLEOPS_H

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Print.H>
#include "ParticleData.H"

namespace spades::particles {

//! Functor for printing particle data on host
struct Print
{
    /**
       @brief Print particle data to screen
       @param p [in] particle whose data is printed
     **/
    template <class P>
    void operator()(P& p) const
    {
        amrex::Print() << "Particle data: \n"
                       << "  id = " << p.id() << "\n"
                       << "  cpu = " << p.cpu() << "\n"
                       << "  type = " << p.idata(CommonIntData::type_id)
                       << "\n";
        for (int i = 0; i < p.NInt; i++) {
            amrex::Print() << "  int comp(" << i << ") = " << p.idata(i)
                           << "\n";
        }
        for (int i = 0; i < p.NReal; i++) {
            amrex::Print() << "  real comp(" << i << ") = " << p.rdata(i)
                           << "\n";
        }
    }
};

//! Functor for printing particle data on device
#ifdef AMREX_USE_GPU
struct DevicePrint
{
    /**
       @brief Print particle data to screen
       @param p [in] particle whose data is printed
     **/
    template <class P>
    AMREX_GPU_DEVICE void operator()(P& p) const
    {
      printf("Particle data: id = %ld, cpu = %ld, type = %d, timestamp = %.8e\n", static_cast<amrex::Long>(p.id()), static_cast<amrex::Long>(p.cpu()), p.idata(CommonIntData::type_id), p.rdata(CommonRealData::timestamp));
    }
};
#endif

//! Functor for copying a particle's data to another
struct Copy
{
    /**
       @brief Copy particle data
       @param p1 [in] first particle
       @param p2 [out] second particle
     **/
    template <class P>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(const P& p1, P& p2) const
    {
        AMREX_ASSERT(p1.NReal == p2.NReal);
        AMREX_ASSERT(p1.NInt == p2.NInt);
        for (int i = 0; i < p1.NReal; i++) {
            p2.rdata(i) = p1.rdata(i);
        }
        for (int i = 0; i < p1.NInt; i++) {
            p2.idata(i) = p1.idata(i);
        }
    }
};

//! Functor for accessing a particle in a cell
template <class P>
struct Get
{
    /**
       @brief Constructor
       @param iv [in] cell index of the particle
       @param counts [in] cell counts of particle types
       @param offsets [in] cell offsets of particle types
       @param pstruct [in] particle array
     **/
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    Get(const amrex::IntVect& iv,
        const amrex::Array4<const int>& counts,
        const amrex::Array4<const int>& offsets,
        P* pstruct)
        : m_iv(iv), m_counts(counts), m_offsets(offsets), m_pstruct(pstruct)
    {}

    /**
       @brief Get a particle in a cell
       @param n [in] nth particle to return
       @param typ [in] particle type that is requested
       @return particle of type \p typ and at position \p n
     **/
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE P&
    operator()(const int n, const int typ) const
    {

        AMREX_ASSERT(m_counts(m_iv, typ) > n);
        const int idx = m_offsets(m_iv, typ) + n;

        P& p = m_pstruct[idx];

        AMREX_ASSERT(p.id() >= 0);
        AMREX_ASSERT(p.idata(CommonIntData::type_id) == typ);

        return p;
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
        P& p = m_pstruct[idx];
        return p.idata(CommonIntData::type_id) == typ;
    }

    /**
       @brief Assert that the particle is of a different but expected type

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
        P& p = m_pstruct[idx];
        AMREX_ASSERT(p.idata(CommonIntData::type_id) == expected_type);
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
    P* m_pstruct;
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

} // namespace spades::particles
#endif
