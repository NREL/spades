#ifndef PARTICLEOPS_H
#define PARTICLEOPS_H

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Print.H>
#include "ParticleData.H"

namespace spades::particles {

//! Object for storing particle data pointers
template <size_t NReal, size_t NInt, class P, class RT, class IT>
struct ParticleData
{
    ParticleData(
        P* aos, std::array<RT, NReal>& rdata, std::array<IT, NInt>& idata)
        : m_aos(aos)
    {
        for (int i = 0; i < NReal; i++) {
            m_rdata[i] = rdata[i].data();
        }
        for (int i = 0; i < NInt; i++) {
            m_idata[i] = idata[i].data();
        }
    }

    P* m_aos = nullptr;
    amrex::GpuArray<amrex::Real*, NReal> m_rdata;
    amrex::GpuArray<int*, NInt> m_idata;
};

//! Functor for printing particle data on host
struct Print
{
    /**
       @brief Print particle data to screen
       @param p [in] particle whose data is printed
       @param n [in] particle index whose data is printed
       @param rdata [in] particle real data
       @param idata [in] particle int data
     **/
    template <class P, class RT, class IT>
    void operator()(P& p, const amrex::Long n, RT& rdata, IT& idata) const
    {
        amrex::Print() << "Particle data: \n"
                       << "  id = " << p.id() << "\n"
                       << "  cpu = " << p.cpu() << "\n"
                       << "  type = " << idata[CommonIntData::type_id][n]
                       << "\n";
        for (int i = 0; i < p.NInt; i++) {
            amrex::Print() << "  int comp(" << i << ") = " << idata[i][n]
                           << "\n";
        }
        for (int i = 0; i < p.NReal; i++) {
            amrex::Print() << "  real comp(" << i << ") = " << rdata[i][n]
                           << "\n";
        }
    }
};

//! Functor for printing particle data on device
struct DevicePrint
{
    /**
         @brief Print particle data to screen
         @param p [in] particle whose data is printed
         @param n [in] particle index whose data is printed
         @param rdata [in] particle real data
         @param idata [in] particle int data
       **/
    template <class P, class RT, class IT>
    AMREX_GPU_DEVICE void
    operator()(P& p, const amrex::Long n, RT& rdata, IT& idata) const
    {
        printf(
            "Particle data: id = %ld, cpu = %ld, type = %d, timestamp = %.8e\n",
            static_cast<amrex::Long>(p.id()), static_cast<amrex::Long>(p.cpu()),
            idata[CommonIntData::type_id][n], rdata[CommonRealData::timestamp]);
    }
};

//! Functor for copying a particle's data to another
struct Copy
{
    /**
         @brief Copy particle data
         @param p1 [in] first particle index
         @param p2 [in] second particle index
         @param nr [in] number of real components
         @param ni [in] number of int components
         @param rdata [in, out] particle real data
         @param idata [in, out] particle int data
       **/
    template <class RT, class IT>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void operator()(
        const amrex::Long p1,
        const amrex::Long p2,
        const int nr,
        const int ni,
        RT& rdata,
        IT& idata) const
    {
        for (int i = 0; i < nr; i++) {
            rdata[i][p2] = rdata[i][p1];
        }
        for (int i = 0; i < ni; i++) {
            idata[i][p2] = idata[i][p1];
        }
    }
};

//! Functor for accessing a particle in a cell
template <class P, class IT>
struct Get
{
    /**
       @brief Constructor
       @param iv [in] cell index of the particle
       @param counts [in] cell counts of particle types
       @param offsets [in] cell offsets of particle types
       @param pstruct [in] particle array
       @param idata [in] particle int data
     **/
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    Get(const amrex::IntVect& iv,
        const amrex::Array4<const int>& counts,
        const amrex::Array4<const int>& offsets,
        P* pstruct,
        IT& idata)
        : m_iv(iv)
        , m_counts(counts)
        , m_offsets(offsets)
        , m_pstruct(pstruct)
        , m_idata(idata)
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

        AMREX_ASSERT(m_pstruct[idx].id() >= 0);
        AMREX_ASSERT(m_idata[CommonIntData::type_id][idx] == typ);

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
        return m_idata[CommonIntData::type_id][idx] == typ;
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
        AMREX_ASSERT(m_idata[CommonIntData::type_id][idx] == expected_type);
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

    //! Int data
    IT& m_idata;
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
     @param rdata [in] particle real data
     @param idata [in] particle int data
     @return boolean of p1 < p2 based on criteria
   **/
    template <class RT, class IT>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool operator()(
        const amrex::Long p1, const amrex::Long p2, RT& rdata, IT& idata) const
    {
        // sort by iv, then particle type then by timestamp
        const amrex::IntVect piv1(AMREX_D_DECL(
            idata[CommonIntData::i][p1], idata[CommonIntData::j][p1],
            idata[CommonIntData::k][p1]));
        const amrex::IntVect piv2(AMREX_D_DECL(
            idata[CommonIntData::i][p2], idata[CommonIntData::j][p2],
            idata[CommonIntData::k][p2]));

        const auto m1 = idata[CommonIntData::type_id][p1];
        const auto m2 = idata[CommonIntData::type_id][p2];
        const auto t1 = rdata[CommonRealData::timestamp][p1];
        const auto t2 = rdata[CommonRealData::timestamp][p2];
        return (piv1 < piv2) ||
               (piv1 == piv2 && ((m1 < m2) || (m1 == m2 && t1 < t2)));
    }
};

} // namespace spades::particles
#endif
