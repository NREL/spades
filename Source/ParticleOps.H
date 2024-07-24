#ifndef PARTICLEOPS_H
#define PARTICLEOPS_H

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Print.H>
#include "ParticleData.H"

namespace spades::particles {

//! Functor for printing particle data
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
                       << "  id: " << p.id() << "\n"
                       << "  cpu: " << p.cpu() << "\n"
                       << "  type: " << p.idata(IntData::type_id) << "\n"
                       << "  sender: " << p.idata(IntData::sender) << "\n"
                       << "  receiver: " << p.idata(IntData::receiver) << "\n"
                       << "  pair: " << p.idata(IntData::pair) << "\n"
                       << "  timestamp: " << p.rdata(RealData::timestamp)
                       << "\n"
                       << "  old_timestamp: "
                       << p.rdata(RealData::old_timestamp) << "\n"
                       << "  creation_time: "
                       << p.rdata(RealData::creation_time) << std::endl;
    }
};

/**
   @brief Functor for comparing particles

   This is used to sort particles based on cell index, message type,
   and timestamp (in that order).
 **/
struct Compare
{
    /**
       @brief Compare particles
       @param p1 [in] first particle
       @param p2 [in] second particle
       @return boolean of p1 < p2 based on criteria
     **/
    template <class P>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
    operator()(const P& p1, const P& p2) const
    {
        // sort by iv, then message type then by timestamp
        const amrex::IntVect piv1(AMREX_D_DECL(
            p1.idata(IntData::i), p1.idata(IntData::j), p1.idata(IntData::k)));
        const amrex::IntVect piv2(AMREX_D_DECL(
            p2.idata(IntData::i), p2.idata(IntData::j), p2.idata(IntData::k)));

        const auto m1 = p1.idata(IntData::type_id);
        const auto m2 = p2.idata(IntData::type_id);
        const auto t1 = p1.rdata(RealData::timestamp);
        const auto t2 = p2.rdata(RealData::timestamp);
        return (piv1 < piv2) ||
               (piv1 == piv2 && ((m1 < m2) || (m1 == m2 && t1 < t2)));
    }
};

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
        for (int i = 0; i < RealData::ncomps; i++) {
            p2.rdata(i) = p1.rdata(i);
        }
        for (int i = 0; i < IntData::ncomps; i++) {
            p2.idata(i) = p1.idata(i);
        }
    }
};

//! Functor for making a conjugate particle out of another
struct MakeConjugate
{
    /**
       @brief Create a conjugate particle
       @param p1 [in] first particle
       @param p2 [out] second particle
     **/
    template <class P>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(const P& p1, P& p2) const
    {
        Copy()(p1, p2);
        p2.idata(IntData::type_id) = MessageTypes::CONJUGATE;
    }
};

//! Functor for populating particle data
struct Create
{
    /**
       @brief Create a conjugate particle
       @param p [in, out] particle
       @param timestamp [in] time stamp
       @param pos [in] position
       @param iv [in] index
       @param sender [in] sender id
       @param receiver [in] receiver id
     **/
    template <class P>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void operator()(
        P& p,
        const amrex::Real timestamp,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& pos,
        const amrex::IntVect& iv,
        const int sender,
        const int receiver) const
    {
        AMREX_ASSERT(p.idata(IntData::type_id) == MessageTypes::UNDEFINED);

        p.idata(IntData::type_id) = MessageTypes::MESSAGE;
        p.rdata(RealData::timestamp) = timestamp;

        AMREX_D_TERM(p.pos(0) = pos[0];, p.pos(1) = pos[1];, p.pos(2) = pos[2];)
        AMREX_D_TERM(p.idata(IntData::i) = iv[0];, p.idata(IntData::j) = iv[1];
                     , p.idata(IntData::k) = iv[2];)
        p.idata(IntData::sender) = sender;
        p.idata(IntData::receiver) = receiver;
    }
};

//! Functor for making a particle undefined
struct MarkUndefined
{
    /**
       @brief Mark a particle as undefined
       @param p [in, out] particle
     **/
    template <class P>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void operator()(P& p) const
    {
        p.idata(IntData::type_id) = MessageTypes::UNDEFINED;
        p.idata(IntData::pair) = -1;
        p.rdata(RealData::timestamp) = 0.0;
        p.rdata(RealData::old_timestamp) = 0.0;
        p.rdata(RealData::creation_time) = 0.0;
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

        AMREX_ALWAYS_ASSERT(m_counts(m_iv, typ) > n);
        const int idx = m_offsets(m_iv, typ) + n;

        P& p = m_pstruct[idx];

        AMREX_ALWAYS_ASSERT(p.id() >= 0);
        AMREX_ALWAYS_ASSERT(p.idata(IntData::type_id) == typ);

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
        AMREX_ALWAYS_ASSERT(m_counts(m_iv, typ) > n);
        const int idx = m_offsets(m_iv, typ) + n;
        P& p = m_pstruct[idx];
        return p.idata(IntData::type_id) == typ;
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
        const int idx = m_offsets(m_iv, typ) + n;
        P& p = m_pstruct[idx];
        AMREX_ALWAYS_ASSERT(p.idata(IntData::type_id) == expected_type);
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

} // namespace spades::particles
#endif
