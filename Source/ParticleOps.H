#ifndef PARTICLEOPS_H
#define PARTICLEOPS_H

#include <AMReX.H>
#include <AMReX_Array4.H>
#include "ParticleData.H"

namespace spades::particles {

struct Compare
{
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

struct Copy
{
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

struct MakeConjugate
{
    template <class P>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(const P& p1, P& p2) const
    {
        Copy()(p1, p2);
        p2.idata(IntData::type_id) = MessageTypes::CONJUGATE;
    }
};

struct Create
{
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

struct MarkUndefined
{
    template <class P>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void operator()(P& p) const
    {
        p.idata(IntData::type_id) = MessageTypes::UNDEFINED;
        p.idata(IntData::pair) = -1;
    }
};

template <class P>
struct Get
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    Get(const amrex::IntVect& iv,
        const amrex::Array4<const int>& counts,
        const amrex::Array4<const int>& offsets,
        P* pstruct)
        : m_iv(iv), m_counts(counts), m_offsets(offsets), m_pstruct(pstruct)
    {}

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE P&
    operator()(const int n, const int typ) const
    {

        AMREX_ALWAYS_ASSERT(m_counts(m_iv, typ) > n);
        const int idx = m_offsets(m_iv, typ) + n;

        P& p = m_pstruct[idx];

        AMREX_ALWAYS_ASSERT(p.idata(IntData::type_id) == typ);

        return p;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
    check(const int n, const int typ) const
    {
        AMREX_ALWAYS_ASSERT(m_counts(m_iv, typ) > n);
        const int idx = m_offsets(m_iv, typ) + n;
        P& p = m_pstruct[idx];
        return p.idata(IntData::type_id) == typ;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    assert_different(const int n, const int typ, const int expected_type) const
    {
        const int idx = m_offsets(m_iv, typ) + n;
        P& p = m_pstruct[idx];
        AMREX_ALWAYS_ASSERT(p.idata(IntData::type_id) == expected_type);
    }

    const amrex::IntVect& m_iv;
    const amrex::Array4<const int>& m_counts;
    const amrex::Array4<const int>& m_offsets;
    P* m_pstruct;
};

} // namespace spades::particles
#endif
