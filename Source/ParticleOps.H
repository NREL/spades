#ifndef PARTICLEOPS_H
#define PARTICLEOPS_H

#include <AMReX.H>
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
            p1.idata(particles::IntData::i), p1.idata(particles::IntData::j),
            p1.idata(particles::IntData::k)));
        const amrex::IntVect piv2(AMREX_D_DECL(
            p2.idata(particles::IntData::i), p2.idata(particles::IntData::j),
            p2.idata(particles::IntData::k)));

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

struct MakeAntiMessage
{
    template <class P>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(const P& p1, P& p2) const
    {
        Copy()(p1, p2);
        p2.idata(IntData::type_id) = MessageTypes::anti_message;
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
        AMREX_ASSERT(p.idata(IntData::type_id) == MessageTypes::undefined);

        p.idata(IntData::type_id) = MessageTypes::message;
        p.rdata(RealData::timestamp) = timestamp;

        AMREX_D_TERM(p.pos(0) = pos[0];, p.pos(1) = pos[1];, p.pos(2) = pos[2];)
        AMREX_D_TERM(p.idata(IntData::i) = iv[0];, p.idata(IntData::j) = iv[1];
                     , p.idata(IntData::k) = iv[2];)
        p.idata(IntData::sender) = sender;
        p.idata(IntData::receiver) = receiver;
    }
};

} // namespace spades::particles
#endif
