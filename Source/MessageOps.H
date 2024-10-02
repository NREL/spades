#ifndef MESSAGEOPS_H
#define MESSAGEOPS_H

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Print.H>
#include "MessageData.H"
#include "ParticleOps.H"

namespace spades::particles {

/**
   @brief Functor for comparing messages

   This is used to sort messages based on cell index, message type,
   and timestamp (in that order).
 **/
struct CompareMessage
{
    /**
       @brief Compare messages
       @param p1 [in] first message
       @param p2 [in] second message
       @return boolean of p1 < p2 based on criteria
     **/
    template <class P>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
    operator()(const P& p1, const P& p2) const
    {
        // sort by iv, then message type then by timestamp
        const amrex::IntVect piv1(AMREX_D_DECL(
            p1.idata(MessageIntData::i), p1.idata(MessageIntData::j),
            p1.idata(MessageIntData::k)));
        const amrex::IntVect piv2(AMREX_D_DECL(
            p2.idata(MessageIntData::i), p2.idata(MessageIntData::j),
            p2.idata(MessageIntData::k)));
        const auto m1 = p1.idata(MessageIntData::type_id);
        const auto m2 = p2.idata(MessageIntData::type_id);
        const auto t1 = p1.rdata(MessageRealData::timestamp);
        const auto t2 = p2.rdata(MessageRealData::timestamp);
        return (piv1 < piv2) ||
               (piv1 == piv2 && ((m1 < m2) || (m1 == m2 && t1 < t2)));
    }
};

//! Functor for making a conjugate message out of another
struct MakeConjugate
{
    /**
       @brief Create a conjugate message
       @param p1 [in] first message
       @param p2 [out] second message
     **/
    template <class P>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(const P& p1, P& p2) const
    {
        Copy()(p1, p2);
        p2.idata(MessageIntData::type_id) = MessageTypes::CONJUGATE;
    }
};

//! Functor for populating message data
struct CreateMessage
{
    /**
       @brief Create a message
       @param p [in, out] message
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
        AMREX_ASSERT(
            p.idata(MessageIntData::type_id) == MessageTypes::UNDEFINED);

        p.idata(MessageIntData::type_id) = MessageTypes::MESSAGE;
        p.rdata(MessageRealData::timestamp) = timestamp;

        AMREX_D_TERM(p.pos(0) = pos[0];, p.pos(1) = pos[1];, p.pos(2) = pos[2];)
        AMREX_D_TERM(p.idata(MessageIntData::i) = iv[0];
                     , p.idata(MessageIntData::j) = iv[1];
                     , p.idata(MessageIntData::k) = iv[2];)
        p.idata(MessageIntData::sender) = sender;
        p.idata(MessageIntData::receiver) = receiver;
    }
};

//! Functor for making a message undefined
struct MarkMessageUndefined
{
    /**
       @brief Mark a message as undefined
       @param p [in, out] message
     **/
    template <class P>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void operator()(P& p) const
    {
        p.idata(MessageIntData::type_id) = MessageTypes::UNDEFINED;
        p.idata(MessageIntData::pair) = -1;
        p.rdata(MessageRealData::timestamp) = 0.0;
        p.rdata(MessageRealData::old_timestamp) = 0.0;
        p.rdata(MessageRealData::creation_time) = 0.0;
    }
};

} // namespace spades::particles
#endif
