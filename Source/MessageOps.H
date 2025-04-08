#ifndef MESSAGEOPS_H
#define MESSAGEOPS_H

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Print.H>
#include "MessageData.H"
#include "ParticleOps.H"

namespace spades::particles {

//! Functor for populating message data
struct CreateMessage
{
    /**
       @brief Create a message
       @param n [in] message index
       @param parrs [in, out] message arrays
       @param timestamp [in] time stamp
       @param pos [in] position
       @param iv [in] index
       @param sender_lp [in] sender LP id
       @param sender_entity [in] entity id
       @param receiver_lp [in] receiver LP id
       @param receiver_entity [in] entity id
     **/
    template <class PArrs>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void operator()(
        const amrex::Long n,
        PArrs& parrs,
        const amrex::Real timestamp,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& pos,
        const amrex::IntVect& iv,
        const int sender_lp,
        const int sender_entity,
        const int receiver_lp,
        const int receiver_entity) const
    {
        AMREX_ASSERT(
            parrs.m_idata[MessageIntData::type_id][n] ==
            MessageTypes::UNDEFINED);

        parrs.m_idata[MessageIntData::type_id][n] = MessageTypes::MESSAGE;
        parrs.m_rdata[MessageRealData::timestamp][n] = timestamp;

        auto& p = parrs.m_aos[n];
        AMREX_D_TERM(p.pos(0) = pos[0];, p.pos(1) = pos[1];, p.pos(2) = pos[2];)
        AMREX_D_TERM(parrs.m_idata[MessageIntData::i][n] = iv[0];
                     , parrs.m_idata[MessageIntData::j][n] = iv[1];
                     , parrs.m_idata[MessageIntData::k][n] = iv[2];)
        parrs.m_idata[MessageIntData::sender_lp][n] = sender_lp;
        parrs.m_idata[MessageIntData::sender_entity][n] = sender_entity;
        parrs.m_idata[MessageIntData::receiver_lp][n] = receiver_lp;
        parrs.m_idata[MessageIntData::receiver_entity][n] = receiver_entity;
    }
};

//! Functor for making a message undefined
struct MarkMessageUndefined
{
    /**
       @brief Mark a message as undefined
       @param n [in] message index
       @param parrs [in, out] message arrays
     **/
    template <class PArrs>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    operator()(const amrex::Long n, PArrs& parrs) const
    {
        parrs.m_idata[MessageIntData::type_id][n] = MessageTypes::UNDEFINED;
        parrs.m_idata[MessageIntData::pair][n] = -1;
        parrs.m_rdata[MessageRealData::timestamp][n] = 0.0;
        parrs.m_rdata[MessageRealData::old_timestamp][n] = 0.0;
        parrs.m_rdata[MessageRealData::creation_time][n] = 0.0;
    }
};

} // namespace spades::particles
#endif
