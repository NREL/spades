#ifndef MESSAGEOPS_H
#define MESSAGEOPS_H

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Print.H>
#include "MessageData.H"
#include "ParticleOps.H"

namespace spades::particles {

//! Functor for making a conjugate message out of another
struct MakeConjugate
{
    /**
         @brief Create a conjugate message
         @param p1 [in] first message index
         @param p2 [out] second message index
         @param nr [in] number of real components
         @param ni [in] number of int components
         @param rdata [in, out] message real data
         @param idata [in, out] message int data
       **/
    template <class RT, class IT>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void operator()(
        const amrex::Long p1,
        amrex::Long p2,
        const int nr,
        const int ni,
        RT& rdata,
        IT& idata) const
    {
        Copy()(p1, p2, nr, ni, rdata, idata);
        idata[MessageIntData::type_id][p2] = MessageTypes::CONJUGATE;
    }
};

//! Functor for populating message data
struct CreateMessage
{
    /**
       @brief Create a message
       @param p [in, out] message
       @param n [in] message index
       @param rdata [in, out] message real data
       @param idata [in, out] message int data
       @param timestamp [in] time stamp
       @param pos [in] position
       @param iv [in] index
       @param sender_lp [in] sender LP id
       @param sender_entity [in] entity id
       @param receiver_lp [in] receiver LP id
       @param receiver_entity [in] entity id
     **/
    template <class P, class RT, class IT>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void operator()(
        P& p,
        const amrex::Long n,
        RT& rdata,
        IT& idata,
        const amrex::Real timestamp,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& pos,
        const amrex::IntVect& iv,
        const int sender_lp,
        const int sender_entity,
        const int receiver_lp,
        const int receiver_entity) const
    {
        AMREX_ASSERT(
            idata[MessageIntData::type_id][n] == MessageTypes::UNDEFINED);

        idata[MessageIntData::type_id][n] = MessageTypes::MESSAGE;
        rdata[MessageRealData::timestamp][n] = timestamp;

        AMREX_D_TERM(p.pos(0) = pos[0];, p.pos(1) = pos[1];, p.pos(2) = pos[2];)
        AMREX_D_TERM(idata[MessageIntData::i][n] = iv[0];
                     , idata[MessageIntData::j][n] = iv[1];
                     , idata[MessageIntData::k][n] = iv[2];)
        idata[MessageIntData::sender_lp][n] = sender_lp;
        idata[MessageIntData::sender_entity][n] = sender_entity;
        idata[MessageIntData::receiver_lp][n] = receiver_lp;
        idata[MessageIntData::receiver_entity][n] = receiver_entity;
    }
};

//! Functor for making a message undefined
struct MarkMessageUndefined
{
    /**
       @brief Mark a message as undefined
       @param n [in] message index
       @param rdata [in, out] message real data
       @param idata [in, out] message int data
     **/
    template <class RT, class IT>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    operator()(const amrex::Long n, RT& rdata, IT& idata) const
    {
        idata[MessageIntData::type_id][n] = MessageTypes::UNDEFINED;
        idata[MessageIntData::pair][n] = -1;
        rdata[MessageRealData::timestamp][n] = 0.0;
        rdata[MessageRealData::old_timestamp][n] = 0.0;
        rdata[MessageRealData::creation_time][n] = 0.0;
    }
};

} // namespace spades::particles
#endif
