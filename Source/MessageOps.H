#ifndef MESSAGEOPS_H
#define MESSAGEOPS_H

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Print.H>
#include "MessageData.H"
#include "ParticleOps.H"

namespace spades::particles {

//! Message operations
template <class P, class RT, class IT>
struct MessageOps
    : public ParticleOps<
          MessageRealData::ncomps,
          MessageIntData::ncomps,
          P,
          RT,
          IT>
{

    /**
       @brief Constructor
       @param aos [in] particle AOS data
       @param rdata [in] particle real data
       @param idata [in] particle int data
     **/
    MessageOps(
        P* aos,
        std::array<RT, MessageRealData::ncomps>& rdata,
        std::array<IT, MessageIntData::ncomps>& idata)
        : ParticleOps<
              MessageRealData::ncomps,
              MessageIntData::ncomps,
              P,
              RT,
              IT>(aos, rdata, idata) {};

    /**
       @brief Mark a message as undefined
       @param n [in] message index
     **/
    void mark_undefined(const amrex::Long n) const
    {
        this->m_idata[MessageIntData::type_id][n] = MessageTypes::UNDEFINED;
        this->m_idata[MessageIntData::pair][n] = -1;
        this->m_rdata[MessageRealData::timestamp][n] = 0.0;
        this->m_rdata[MessageRealData::old_timestamp][n] = 0.0;
        this->m_rdata[MessageRealData::creation_time][n] = 0.0;
    }

    /**
       @brief Create a message
       @param n [in] message index
       @param timestamp [in] time stamp
       @param pos [in] position
       @param iv [in] index
       @param sender_lp [in] sender LP id
       @param sender_entity [in] entity id
       @param receiver_lp [in] receiver LP id
       @param receiver_entity [in] entity id
     **/
    void create_message(
        const amrex::Long n,
        const amrex::Real timestamp,
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& pos,
        const amrex::IntVect& iv,
        const int sender_lp,
        const int sender_entity,
        const int receiver_lp,
        const int receiver_entity) const
    {
        AMREX_ASSERT(
            this->m_idata[MessageIntData::type_id][n] ==
            MessageTypes::UNDEFINED);

        this->m_idata[MessageIntData::type_id][n] = MessageTypes::MESSAGE;
        this->m_rdata[MessageRealData::timestamp][n] = timestamp;

        auto& p = this->m_aos[n];
        AMREX_D_TERM(p.pos(0) = pos[0];, p.pos(1) = pos[1];, p.pos(2) = pos[2];)
        AMREX_D_TERM(this->m_idata[MessageIntData::i][n] = iv[0];
                     , this->m_idata[MessageIntData::j][n] = iv[1];
                     , this->m_idata[MessageIntData::k][n] = iv[2];)
        this->m_idata[MessageIntData::sender_lp][n] = sender_lp;
        this->m_idata[MessageIntData::sender_entity][n] = sender_entity;
        this->m_idata[MessageIntData::receiver_lp][n] = receiver_lp;
        this->m_idata[MessageIntData::receiver_entity][n] = receiver_entity;
    }
};

} // namespace spades::particles
#endif
