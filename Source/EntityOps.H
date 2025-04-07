#ifndef ENTITYOPS_H
#define ENTITYOPS_H

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Print.H>
#include "ParticleOps.H"
#include "EntityData.H"

namespace spades::particles {

//! Entity operations
template <class P, class RT, class IT>
struct EntityOps
    : public ParticleOps<
          EntityRealData::ncomps,
          EntityIntData::ncomps,
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
    EntityOps(
        P* aos,
        std::array<RT, EntityRealData::ncomps>& rdata,
        std::array<IT, EntityIntData::ncomps>& idata)
        : ParticleOps<EntityRealData::ncomps, EntityIntData::ncomps, P, RT, IT>(
              aos, rdata, idata) {};

    /**
       @brief Mark an entity as undefined
       @param n [in] message index
     **/
    void mark_undefined(const amrex::Long n) const
    {
        this->m_idata[EntityIntData::type_id][n] = EntityTypes::UNDEFINED;
        this->m_rdata[EntityRealData::timestamp][n] = 0.0;
    }
};

} // namespace spades::particles
#endif
