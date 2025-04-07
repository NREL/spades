#ifndef ENTITYOPS_H
#define ENTITYOPS_H

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Print.H>
#include "ParticleOps.H"
#include "EntityData.H"

namespace spades::particles {

// Entity data
template <class P, class RT, class IT>
struct EntityOps
    : public ParticleOps<
          EntityRealData::ncomps,
          EntityIntData::ncomps,
          P,
          RT,
          IT>
{

    EntityOps(
        P* aos,
        std::array<RT, EntityRealData::ncomps>& rdata,
        std::array<IT, EntityIntData::ncomps>& idata)
        : ParticleOps<EntityRealData::ncomps, EntityIntData::ncomps, P, RT, IT>(
              aos, rdata, idata) {};

    void mark_undefined(const amrex::Long n) const
    {
        this->m_idata[EntityIntData::type_id][n] = EntityTypes::UNDEFINED;
        this->m_rdata[EntityRealData::timestamp][n] = 0.0;
    }
};

} // namespace spades::particles
#endif
