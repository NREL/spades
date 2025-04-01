#ifndef MESSAGEOPS_H
#define MESSAGEOPS_H

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Print.H>
#include "ParticleOps.H"
#include "EntityData.H"

namespace spades::particles {

//! Functor for making a entity undefined
struct MarkEntityUndefined
{
    /**
       @brief Mark a entity as undefined
       @param n [in] entity index
       @param rdata [in, out] entity real data
       @param idata [in, out] entity int data
     **/
    template <class RT, class IT>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    operator()(const amrex::Long n, RT& rdata, IT& idata) const
    {
        idata[EntityIntData::type_id][n] = EntityTypes::UNDEFINED;
        rdata[EntityRealData::timestamp][n] = 0.0;
    }
};

} // namespace spades::particles
#endif
