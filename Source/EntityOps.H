#ifndef ENTITYOPS_H
#define ENTITYOPS_H

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
       @param parrs [in, out] entity arrays
     **/
    template <class PArrs>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    operator()(const amrex::Long n, PArrs& parrs) const
    {
        parrs.m_idata[CommonIntData::type_id][n] = EntityTypes::UNDEFINED;
        parrs.m_rdata[CommonRealData::timestamp][n] = 0.0;
    }
};

} // namespace spades::particles
#endif
