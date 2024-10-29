#ifndef MESSAGEOPS_H
#define MESSAGEOPS_H

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Print.H>
#include "ParticleOps.H"
#include "EntityData.H"

namespace spades::particles {

/**
   @brief Functor for comparing entities

   This is used to sort entities based on cell index, entity type,
   and timestamp (in that order).
 **/
struct CompareEntity
{
    /**
       @brief Compare entities
       @param p1 [in] first entity
       @param p2 [in] second entity
       @return boolean of p1 < p2 based on criteria
     **/
    template <class P>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool
    operator()(const P& p1, const P& p2) const
    {
        // sort by iv, then entity type then by timestamp
        const amrex::IntVect piv1(AMREX_D_DECL(
            p1.idata(EntityIntData::i), p1.idata(EntityIntData::j),
            p1.idata(EntityIntData::k)));
        const amrex::IntVect piv2(AMREX_D_DECL(
            p2.idata(EntityIntData::i), p2.idata(EntityIntData::j),
            p2.idata(EntityIntData::k)));

        const auto m1 = p1.idata(EntityIntData::type_id);
        const auto m2 = p2.idata(EntityIntData::type_id);
        const auto t1 = p1.rdata(EntityRealData::timestamp);
        const auto t2 = p2.rdata(EntityRealData::timestamp);
        return (piv1 < piv2) ||
               (piv1 == piv2 && ((m1 < m2) || (m1 == m2 && t1 < t2)));
    }
};

//! Functor for making a entity undefined
struct MarkEntityUndefined
{
    /**
       @brief Mark a entity as undefined
       @param p [in, out] entity
     **/
    template <class P>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void operator()(P& p) const
    {
        p.idata(EntityIntData::type_id) = EntityTypes::UNDEFINED;
        p.rdata(EntityRealData::timestamp) = 0.0;
    }
};

} // namespace spades::particles
#endif
