#ifndef ENTITYDATA_H
#define ENTITYDATA_H

#include <AMReX.H>
#include "ParticleData.H"

namespace spades::particles {

//! Entity real data
struct EntityRealData
{
    enum : std::uint8_t { ncomps = CommonRealData::ncomps };
};

//! Entity int data
struct EntityIntData
{
    enum : std::uint8_t { owner = CommonIntData::ncomps, ncomps };
};

//! Entity types
struct EntityTypes
{
    //! Entity
    static constexpr int ENTITY = 0;
    //! Processed entity
    static constexpr int BACKUP = 1;
    //! Undefined entity
    static constexpr int UNDEFINED = 2;
    //! Number of different entity types
    static constexpr int NTYPES = 3;
};
} // namespace spades::particles
#endif
