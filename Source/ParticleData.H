#ifndef PARTICLEDATA_H
#define PARTICLEDATA_H

#include <AMReX.H>
#include <AMReX_SPACE.H>

namespace spades::particles {

//! Common particle real data
struct CommonRealData
{
    enum : std::uint8_t { timestamp, ncomps };
};

//! Common particle int data
struct CommonIntData
{
    enum : std::uint8_t { AMREX_D_DECL(i, j, k), type_id, ncomps };
};

} // namespace spades::particles
#endif
