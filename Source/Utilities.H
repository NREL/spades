#ifndef UTILITIES_H
#define UTILITIES_H

#include <AMReX.H>
#include <AMReX_RealVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "Constants.H"

namespace spades {
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
random_exponential(const amrex::Real lambda, amrex::RandomEngine const& engine)
{

    return std::log(1 - amrex::Random(engine)) / (-lambda);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
random_exponential(const amrex::Real lambda)
{

    return std::log(1 - amrex::Random()) / (-lambda);
}
} // namespace spades
#endif
