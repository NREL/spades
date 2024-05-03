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

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Long
pairing_function(const amrex::Long k1, const amrex::Long k2)
{
    // Cantor pairing function
    return (k1 + k2) * (k1 + k2 + 1) / 2 + k2;
}

// skip to next line in Header
void goto_next_line(std::istream& is);

// equivalent to amrex::ParallelDescriptor::ReadAndBcastFile but for all ranks
void read_file(
    const std::string& filename,
    amrex::Vector<char>& charBuf,
    bool bExitOnError);

} // namespace spades
#endif
