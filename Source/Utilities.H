#ifndef UTILITIES_H
#define UTILITIES_H

#include <AMReX.H>
#include <AMReX_RealVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Random.H>
#include "Constants.H"

namespace spades {
/**
   @brief Cantor pairing function
   @param k1 [in] first number
   @param k2 [in] second number
   @return value of the Cantor pairing function
 **/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Long
pairing_function(const amrex::Long k1, const amrex::Long k2)
{
    return (k1 + k2) * (k1 + k2 + 1) / 2 + k2;
}

/**
   @brief Exponential distribution
   @param lambda [in] rate parameter
   @param engine [in] random number generator
   @return sample from an exponential distribution
 **/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
random_exponential(const amrex::Real lambda, amrex::RandomEngine const& engine)
{
    return std::log(1 - amrex::Random(engine)) / (-lambda);
}

/**
   @brief Exponential distribution
   @param lambda [in] rate parameter
   @return sample from an exponential distribution
 **/
amrex::Real random_exponential(const amrex::Real lambda);

/**
   @brief Maximum unsigned integer representation given a number of bits
   @param nbits [in] number of bits
   @return maximum integer representation
 **/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE std::uint64_t
max_representation(const int nbits)
{
    AMREX_ASSERT(nbits < 64);
    std::uint64_t result = 2;
    for (int i = 0; i < nbits; i++) {
        result *= 2;
    }
    return result;
}

/**
   @brief bitmask given a number of bits
   @param nbits [in] number of bits
   @return bitmask
 **/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE std::uint64_t bitmask(const int nbits)
{
    AMREX_ASSERT(nbits < 64);
    // Shift 1 by nbits (create power of 2)
    // Subtract 1 to fill all lower nbits with 1s
    // bitmask(16): 0xFFFF
    return (1ULL << nbits) - 1;
}

/**
   @brief Skip to the next line in stream
   @param is [in, out] stream
 **/
void goto_next_line(std::istream& is);

/**
   @brief Read and broadcast file to all ranks

   Equivalent to amrex::ParallelDescriptor::ReadAndBcastFile but for all ranks.

   @param filename [in] file to read
   @param charBuf [out] character buffer
   @param bExitOnError [in] exit if error
 **/
void read_file(
    const std::string& filename,
    amrex::Vector<char>& charBuf,
    bool bExitOnError);

}
} // namespace spades
#endif
