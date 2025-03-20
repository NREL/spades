#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <AMReX_REAL.H>
/**
   @brief SPADES constants
 **/
namespace spades::constants {
//! Number of states
static constexpr int N_STATES = 2;

//! Index of the local virtual time
static constexpr int LVT_IDX = 0;

//! Index of the rollback counts
static constexpr int RLB_IDX = 1;

//! A large negative number
static constexpr amrex::Real LOW_NUM = static_cast<amrex::Real>(
    std::numeric_limits<amrex::Real>::lowest() * 1e-10);

//! A number close to zero
static constexpr amrex::Real SMALL_NUM = static_cast<amrex::Real>(
    std::numeric_limits<amrex::Real>::epsilon() * 1e10);

//! A number very close to zero
static constexpr amrex::Real EPS =
    static_cast<amrex::Real>(std::numeric_limits<amrex::Real>::epsilon());

//! A tolerance
static constexpr amrex::Real TOL = 1e-6;

//! A large number
static constexpr amrex::Real LARGE_NUM =
    static_cast<amrex::Real>(std::numeric_limits<amrex::Real>::max() * 1e-10);

static constexpr amrex::Real HALF = 0.5;

static constexpr int TWO_TO_EIGHT = 256;

//! Number of digits to use in output file names
static constexpr int FILE_NAME_DIGITS = 5;

//! Number of digits to use in random seed file names
static constexpr int RNG_FILE_NAME_DIGITS = 10;

//! Precision to use in header files
static constexpr int HEADER_FILE_PRECISION = 17;

//! Maximum number of rollback iterations allowed
static constexpr int MAX_ROLLBACK_ITERS = 10;

#if (AMREX_SPACEDIM == 2)
//! Bits to encode the cell i value
static constexpr int I_NBITS = 16;

//! Bits to encode the cell j value
static constexpr int J_NBITS = 16;

//! Bits to encode the cell k value
static constexpr int K_NBITS = 0;

//! Bits to encode the type
static constexpr int TYPE_NBITS = 4;

//! Bits to encode the timestamp
static constexpr int TIMESTAMP_NBITS = 28;
#endif

//! Bits to encode the timestamp
static constexpr int TOTAL_NBITS = 64;

} // namespace spades::constants

#endif
