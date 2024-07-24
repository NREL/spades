#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <AMReX_REAL.H>
/**
   @brief SPADES constants
 **/
namespace spades::constants {

static constexpr int N_STATES = 2;
static constexpr int LVT_IDX = 0;
static constexpr int RLB_IDX = 1;

static constexpr amrex::Real LOW_NUM = static_cast<amrex::Real>(
    std::numeric_limits<amrex::Real>::lowest() * 1e-10);

static constexpr amrex::Real SMALL_NUM = static_cast<amrex::Real>(
    std::numeric_limits<amrex::Real>::epsilon() * 1e10);

static constexpr amrex::Real EPS =
    static_cast<amrex::Real>(std::numeric_limits<amrex::Real>::epsilon());

static constexpr amrex::Real TOL = 1e-6;

static constexpr amrex::Real LARGE_NUM =
    static_cast<amrex::Real>(std::numeric_limits<amrex::Real>::max() * 1e-10);

static constexpr amrex::Real HALF = 0.5;

static constexpr int FILE_NAME_DIGITS = 5;
static constexpr int RNG_FILE_NAME_DIGITS = 10;

static constexpr int HEADER_FILE_PRECISION = 17;

static constexpr int MAX_ROLLBACK_ITERS = 10;

} // namespace spades::constants

#endif
