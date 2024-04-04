#ifndef PARTICLEDATA_H
#define PARTICLEDATA_H

#include <AMReX.H>

namespace spades::particles {

struct RealData
{
    enum { timestamp, ncomps };
};

struct IntData
{
    enum { AMREX_D_DECL(i, j, k), type_id, sender, receiver, ncomps };
};

struct MessageTypes
{
    static constexpr int anti_message = 0;
    static constexpr int undefined = 1;
    static constexpr int message = 2;
    static constexpr int processed = 3;
    static constexpr int ntypes = 4;
};
} // namespace spades::particles
#endif
