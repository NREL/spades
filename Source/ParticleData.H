#ifndef PARTICLEDATA_H
#define PARTICLEDATA_H

#include <AMReX.H>

namespace spades::particles {

struct RealData
{
    enum { timestamp, old_timestamp, ncomps };
};

struct IntData
{
    enum { AMREX_D_DECL(i, j, k), type_id, sender, receiver, ncomps };
};

struct MessageTypes
{
    static constexpr int MESSAGE = 0;
    static constexpr int PROCESSED = 1;
    static constexpr int ANTI_MESSAGE = 2;
    static constexpr int CONJUGATE = 3;
    static constexpr int UNDEFINED = 4;
    static constexpr int NTYPES = 5;
};
} // namespace spades::particles
#endif
