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

enum struct MessageType : int {
    anti_message = -1,
    undefined,
    message,
    processed
};
} // namespace spades::particles
#endif
