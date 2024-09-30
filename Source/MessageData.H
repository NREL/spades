#ifndef MESSAGEDATA_H
#define MESSAGEDATA_H

#include <AMReX.H>

namespace spades::particles {
//! Message real data
struct RealData
{
    enum : std::uint8_t { timestamp, old_timestamp, creation_time, ncomps };
};

//! Message int data
struct IntData
{
    enum : std::uint8_t {
        AMREX_D_DECL(i, j, k),
        type_id,
        sender,
        receiver,
        pair,
        ncomps
    };
};

//! Message types
struct MessageTypes
{
    //! Message
    static constexpr int MESSAGE = 0;
    //! Processed message
    static constexpr int PROCESSED = 1;
    //! Anti-message
    static constexpr int ANTI = 2;
    //! Conjugate message (unsent anti-message)
    static constexpr int CONJUGATE = 3;
    //! Undefined message
    static constexpr int UNDEFINED = 4;
    //! Number of different message types
    static constexpr int NTYPES = 5;
};
} // namespace spades::particles
#endif
