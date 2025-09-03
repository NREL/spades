#ifndef MESSAGEDATA_H
#define MESSAGEDATA_H

#include <AMReX.H>
#include "ParticleData.H"

namespace spades::particles {
//! Message real data
struct MessageRealData
{
    enum : std::uint8_t {
        old_timestamp = CommonRealData::ncomps,
        creation_time,
        ncomps
    };
};

//! Message int data
struct MessageIntData
{
    enum : std::uint8_t {
        sender_lp = CommonIntData::ncomps,
        receiver_lp,
        sender_entity,
        receiver_entity,
        pair_cpu,
        pair_id,
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
