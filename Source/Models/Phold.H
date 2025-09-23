#ifndef PHOLD_H
#define PHOLD_H
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include "Source/Constants.H"
#include "Source/Utilities.H"
#include "Source/MessageOps.H"
#include "Source/MessageData.H"
#include "Source/EntityData.H"

namespace spades::models {

//! PHOLD model
class Phold
{
public:
    /**
       @brief Class identifier name
       @return class identifier
     **/
    static std::string identifier() { return "phold"; }

    //! Constructor
    explicit Phold();

    //! Read user parameters
    void read_parameters();

    [[nodiscard]] const amrex::Geometry& geom() const { return m_geom; }

    /**
     * @struct ProcessOp
     * @brief A struct that defines an operator for processing messages
     */
    struct ProcessOp
    {
        /**
         * @brief Callable operator to process a message
         *
         * @tparam MParrs type of the message particle array
         * @tparam EParrs type of the entity particle array
         * @param msg_parrs message particle array
         * @param ent_parrs entity particle array
         * @param iv cell index of the particle
         * @param rcv_msg index of the received message
         * @param snd_msg index of the sending message particle
         * @param ent global ID of the entity
         * @param ent_idx index of the entity
         * @param engine random engine
         */
        template <class MParrs, class EParrs>
        AMREX_GPU_DEVICE void operator()(
            const MParrs& msg_parrs,
            const EParrs& ent_parrs,
            const amrex::IntVect& iv,
            const int rcv_msg,
            const int snd_msg,
            const int ent,
            const int ent_idx,
            amrex::RandomEngine const& engine) const
        {
            const auto& dlo = m_dom.smallEnd();
            const auto& dhi = m_dom.bigEnd();
            const auto ts =
                msg_parrs
                    .m_rdata[particles::CommonRealData::timestamp][rcv_msg];

            msg_parrs
                .m_rdata[particles::MessageRealData::old_timestamp][rcv_msg] =
                ent_parrs
                    .m_rdata[particles::CommonRealData::timestamp][ent_idx];
            ent_parrs.m_rdata[particles::CommonRealData::timestamp][ent_idx] =
                ts;
            msg_parrs.m_idata[particles::CommonIntData::type_id][rcv_msg] =
                particles::MessageTypes::PROCESSED;

            const auto ent_lvt =
                ent_parrs
                    .m_rdata[particles::CommonRealData::timestamp][ent_idx];
            const amrex::IntVect iv_dest(AMREX_D_DECL(
                amrex::Random_int(dhi[0] - dlo[0] + 1, engine) + dlo[0],
                amrex::Random_int(dhi[1] - dlo[1] + 1, engine) + dlo[1],
                amrex::Random_int(dhi[2] - dlo[2] + 1, engine) + dlo[2]));
            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> pos = {
                AMREX_D_DECL(
                    m_plo[0] + (iv_dest[0] + constants::HALF) * m_dx[0],
                    m_plo[1] + (iv_dest[1] + constants::HALF) * m_dx[1],
                    m_plo[2] + (iv_dest[2] + constants::HALF) * m_dx[2])};
            const int rcv_ent =
                static_cast<int>(amrex::Random_int(m_entities_per_lp, engine));
            const amrex::Real next_ts =
                ent_lvt + random_exponential(m_lambda, engine) + m_lookahead;

            particles::CreateMessage()(
                snd_msg, msg_parrs, next_ts, pos, iv_dest,
                static_cast<int>(m_dom.index(iv)), ent,
                static_cast<int>(m_dom.index(iv_dest)), rcv_ent);
        }

        //! Lookahead value (optional user input)
        amrex::Real m_lookahead{1.0};

        //! Width of exponential distribution (optional user input)
        //! smaller values of lambda = larger variance in random values
        //! larger values of lambda = smaller variance in random values
        amrex::Real m_lambda{1.0};

        //! Number of entities per logical process (optional user input)
        int m_entities_per_lp{1};

        //! Spatial discretization
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_dx;

        //! Domain
        amrex::Box m_dom;

        //! Problem low bounds
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_plo;
    };

    using ProcessType = ProcessOp;

    /**
       @brief Message processing operator
       @return ProcessType message processing operator
     **/
    [[nodiscard]] ProcessType process_op() const { return m_process_op; }

    // Message processing operator
    ProcessOp m_process_op;

    /**
     * @struct InitEntityOp
     * @brief A struct that defines an operator for initializing entities
     */
    struct InitEntityOp
    {
        /**
         * @brief Callable operator to initialize an entity
         *
         * @tparam EParrs type of the entity particle array
         * @param ent_parrs entity particle array
         * @param n index of the particle
         */
        template <class EParrs>
        AMREX_GPU_DEVICE void
        operator()(const EParrs& ent_parrs, const int n) const
        {
            const amrex::Real ts = 0.0;

            ent_parrs.m_rdata[particles::CommonRealData::timestamp][n] = ts;
            ent_parrs.m_idata[particles::CommonIntData::type_id][n] =
                particles::EntityTypes::ENTITY;
        }

        //! Number of entities per logical process (optional user input)
        int m_entities_per_lp{1};
    };

    using InitEntityType = InitEntityOp;

    /**
       @brief Initialize entity operator
       @return InitEntityType initialize entity operator
     **/
    [[nodiscard]] InitEntityType init_entity_op() const
    {
        return m_init_entity_op;
    }

    // Initialize entity operator
    InitEntityOp m_init_entity_op;

    /**
     * @struct InitMessageOp
     * @brief A struct that defines an operator for initializing entities
     */
    struct InitMessageOp
    {
        /**
         * @brief Callable operator to initialize an message
         *
         * @tparam MParrs type of the message particle array
         * @param ent_parrs message particle array
         * @param n index of the particle
         */
        template <class MParrs>
        AMREX_GPU_DEVICE void operator()(
            const MParrs& msg_parrs,
            const int n,
            amrex::RandomEngine const& engine) const
        {
            const amrex::Real ts =
                random_exponential(m_lambda, engine) + m_lookahead;
            msg_parrs.m_idata[particles::CommonIntData::type_id][n] =
                particles::MessageTypes::MESSAGE;
            msg_parrs.m_rdata[particles::CommonRealData::timestamp][n] = ts;
            msg_parrs.m_rdata[particles::MessageRealData::creation_time][n] =
                0.0;
        }

        //! Lookahead value (optional user input)
        amrex::Real m_lookahead{1.0};

        //! Initial messages per LP
        int m_messages_per_lp{1};

        //! Width of exponential distribution (optional user input)
        //! smaller values of lambda = larger variance in random values
        //! larger values of lambda = smaller variance in random values
        amrex::Real m_lambda{1.0};
    };

    using InitMessageType = InitMessageOp;

    /**
       @brief Initialize message operator
       @return InitMessageType initialize message operator
     **/
    [[nodiscard]] InitMessageType init_message_op() const
    {
        return m_init_message_op;
    }

    // Initialize message operator
    InitMessageOp m_init_message_op;

private:
    // Domain geometry
    amrex::Geometry m_geom;
};

} // namespace spades::models
#endif /* PHOLD_H */
