#ifndef MESSAGEPARTICLECONTAINER_H
#define MESSAGEPARTICLECONTAINER_H
#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_Random.H>
#include "SpadesParticleContainer.H"
#include "MessageData.H"
#include "MessageOps.H"
#include "Utilities.H"

/**
   @brief SPADES particles
 **/
namespace spades::particles {

//! Main SPADES message container
class MessageParticleContainer
    : public SpadesParticleContainer<
          MessageTypes::NTYPES,
          0,
          0,
          MessageRealData::ncomps,
          MessageIntData::ncomps>
{
public:
    /**
       @brief Class identifier name
       @return class identifier
     **/
    static std::string identifier() { return "messages"; }

    /**
       @brief Constructor
       @param par_gdb [in] particle database
       @param ngrow [in] number of grow cells
     **/
    explicit MessageParticleContainer(amrex::AmrParGDB* par_gdb, int ngrow = 0);

    /**
       @brief Constructor
       @param geom [in] geometry
       @param dmap [in] distribution map
       @param ba [in] box array
       @param ngrow [in] number of grow cells
     **/
    explicit MessageParticleContainer(
        const amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<amrex::DistributionMapping>& dmap,
        const amrex::Vector<amrex::BoxArray>& ba,
        int ngrow = 0);

    /**
       @brief Initialize the messages
       @param model [in] model
     **/
    template <typename Model>
    void initialize_messages(const Model& model);

    //! Sort the messages
    void sort() override;

    //! Update the undefined messages
    void update_undefined();

    //! Resolve message pairs (remove message/anti-message pairs)
    void resolve_pairs();

    //! Compute the minimum time stamp of the messages
    amrex::Real compute_gvt();

    /**
       @brief Perform garbage collection
       @param gvt [in] global virtual time
     **/
    void garbage_collect(const amrex::Real gvt);

    void write_plot_file(const std::string& plt_filename) override
    {
        write_plot_file_impl(plt_filename, identifier());
    };

    //! Read user parameters
    void read_parameters() override;

    //! Initialize variable names
    void initialize_variable_names() override;
};

template <typename Model>
void MessageParticleContainer::initialize_messages(const Model& model)
{
    BL_PROFILE("spades::MessageParticleContainer::initialize_messages()");

    const auto& plo = Geom(LEV).ProbLoArray();
    const auto& dx = Geom(LEV).CellSizeArray();
    const auto& dom = Geom(LEV).Domain();
    const auto init_message_op = model.init_message_op();
    const auto messages_per_lp = init_message_op.m_messages_per_lp;
    const int total_messages_per_lp = 3 * messages_per_lp;
    AMREX_ALWAYS_ASSERT(total_messages_per_lp > messages_per_lp);

    for (amrex::MFIter mfi = MakeMFIter(LEV); mfi.isValid(); ++mfi) {
        DefineAndReturnParticleTile(LEV, mfi);
    }

    amrex::iMultiFab num_particles(
        ParticleBoxArray(LEV), ParticleDistributionMap(LEV), 1, 0,
        amrex::MFInfo());
    amrex::iMultiFab init_offsets(
        ParticleBoxArray(LEV), ParticleDistributionMap(LEV), 1, 0,
        amrex::MFInfo());
    num_particles.setVal(total_messages_per_lp);
    init_offsets.setVal(0);

    for (amrex::MFIter mfi = MakeMFIter(LEV); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();

        const auto ncells = static_cast<int>(box.numPts());
        const int* in = num_particles[mfi].dataPtr();
        int* out = init_offsets[mfi].dataPtr();
        const auto np = amrex::Scan::PrefixSum<int>(
            ncells, [=] AMREX_GPU_DEVICE(int i) -> int { return in[i]; },
            [=] AMREX_GPU_DEVICE(int i, int const& xi) { out[i] = xi; },
            amrex::Scan::Type::exclusive, amrex::Scan::retSum);

        const amrex::Long pid = ParticleType::NextID();
        ParticleType::NextID(pid + np);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            static_cast<amrex::Long>(pid + np) < amrex::LastParticleID,
            "Error: overflow on particle id numbers!");

        const auto my_proc = amrex::ParallelDescriptor::MyProc();
        const auto& offset_arr = init_offsets[mfi].const_array();
        const auto& num_particles_arr = num_particles[mfi].const_array();
        const auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
        auto& pti = GetParticles(LEV)[index];
        pti.resize(np);
        const auto parrs = particle_arrays(pti);
        amrex::ParallelForRNG(
            box, [=] AMREX_GPU_DEVICE(
                     int i, int j, int AMREX_D_PICK(, , k),
                     amrex::RandomEngine const& engine) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const int start = offset_arr(iv);
                for (int n = start; n < start + num_particles_arr(iv); n++) {
                    auto& p = parrs.m_aos[n];
                    p.id() = pid + n;
                    p.cpu() = my_proc;

                    MarkMessageUndefined()(n, parrs);
                    parrs.m_idata[MessageIntData::sender_lp][n] =
                        static_cast<int>(dom.index(iv));
                    parrs.m_idata[MessageIntData::sender_entity][n] = 0;
                    parrs.m_idata[MessageIntData::receiver_lp][n] =
                        static_cast<int>(dom.index(iv));
                    parrs.m_idata[MessageIntData::receiver_entity][n] = 0;

                    AMREX_D_TERM(
                        p.pos(0) = plo[0] + (iv[0] + constants::HALF) * dx[0];
                        , p.pos(1) = plo[1] + (iv[1] + constants::HALF) * dx[1];
                        ,
                        p.pos(2) = plo[2] + (iv[2] + constants::HALF) * dx[2];)

                    AMREX_D_TERM(parrs.m_idata[CommonIntData::i][n] = iv[0];
                                 , parrs.m_idata[CommonIntData::j][n] = iv[1];
                                 , parrs.m_idata[CommonIntData::k][n] = iv[2];)
                }

                for (int n = start; n < start + messages_per_lp; n++) {
                    init_message_op(parrs, n, engine);
                    parrs.m_idata[MessageIntData::pair_id][n] = -1;
                    parrs.m_idata[MessageIntData::pair_cpu][n] = -1;
                }
            });

        // This is necessary
        amrex::Gpu::streamSynchronize();
    }
    Redistribute();

    // Sanity check all initial particles
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MyParIter pti(*this, LEV); pti.isValid(); ++pti) {
        const size_t np = pti.numParticles();
        const auto parrs = particle_arrays(pti.GetParticleTile());

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pidx) noexcept {
            bool valid_type = false;
            for (int typ = 0; typ < MessageTypes::NTYPES; typ++) {
                valid_type = parrs.m_idata[CommonIntData::type_id][pidx] == typ;
                if (valid_type) {
                    break;
                }
            }
            AMREX_ASSERT(valid_type);
            AMREX_ASSERT(parrs.m_aos[pidx].id() >= 0);
        });
    }
}

} // namespace spades::particles
#endif
