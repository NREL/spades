#ifndef ENTITYPARTICLECONTAINER_H
#define ENTITYPARTICLECONTAINER_H
#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_Random.H>
#include "SpadesParticleContainer.H"
#include "EntityData.H"
#include "EntityOps.H"
#include "Utilities.H"

/**
   @brief SPADES particles
 **/
namespace spades::particles {

//! Main SPADES entity container
class EntityParticleContainer
    : public SpadesParticleContainer<
          EntityTypes::NTYPES,
          0,
          0,
          EntityRealData::ncomps,
          EntityIntData::ncomps>
{
public:
    /**
       @brief Class identifier name
       @return class identifier
     **/
    static std::string identifier() { return "entities"; }

    /**
       @brief Constructor
       @param par_gdb [in] particle database
       @param ngrow [in] number of grow cells
     **/
    explicit EntityParticleContainer(amrex::AmrParGDB* par_gdb, int ngrow = 0);

    /**
       @brief Constructor
       @param geom [in] geometry
       @param dmap [in] distribution map
       @param ba [in] box array
       @param ngrow [in] number of grow cells
     **/
    explicit EntityParticleContainer(
        const amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<amrex::DistributionMapping>& dmap,
        const amrex::Vector<amrex::BoxArray>& ba,
        int ngrow = 0);

    //! Initialize the entities
    template <typename Model>
    void initialize_entities(const Model& model);

    //! Sort the entities
    void sort() override;

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
void EntityParticleContainer::initialize_entities(const Model& model)
{
    BL_PROFILE("spades::EntityParticleContainer::initialize_entities()");

    const auto& plo = Geom(LEV).ProbLoArray();
    const auto& dx = Geom(LEV).CellSizeArray();
    const auto& dom = Geom(LEV).Domain();
    const auto init_entity_op = model.init_entity_op();
    const auto entities_per_lp = init_entity_op.m_entities_per_lp;

    for (amrex::MFIter mfi = MakeMFIter(LEV); mfi.isValid(); ++mfi) {
        DefineAndReturnParticleTile(LEV, mfi);
    }

    amrex::iMultiFab num_particles(
        ParticleBoxArray(LEV), ParticleDistributionMap(LEV), 1, 0,
        amrex::MFInfo());
    amrex::iMultiFab init_offsets(
        ParticleBoxArray(LEV), ParticleDistributionMap(LEV), 1, 0,
        amrex::MFInfo());
    num_particles.setVal(entities_per_lp);
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

        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(
                     int i, int j, int AMREX_D_PICK(, , k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const int start = offset_arr(iv);
                for (int n = start; n < start + num_particles_arr(iv); n++) {
                    auto& p = parrs.m_aos[n];
                    p.id() = pid + n;
                    p.cpu() = my_proc;

                    MarkEntityUndefined()(n, parrs);
                    parrs.m_idata[EntityIntData::owner][n] =
                        static_cast<int>(dom.index(iv));

                    AMREX_D_TERM(
                        p.pos(0) = plo[0] + (iv[0] + constants::HALF) * dx[0];
                        , p.pos(1) = plo[1] + (iv[1] + constants::HALF) * dx[1];
                        ,
                        p.pos(2) = plo[2] + (iv[2] + constants::HALF) * dx[2];)

                    AMREX_D_TERM(parrs.m_idata[CommonIntData::i][n] = iv[0];
                                 , parrs.m_idata[CommonIntData::j][n] = iv[1];
                                 , parrs.m_idata[CommonIntData::k][n] = iv[2];)
                }

                for (int n = start; n < start + entities_per_lp; n++) {
                    init_entity_op(parrs, n);
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
            for (int typ = 0; typ < EntityTypes::NTYPES; typ++) {
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
