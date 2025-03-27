#include <utility>
#include "EntityParticleContainer.H"

namespace spades::particles {

EntityParticleContainer::EntityParticleContainer(
    amrex::AmrParGDB* par_gdb, int ngrow)
    : SpadesParticleContainer<
          EntityTypes,
          EntityRealData::ncomps,
          EntityIntData::ncomps>(par_gdb, ngrow)
{
    read_parameters();
}

EntityParticleContainer::EntityParticleContainer(
    const amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<amrex::DistributionMapping>& dmap,
    const amrex::Vector<amrex::BoxArray>& ba,
    int ngrow)
    : SpadesParticleContainer<
          EntityTypes,
          EntityRealData::ncomps,
          EntityIntData::ncomps>(geom, dmap, ba, ngrow)
{
    read_parameters();
}

void EntityParticleContainer::read_parameters()
{
    SpadesParticleContainer::read_parameters();
    {
        amrex::ParmParse pp("spades");
        pp.query("entities_per_lp", m_entities_per_lp);
    }
}

void EntityParticleContainer::initialize_variable_names()
{
    BL_PROFILE("spades::EntityParticleContainer::initialize_variable_names()");

    m_real_data_names.resize(EntityRealData::ncomps, "");
    m_writeflags_real.resize(EntityRealData::ncomps, 0);
    m_int_data_names.resize(EntityIntData::ncomps, "");
    m_writeflags_int.resize(EntityIntData::ncomps, 0);

    m_real_data_names[EntityRealData::timestamp] = "timestamp";
    m_writeflags_real[EntityRealData::timestamp] = 1;

    m_int_data_names[EntityIntData::type_id] = "type_id";
    m_writeflags_int[EntityIntData::type_id] = 1;
    m_int_data_names[EntityIntData::owner] = "owner";
    m_writeflags_int[EntityIntData::owner] = 1;
}

void EntityParticleContainer::initialize_entities()
{
    BL_PROFILE("spades::EntityParticleContainer::initialize_entities()");

    const auto& plo = Geom(LEV).ProbLoArray();
    const auto& dx = Geom(LEV).CellSizeArray();
    const auto& dom = Geom(LEV).Domain();
    const auto entities_per_lp = m_entities_per_lp;

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
        auto& pti = DefineAndReturnParticleTile(LEV, mfi);
        pti.resize(np);
        auto* aos = pti.GetArrayOfStructs().dataPtr();
        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(
                     int i, int j, int AMREX_D_PICK(, , k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const int start = offset_arr(iv);
                for (int n = start; n < start + num_particles_arr(iv); n++) {
                    auto& p = aos[n];
                    p.id() = pid + n;
                    p.cpu() = my_proc;

                    MarkEntityUndefined()(p);
                    p.idata(EntityIntData::owner) =
                        static_cast<int>(dom.index(iv));

                    AMREX_D_TERM(
                        p.pos(0) = plo[0] + (iv[0] + constants::HALF) * dx[0];
                        , p.pos(1) = plo[1] + (iv[1] + constants::HALF) * dx[1];
                        ,
                        p.pos(2) = plo[2] + (iv[2] + constants::HALF) * dx[2];)

                    AMREX_D_TERM(p.idata(EntityIntData::i) = iv[0];
                                 , p.idata(EntityIntData::j) = iv[1];
                                 , p.idata(EntityIntData::k) = iv[2];)
                }

                for (int n = start; n < start + entities_per_lp; n++) {
                    auto& pent = aos[n];
                    const amrex::Real ts = 0.0;

                    pent.rdata(EntityRealData::timestamp) = ts;
                    pent.idata(EntityIntData::type_id) = EntityTypes::ENTITY;
                }
            });
    }
    Redistribute();

    // Sanity check all initial particles
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MyParIter pti(*this, LEV); pti.isValid(); ++pti) {
        auto index = std::make_pair(pti.index(), pti.LocalTileIndex());

        auto& particle_tile = GetParticles(LEV)[index];
        const size_t np = particle_tile.numParticles();
        auto& particles = particle_tile.GetArrayOfStructs();
        auto* pstruct = particles().dataPtr();

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pindex) noexcept {
            auto& p = pstruct[pindex];
            bool valid_type = false;
            for (int typ = 0; typ < EntityTypes::NTYPES; typ++) {
                valid_type = p.idata(EntityIntData::type_id) == typ;
                if (valid_type) {
                    break;
                }
            }
            AMREX_ASSERT(valid_type);
            AMREX_ASSERT(p.id() >= 0);
        });
    }
}

void EntityParticleContainer::sort()
{
    BL_PROFILE("spades::EntityParticleContainer::sort()");

    sort_impl(CompareEntity());
}

} // namespace spades::particles
