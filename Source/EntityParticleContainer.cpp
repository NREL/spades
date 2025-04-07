#include <utility>
#include "EntityParticleContainer.H"

namespace spades::particles {

EntityParticleContainer::EntityParticleContainer(
    amrex::AmrParGDB* par_gdb, int ngrow)
    : SpadesParticleContainer<
          EntityTypes,
          0,
          0,
          EntityRealData::ncomps,
          EntityIntData::ncomps>(par_gdb, ngrow)
{}

EntityParticleContainer::EntityParticleContainer(
    const amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<amrex::DistributionMapping>& dmap,
    const amrex::Vector<amrex::BoxArray>& ba,
    int ngrow)
    : SpadesParticleContainer<
          EntityTypes,
          0,
          0,
          EntityRealData::ncomps,
          EntityIntData::ncomps>(geom, dmap, ba, ngrow)
{}

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

    AMREX_D_TERM(m_int_data_names[EntityIntData::i] = "i";
                 m_writeflags_int[EntityIntData::i] = 1;
                 , m_int_data_names[EntityIntData::j] = "j";
                 m_writeflags_int[EntityIntData::j] = 1;
                 , m_int_data_names[EntityIntData::k] = "k";
                 m_writeflags_int[EntityIntData::k] = 1;)
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
        auto& aos = pti.GetArrayOfStructs();
        auto& soa = pti.GetStructOfArrays();
        const auto pdata = EntityData<ParticleType, RealVector, IntVector>(
            aos.dataPtr(), soa.GetRealData(), soa.GetIntData());

        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(
                     int i, int j, int AMREX_D_PICK(, , k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const int start = offset_arr(iv);
                for (int n = start; n < start + num_particles_arr(iv); n++) {
                    auto& p = pdata.m_aos[n];
                    p.id() = pid + n;
                    p.cpu() = my_proc;

                    pdata.mark_undefined(n);
                    pdata.m_idata[EntityIntData::owner][n] =
                        static_cast<int>(dom.index(iv));

                    AMREX_D_TERM(
                        p.pos(0) = plo[0] + (iv[0] + constants::HALF) * dx[0];
                        , p.pos(1) = plo[1] + (iv[1] + constants::HALF) * dx[1];
                        ,
                        p.pos(2) = plo[2] + (iv[2] + constants::HALF) * dx[2];)

                    AMREX_D_TERM(pdata.m_idata[EntityIntData::i][n] = iv[0];
                                 , pdata.m_idata[EntityIntData::j][n] = iv[1];
                                 , pdata.m_idata[EntityIntData::k][n] = iv[2];)
                }

                for (int n = start; n < start + entities_per_lp; n++) {
                    const amrex::Real ts = 0.0;

                    pdata.m_rdata[EntityRealData::timestamp][n] = ts;
                    pdata.m_idata[EntityIntData::type_id][n] =
                        EntityTypes::ENTITY;
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
        auto& aos = pti.GetArrayOfStructs();
        auto& soa = pti.GetStructOfArrays();
        const auto pdata = EntityData<ParticleType, RealVector, IntVector>(
            aos.dataPtr(), soa.GetRealData(), soa.GetIntData());

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pidx) noexcept {
            auto& p = pdata.m_aos[pidx];
            bool valid_type = false;
            for (int typ = 0; typ < EntityTypes::NTYPES; typ++) {
                valid_type = pdata.m_idata[EntityIntData::type_id][pidx] == typ;
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

    sort_impl(CompareParticle());
}

} // namespace spades::particles
