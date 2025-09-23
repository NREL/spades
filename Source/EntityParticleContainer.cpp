#include <utility>
#include "EntityParticleContainer.H"

namespace spades::particles {

EntityParticleContainer::EntityParticleContainer(
    amrex::AmrParGDB* par_gdb, int ngrow)
    : SpadesParticleContainer<
          EntityTypes::NTYPES,
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
          EntityTypes::NTYPES,
          0,
          0,
          EntityRealData::ncomps,
          EntityIntData::ncomps>(geom, dmap, ba, ngrow)
{}

void EntityParticleContainer::read_parameters()
{
    SpadesParticleContainer::read_parameters();
}

void EntityParticleContainer::initialize_variable_names()
{
    BL_PROFILE("spades::EntityParticleContainer::initialize_variable_names()");

    m_real_data_names.resize(EntityRealData::ncomps, "");
    m_writeflags_real.resize(EntityRealData::ncomps, 0);
    m_int_data_names.resize(EntityIntData::ncomps, "");
    m_writeflags_int.resize(EntityIntData::ncomps, 0);

    m_real_data_names[CommonRealData::timestamp] = "timestamp";
    m_writeflags_real[CommonRealData::timestamp] = 1;

    AMREX_D_TERM(
        m_int_data_names[CommonIntData::i] = "i";
        m_writeflags_int[CommonIntData::i] = 1;
        , m_int_data_names[CommonIntData::j] = "j";
        m_writeflags_int[CommonIntData::j] = 1;
        , m_int_data_names[CommonIntData::k] = "k";
        m_writeflags_int[CommonIntData::k] = 1;)
    m_int_data_names[CommonIntData::type_id] = "type_id";
    m_writeflags_int[CommonIntData::type_id] = 1;
    m_int_data_names[EntityIntData::owner] = "owner";
    m_writeflags_int[EntityIntData::owner] = 1;
}

void EntityParticleContainer::sort()
{
    BL_PROFILE("spades::EntityParticleContainer::sort()");

    sort_impl(CompareParticle());
}
} // namespace spades::particles
