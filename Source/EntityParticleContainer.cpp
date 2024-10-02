#include <utility>
#include "EntityParticleContainer.H"

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#endif

namespace spades::particles {

EntityParticleContainer::EntityParticleContainer(
    amrex::AmrParGDB* par_gdb, int ngrow)
    : amrex::NeighborParticleContainer<
          EntityRealData::ncomps,
          EntityIntData::ncomps>(par_gdb, ngrow)
    , m_info(identifier())
    , m_ngrow(ngrow)
{
    const int nlevs_max = par_gdb->maxLevel() + 1;

    if (nlevs_max > 1) {
        amrex::Abort(
            "spades::SPADES::EntityParticleContainer::"
            "EntityParticleContainer(): not supporting multilevel right "
            "now");
    }

    initialize_vectors();
}

EntityParticleContainer::EntityParticleContainer(
    const amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<amrex::DistributionMapping>& dmap,
    const amrex::Vector<amrex::BoxArray>& ba,
    int ngrow)
    : amrex::NeighborParticleContainer<
          EntityRealData::ncomps,
          EntityIntData::ncomps>(geom, dmap, ba, {2}, ngrow)
    , m_info(identifier())
    , m_ngrow(ngrow)
{
    if (geom.size() > 1) {
        amrex::Abort(
            "spades::SPADES::EntityParticleContainer::"
            "EntityParticleContainer(): not supporting multilevel right "
            "now");
    }

    initialize_vectors();
}

void EntityParticleContainer::initialize_vectors()
{
    BL_PROFILE("spades::EntityParticleContainer::initialize_vectors()");

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

void EntityParticleContainer::initialize_state()
{
    BL_PROFILE("spades::EntityParticleContainer::initialize_state()");

    m_entity_counts.define(
        ParticleBoxArray(LEV), ParticleDistributionMap(LEV),
        EntityTypes::NTYPES, m_ngrow, amrex::MFInfo());

    m_offsets.define(
        ParticleBoxArray(LEV), ParticleDistributionMap(LEV),
        EntityTypes::NTYPES, m_ngrow, amrex::MFInfo());

    m_entity_counts.setVal(0);
    m_offsets.setVal(0);
}

void EntityParticleContainer::clear_state()
{
    BL_PROFILE("spades::EntityParticleContainer::clear_state()");

    m_entity_counts.clear();
    m_offsets.clear();
}

void EntityParticleContainer::update_counts()
{
    BL_PROFILE("spades::EntityParticleContainer::update_counts()");
    count_entities();
    count_offsets();
}

void EntityParticleContainer::count_entities()
{
    BL_PROFILE("spades::EntityParticleContainer::count_entities()");

    m_entity_counts.setVal(0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = MakeMFIter(LEV); mfi.isValid(); ++mfi) {

        const amrex::Box& box = mfi.tilebox();
        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        const auto& cnt_arr = m_entity_counts.array(mfi);
        auto& pti = GetParticles(LEV)[std::make_pair(gid, tid)];
        const auto& particles = pti.GetArrayOfStructs();
        const auto* pstruct = particles().dataPtr();
        const int np = pti.numParticles();

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pindex) noexcept {
            const auto& p = pstruct[pindex];
            const amrex::IntVect iv(AMREX_D_DECL(
                p.idata(EntityIntData::i), p.idata(EntityIntData::j),
                p.idata(EntityIntData::k)));

            if (box.contains(iv)) {
                amrex::Gpu::Atomic::AddNoRet(
                    &cnt_arr(iv, p.idata(EntityIntData::type_id)), 1);
            }
        });
    }
}

void EntityParticleContainer::count_offsets()
{
    BL_PROFILE("spades::EntityParticleContainer::count_offsets()");

    m_offsets.setVal(0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = MakeMFIter(LEV); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        const auto ncell = box.numPts();
        const auto& cnt_arr = m_entity_counts.const_array(mfi);
        const auto& offsets_arr = m_offsets.array(mfi);
        int* p_offsets = offsets_arr.dataPtr();
        amrex::Scan::PrefixSum<int>(
            ncell,
            [=] AMREX_GPU_DEVICE(int i) -> int {
                const auto iv = box.atOffset(i);
                int total_entities = 0;
                for (int typ = 0; typ < EntityTypes::NTYPES; typ++) {
                    total_entities += cnt_arr(iv, typ);
                }
                return total_entities;
            },
            [=] AMREX_GPU_DEVICE(int i, const int& xi) { p_offsets[i] = xi; },
            amrex::Scan::Type::exclusive, amrex::Scan::noRetSum);

        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(
                     int i, int j, int AMREX_D_PICK(, , k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                for (int typ = 1; typ < EntityTypes::NTYPES; typ++) {
                    offsets_arr(iv, typ) =
                        offsets_arr(iv, typ - 1) + cnt_arr(iv, typ - 1);
                }
            });
    }
}

void EntityParticleContainer::initialize_entities()
{
    BL_PROFILE("spades::EntityParticleContainer::initialize_entities()");

    const auto& plo = Geom(LEV).ProbLoArray();
    const auto& dx = Geom(LEV).CellSizeArray();
    const auto& dom = Geom(LEV).Domain();

    const int np_per_cell = 100;
    const int ent_per_cell = 1;
    AMREX_ALWAYS_ASSERT(np_per_cell > ent_per_cell);

    amrex::iMultiFab num_particles(
        ParticleBoxArray(LEV), ParticleDistributionMap(LEV), 1, 0,
        amrex::MFInfo());
    amrex::iMultiFab init_offsets(
        ParticleBoxArray(LEV), ParticleDistributionMap(LEV), 1, 0,
        amrex::MFInfo());
    num_particles.setVal(np_per_cell);
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
        amrex::ParallelForRNG(
            box, [=] AMREX_GPU_DEVICE(
                     int i, int j, int AMREX_D_PICK(, , k),
                     amrex::RandomEngine const& engine) noexcept {
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

                for (int n = start; n < start + ent_per_cell; n++) {
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
            AMREX_ALWAYS_ASSERT(valid_type);
            AMREX_ALWAYS_ASSERT(p.id() >= 0);
        });
    }
}

void EntityParticleContainer::sort_entities()
{
    // Taking inspiration from AMReX's SortParticlesByBin
    BL_PROFILE("spades::EntityParticleContainer::sort_entities()");

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = MakeMFIter(LEV); mfi.isValid(); ++mfi) {
        auto& particle_tile = ParticlesAt(LEV, mfi);
        const size_t np = particle_tile.numParticles();

        if (np == 0) {
            continue;
        }

        // BL_PROFILE_VAR(
        //     "spades::EntityParticleContainer::sort_entities::sort_prep",
        //     prep);
        amrex::Gpu::DeviceVector<amrex::Long> cell_list(np);
        auto* p_cell_list = cell_list.data();
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pindex) noexcept {
            p_cell_list[pindex] = pindex;
        });
        // amrex::Gpu::Device::synchronize();
        // BL_PROFILE_VAR_STOP(prep);

        // Sort particle indices based on the cell index
        // BL_PROFILE_VAR(
        //     "spades::EntityParticleContainer::sort_entities::sort",
        //     sort);
        const auto& particles = particle_tile.GetArrayOfStructs();
        const auto* pstruct = particles().dataPtr();
#ifdef AMREX_USE_GPU
#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
        thrust::sort(
            thrust::device, cell_list.begin(), cell_list.end(),
            [=] AMREX_GPU_DEVICE(
                const amrex::Long xi, const amrex::Long yi) noexcept {
                const auto& p1 = pstruct[xi];
                const auto& p2 = pstruct[yi];
                return CompareEntity()(p1, p2);
            });
#else
        // Perform sort on CPU, then copy back to device (not good)
        amrex::Vector<amrex::Long> h_cell_list(np, 0);
        amrex::Gpu::copy(
            amrex::Gpu::deviceToHost, cell_list.begin(), cell_list.end(),
            h_cell_list.begin());
        std::sort(
            h_cell_list.begin(), h_cell_list.end(),
            [=](const amrex::Long xi, const amrex::Long yi) {
                const auto& p1 = pstruct[xi];
                const auto& p2 = pstruct[yi];
                return CompareEntity()(p1, p2);
            });
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, h_cell_list.begin(), h_cell_list.end(),
            cell_list.begin());
#endif
#else
        std::sort(
            cell_list.begin(), cell_list.end(),
            [=](const amrex::Long xi, const amrex::Long yi) {
                const auto& p1 = pstruct[xi];
                const auto& p2 = pstruct[yi];
                return CompareEntity()(p1, p2);
            });
#endif
        // amrex::Gpu::Device::synchronize();
        // BL_PROFILE_VAR_STOP(sort);

        // Reorder the particles in memory
        // BL_PROFILE_VAR(
        //     "spades::EntityParticleContainer::sort_entities::"
        //     "ReorderParticles",
        //     reorder);
        ReorderParticles(LEV, mfi, cell_list.data());
        // amrex::Gpu::Device::synchronize();
        // BL_PROFILE_VAR_STOP(reorder);
    }
    update_counts();
}

void EntityParticleContainer::update_undefined()
{
    BL_PROFILE("spades::EntityParticleContainer::update_undefined()");

    const auto& plo = Geom(LEV).ProbLoArray();
    const auto& dx = Geom(LEV).CellSizeArray();
    const auto& dom = Geom(LEV).Domain();
    const int lower_count = m_lower_undefined_count;
    const int upper_count = m_upper_undefined_count;
    const int reset_count = m_reset_undefined_count;
    AMREX_ALWAYS_ASSERT(reset_count < upper_count);
    AMREX_ALWAYS_ASSERT(lower_count < reset_count);
    int n_removals = 0;

    EntityParticleContainer pc_adds(
        m_gdb->Geom(), m_gdb->DistributionMap(), m_gdb->boxArray(), ngrow());

    for (amrex::MFIter mfi = MakeMFIter(LEV); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        const auto& cnt_arr = m_entity_counts.const_array(mfi);
        const auto& offsets_arr = m_offsets.const_array(mfi);
        auto& particle_tile = GetParticles(LEV)[std::make_pair(gid, tid)];
        auto& particles = particle_tile.GetArrayOfStructs();
        auto* pstruct = particles().dataPtr();

        // remove particles
        // BL_PROFILE_VAR(
        //     "spades::EntityParticleContainer::update_undefined::remove",
        //     remove);
        const auto ncells = static_cast<int>(box.numPts());
        amrex::Gpu::DeviceVector<int> removals(ncells, 0);
        auto* p_removals = removals.data();
        amrex::ParallelFor(ncells, [=] AMREX_GPU_DEVICE(long icell) noexcept {
            const auto iv = box.atOffset(icell);
            const int current_count = cnt_arr(iv, EntityTypes::UNDEFINED);
            if (current_count > upper_count) {
                p_removals[icell] = current_count - reset_count;
            }
        });

        n_removals += amrex::Reduce::Sum(removals.size(), p_removals);

        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(
                     int i, int j, int AMREX_D_PICK(, , k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const auto idx = box.index(iv);
                const auto np = p_removals[idx];
                const auto getter = Get(iv, cnt_arr, offsets_arr, pstruct);
                for (int m = 0; m < np; m++) {
                    auto& p = getter(m, EntityTypes::UNDEFINED);
                    p.id() = -1;
                }
            });
        // amrex::Gpu::Device::synchronize();
        // BL_PROFILE_VAR_STOP(remove);

        // BL_PROFILE_VAR(
        //     "spades::EntityParticleContainer::update_undefined::compute_add",
        //     compute_add);
        amrex::Gpu::DeviceVector<int> additions(ncells, 0);
        auto* p_additions = additions.data();
        amrex::ParallelFor(ncells, [=] AMREX_GPU_DEVICE(long icell) noexcept {
            const auto iv = box.atOffset(icell);
            const int current_count = cnt_arr(iv, EntityTypes::UNDEFINED);
            if (lower_count > current_count) {
                p_additions[icell] = reset_count - current_count;
            }
        });

        amrex::Gpu::DeviceVector<int> update_offsets(ncells, 0);
        auto* p_update_offsets = update_offsets.data();
        const auto np = amrex::Scan::PrefixSum<int>(
            ncells,
            [=] AMREX_GPU_DEVICE(int i) -> int { return p_additions[i]; },
            [=] AMREX_GPU_DEVICE(int i, int const& xi) {
                p_update_offsets[i] = xi;
            },
            amrex::Scan::Type::exclusive, amrex::Scan::retSum);

        const amrex::Long pid = ParticleType::NextID();
        ParticleType::NextID(pid + np);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            static_cast<amrex::Long>(pid + np) < amrex::LastParticleID,
            "Error: overflow on particle id numbers!");

        auto& ptile_adds = pc_adds.DefineAndReturnParticleTile(LEV, mfi);
        ptile_adds.resize(np);
        const auto my_proc = amrex::ParallelDescriptor::MyProc();
        auto* aos = ptile_adds.GetArrayOfStructs().dataPtr();
        amrex::ParallelFor(ncells, [=] AMREX_GPU_DEVICE(long icell) noexcept {
            const int start = p_update_offsets[icell];
            const auto iv = box.atOffset(icell);
            for (int n = start; n < start + p_additions[icell]; n++) {
                auto& p = aos[n];
                p.id() = pid + n;
                p.cpu() = my_proc;

                MarkEntityUndefined()(p);
                p.idata(EntityIntData::owner) = static_cast<int>(dom.index(iv));

                AMREX_D_TERM(
                    p.pos(0) = plo[0] + (iv[0] + constants::HALF) * dx[0];
                    , p.pos(1) = plo[1] + (iv[1] + constants::HALF) * dx[1];
                    , p.pos(2) = plo[2] + (iv[2] + constants::HALF) * dx[2];)

                AMREX_D_TERM(p.idata(EntityIntData::i) = iv[0];
                             , p.idata(EntityIntData::j) = iv[1];
                             , p.idata(EntityIntData::k) = iv[2];)
            }
        });
        // amrex::Gpu::Device::synchronize();
        // BL_PROFILE_VAR_STOP(compute_add);
    }

    addParticles(pc_adds, true);

    amrex::ParallelAllReduce::Sum<int>(
        n_removals, amrex::ParallelContext::CommunicatorSub());
    if (n_removals > 0) {
        Redistribute(); // kind of a bummer
    }
}

void EntityParticleContainer::garbage_collect(const amrex::Real gvt)
{
    BL_PROFILE("spades::EntityParticleContainer::garbage_collect()");
    amrex::Abort("Not implemented");
}

void EntityParticleContainer::reposition_entities()
{
    BL_PROFILE("spades::EntityParticleContainer::reposition_entities()");

    const auto& plo = Geom(LEV).ProbLoArray();
    const auto& dx = Geom(LEV).CellSizeArray();
    const auto& dxi = Geom(LEV).InvCellSizeArray();
    const auto& dom = Geom(LEV).Domain();
    const int nbins = 500;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi = MakeMFIter(LEV); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        const auto& cnt_arr = m_entity_counts.const_array(mfi);
        const auto& offsets_arr = m_offsets.const_array(mfi);
        const auto index = std::make_pair(gid, tid);
        auto& pti = GetParticles(LEV)[index];
        auto& particles = pti.GetArrayOfStructs();
        auto* pstruct = particles().dataPtr();

        amrex::ParallelFor(
            box, [=] AMREX_GPU_DEVICE(
                     int i, int j, int AMREX_D_PICK(, , k)) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                const auto getter = Get(iv, cnt_arr, offsets_arr, pstruct);

                for (int typ = 0; typ < EntityTypes::NTYPES; typ++) {
                    AMREX_ALWAYS_ASSERT(cnt_arr(iv, typ) < nbins);
                    for (int n = 0; n < cnt_arr(iv, typ); n++) {
                        auto& p = getter(n, typ);

                        const amrex::IntVect piv(AMREX_D_DECL(
                            p.idata(EntityIntData::i),
                            p.idata(EntityIntData::j),
                            p.idata(EntityIntData::k)));
                        AMREX_ALWAYS_ASSERT(piv == iv);

                        AMREX_D_TERM(
                            p.pos(0) =
                                plo[0] + iv[0] * dx[0] +
                                (typ + 1) * dx[0] / (EntityTypes::NTYPES + 1);
                            , p.pos(1) = plo[1] + iv[1] * dx[1] +
                                         (n + 1) * dx[1] / nbins;
                            , p.pos(2) =
                                  plo[2] + (iv[2] + constants::HALF) * dx[2];)

                        // ensure the particle didn't change cells
                        AMREX_ALWAYS_ASSERT(
                            piv == getParticleCell(p, plo, dxi, dom));
                    }
                }
            });
    }
}

void EntityParticleContainer::write_plot_file(const std::string& plt_filename)
{
    BL_PROFILE("spades::EntityParticleContainer::write_plot_file()");
    reposition_entities();
    WritePlotFile(
        plt_filename, identifier(), m_writeflags_real, m_writeflags_int,
        m_real_data_names, m_int_data_names);
}
} // namespace spades::particles
