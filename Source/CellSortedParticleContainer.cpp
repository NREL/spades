#include <utility>

#include "CellSortedParticleContainer.H"

namespace spades::particles {

ParticleContainerInfo::ParticleContainerInfo(std::string basename)
    : m_basename(std::move(basename))
{}

ParticleContainerInfo::~ParticleContainerInfo() = default;

CellSortedParticleContainer::CellSortedParticleContainer(
    amrex::AmrParGDB* par_gdb, int ngrow)
    : amrex::NeighborParticleContainer<RealData::ncomps, IntData::ncomps>(
          par_gdb, ngrow)
    , m_info("message")
    , m_ngrow(ngrow)
{

    const int lev = 0;

    m_real_data_names.resize(RealData::ncomps, "");
    m_writeflags_real.resize(RealData::ncomps, 0);
    m_int_data_names.resize(IntData::ncomps, "");
    m_writeflags_int.resize(IntData::ncomps, 0);

    m_real_data_names[RealData::timestamp] = "timestamp";
    m_writeflags_real[RealData::timestamp] = 1;

    m_int_data_names[IntData::type_id] = "type_id";
    m_writeflags_int[IntData::type_id] = 1;
    m_int_data_names[IntData::sender] = "sender";
    m_writeflags_int[IntData::sender] = 1;
    m_int_data_names[IntData::receiver] = "receiver";
    m_writeflags_int[IntData::receiver] = 1;
}

void CellSortedParticleContainer::initialize(const amrex::Real lookahead)
{
    BL_PROFILE("spades::CellSortedParticleContainer::initialize");

    const int lev = 0;
    const auto& plo = Geom(lev).ProbLoArray();
    const auto& dx = Geom(lev).CellSizeArray();
    const auto& dlo = Geom(lev).Domain().smallEnd();
    const auto& dhi = Geom(lev).Domain().bigEnd();

    // Some test particles
    // #ifdef _OPENMP
    // #pragma omp parallel
    // #endif
    //     for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
    //         const amrex::Box& box = mfi.tilebox();
    //         const int gid = mfi.index();
    //         const int tid = mfi.LocalTileIndex();
    //         auto& pti = GetParticles(lev)[std::make_pair(gid, tid)];

    //         for (amrex::IntVect iv = box.smallEnd(); iv <= box.bigEnd();
    //              box.next(iv)) {
    //             const auto test = amrex::Random_int(dhi[0] - dlo[0] + 1) +
    //             dlo[0]; amrex::IntVect iv_src(AMREX_D_DECL(0, 0, 0));
    //             amrex::IntVect iv_dest(AMREX_D_DECL(
    //                 amrex::Random_int(dhi[0] - dlo[0] + 1) + dlo[0],
    //                 amrex::Random_int(dhi[1] - dlo[1] + 1) + dlo[1],
    //                 amrex::Random_int(dhi[2] - dlo[2] + 1) + dlo[2]));
    //             amrex::IntVect iv_dest2(AMREX_D_DECL(
    //                 amrex::Random_int(dhi[0] - dlo[0] + 1) + dlo[0],
    //                 amrex::Random_int(dhi[1] - dlo[1] + 1) + dlo[1],
    //                 amrex::Random_int(dhi[2] - dlo[2] + 1) + dlo[2]));
    //             if (iv == iv_src) {
    //                 {
    //                     particles::CellSortedParticleContainer::ParticleType
    //                     p; p.id() = particles::CellSortedParticleContainer::
    //                         ParticleType::NextID();
    //                     p.cpu() = amrex::ParallelDescriptor::MyProc();

    //                     p.idata(particles::IntData::type_id) =
    //                         static_cast<int>(particles::MessageType::message);
    //                     p.idata(particles::IntData::sender) =
    //                     box.index(iv_src);
    //                     p.idata(particles::IntData::receiver) =
    //                     box.index(iv_dest);
    //                     p.rdata(particles::RealData::timestamp) =
    //                         random_exponential(1.0) + lookahead;

    //                     AMREX_D_TERM(
    //                         p.pos(0) = plo[0] + (iv_dest[0] + 0.5) * dx[0];
    //                         , p.pos(1) = plo[1] + (iv_dest[1] + 0.5) * dx[1];
    //                         , p.pos(2) = plo[2] + (iv_dest[2] + 0.5) *
    //                         dx[2];)

    //                     AMREX_D_TERM(p.idata(particles::IntData::i) =
    //                     iv_dest[0];
    //                                  , p.idata(particles::IntData::j) =
    //                                  iv_dest[1]; ,
    //                                  p.idata(particles::IntData::k) =
    //                                  iv_dest[2];)

    //                     pti.push_back(p);
    //                 }
    //                 {
    //                     particles::CellSortedParticleContainer::ParticleType
    //                     p; p.id() = particles::CellSortedParticleContainer::
    //                         ParticleType::NextID();
    //                     p.cpu() = amrex::ParallelDescriptor::MyProc();

    //                     p.idata(particles::IntData::type_id) =
    //                         static_cast<int>(particles::MessageType::message);
    //                     p.idata(particles::IntData::sender) =
    //                     box.index(iv_src);
    //                     p.idata(particles::IntData::receiver) =
    //                     box.index(iv_dest2);
    //                     p.rdata(particles::RealData::timestamp) =
    //                         random_exponential(1.0) + lookahead;

    //                     AMREX_D_TERM(
    //                         p.pos(0) = plo[0] + (iv_dest2[0] + 0.5) * dx[0];
    //                         , p.pos(1) = plo[1] + (iv_dest2[1] + 0.5) *
    //                         dx[1]; , p.pos(2) = plo[2] + (iv_dest2[2] + 0.5)
    //                         * dx[2];)

    //                     AMREX_D_TERM(p.idata(particles::IntData::i) =
    //                     iv_dest2[0];
    //                                  , p.idata(particles::IntData::j) =
    //                                  iv_dest2[1];
    //                                  ,
    //                                  p.idata(particles::IntData::k) =
    //                                  iv_dest2[2];)

    //                     pti.push_back(p);
    //                 }
    //                 { // creating another particle
    //                     particles::CellSortedParticleContainer::ParticleType
    //                     p; p.id() = particles::CellSortedParticleContainer::
    //                         ParticleType::NextID();
    //                     p.cpu() = amrex::ParallelDescriptor::MyProc();

    //                     p.idata(particles::IntData::type_id) =
    //                         static_cast<int>(particles::MessageType::message);
    //                     p.idata(particles::IntData::sender) =
    //                     box.index(iv_src);
    //                     p.idata(particles::IntData::receiver) =
    //                     box.index(iv_dest);
    //                     p.rdata(particles::RealData::timestamp) =
    //                         random_exponential(1.0) + lookahead;

    //                     AMREX_D_TERM(
    //                         p.pos(0) = plo[0] + (iv_dest[0] + 0.5) * dx[0];
    //                         , p.pos(1) = plo[1] + (iv_dest[1] + 0.5) * dx[1];
    //                         , p.pos(2) = plo[2] + (iv_dest[2] + 0.5) *
    //                         dx[2];)

    //                     AMREX_D_TERM(p.idata(particles::IntData::i) =
    //                     iv_dest[0];
    //                                  , p.idata(particles::IntData::j) =
    //                                  iv_dest[1]; ,
    //                                  p.idata(particles::IntData::k) =
    //                                  iv_dest[2];)

    //                     pti.push_back(p);
    //                 }
    //                 { // creating another particle
    //                     particles::CellSortedParticleContainer::ParticleType
    //                     p; p.id() = particles::CellSortedParticleContainer::
    //                         ParticleType::NextID();
    //                     p.cpu() = amrex::ParallelDescriptor::MyProc();

    //                     p.idata(particles::IntData::type_id) =
    //                         static_cast<int>(particles::MessageType::message);
    //                     p.idata(particles::IntData::sender) =
    //                     box.index(iv_src);
    //                     p.idata(particles::IntData::receiver) =
    //                     box.index(iv_dest2);
    //                     p.rdata(particles::RealData::timestamp) =
    //                         random_exponential(1.0) + lookahead;

    //                     AMREX_D_TERM(
    //                         p.pos(0) = plo[0] + (iv_dest2[0] + 0.5) * dx[0];
    //                         , p.pos(1) = plo[1] + (iv_dest2[1] + 0.5) *
    //                         dx[1]; , p.pos(2) = plo[2] + (iv_dest2[2] + 0.5)
    //                         * dx[2];)

    //                     AMREX_D_TERM(p.idata(particles::IntData::i) =
    //                     iv_dest2[0];
    //                                  , p.idata(particles::IntData::j) =
    //                                  iv_dest2[1];
    //                                  ,
    //                                  p.idata(particles::IntData::k) =
    //                                  iv_dest2[2];)

    //                     pti.push_back(p);
    //                 }
    //             }
    //         }
    //     }

    const int np_per_cell = 100;
    const int msg_per_cell = 10;
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.tilebox();
        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        auto& pti = GetParticles(lev)[std::make_pair(gid, tid)];

        for (amrex::IntVect iv = box.smallEnd(); iv <= box.bigEnd();
             box.next(iv)) {

            for (int i = 0; i < np_per_cell; i++) {
                particles::CellSortedParticleContainer::ParticleType p;
                p.id() = particles::CellSortedParticleContainer::ParticleType::
                    NextID();
                p.cpu() = amrex::ParallelDescriptor::MyProc();

                p.idata(particles::IntData::type_id) =
                    static_cast<int>(particles::MessageType::undefined);
                p.idata(particles::IntData::sender) = box.index(iv);
                p.idata(particles::IntData::receiver) = box.index(iv);

                AMREX_D_TERM(p.pos(0) = plo[0] + (iv[0] + 0.5) * dx[0];
                             , p.pos(1) = plo[1] + (iv[1] + 0.5) * dx[1];
                             , p.pos(2) = plo[2] + (iv[2] + 0.5) * dx[2];)

                AMREX_D_TERM(p.idata(particles::IntData::i) = iv[0];
                             , p.idata(particles::IntData::j) = iv[1];
                             , p.idata(particles::IntData::k) = iv[2];)

                if (i < msg_per_cell) {
                    p.rdata(particles::RealData::timestamp) =
                        random_exponential(1.0) + lookahead;
                    p.idata(particles::IntData::type_id) =
                        static_cast<int>(particles::MessageType::message);
                }

                pti.push_back(p);
            }
        }
    }
}

void CellSortedParticleContainer::assign_cell_lists()
{
    BL_PROFILE("spades::CellSortedParticleContainer::assign_cell_lists");

    const int lev = 0;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
        auto index = std::make_pair(pti.index(), pti.LocalTileIndex());

        auto& particle_tile = GetParticles(lev)[index];

        if (particle_tile.numParticles() == 0) {
            continue;
        }

        const amrex::Box& box = pti.tilebox();

        const amrex::Box& gbox = amrex::grow(box, m_ngrow);

        m_cell_lists[lev].at(index).build(particle_tile, gbox);
    }
}

void CellSortedParticleContainer::garbage_collect(const amrex::Real gvt)
{
    BL_PROFILE("spades::CellSortedParticleContainer::garbage_collect");

    const int lev = 0;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
        auto index = std::make_pair(pti.index(), pti.LocalTileIndex());

        auto& partice_tile = GetParticles(lev)[index];
        const size_t np = partice_tile.numParticles();
        auto& particles = partice_tile.GetArrayOfStructs();
        auto* pstruct = particles().dataPtr();

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE(long pindex) noexcept {
            auto& p = pstruct[pindex];
            if (p.rdata(RealData::timestamp) < gvt) {
                p.idata(IntData::type_id) =
                    static_cast<int>(MessageType::undefined);
            }
        });
    }
}

void CellSortedParticleContainer::update_cell_lists()
{
    BL_PROFILE("spades::CellSortedParticleContainer::update_cell_lists");

    const int lev = 0;

    bool needs_update = false;
    if (!m_cell_lists_initialized) {
        // this is the first call, so we must update
        m_cell_lists_initialized = true;
        needs_update = true;
    } else if (
        (m_BARef != this->ParticleBoxArray(lev).getRefID()) ||
        (m_DMRef != this->ParticleDistributionMap(lev).getRefID())) {
        // the grids have changed, so we must update
        m_BARef = this->ParticleBoxArray(lev).getRefID();
        m_DMRef = this->ParticleDistributionMap(lev).getRefID();
        needs_update = true;
    }

    if (!needs_update) {
        return;
    }

    if (static_cast<int>(m_cell_lists.size()) <= numLevels()) {
        m_cell_lists.resize(numLevels());
    }
    m_cell_lists[lev].clear();
    for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const int gid = mfi.index();
        const int tid = mfi.LocalTileIndex();
        const auto index = std::make_pair(gid, tid);
        m_cell_lists[lev].insert({index, CellList<ParticleType>()});
    }

    assign_cell_lists();
}

void CellSortedParticleContainer::write_plot_file(
    const std::string& plt_filename)
{
    WritePlotFile(
        plt_filename, "particles", m_writeflags_real, m_writeflags_int,
        m_real_data_names, m_int_data_names);
}
} // namespace spades::particles
