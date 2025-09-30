#include "Phold.H"

namespace spades::models {

Phold::Phold() { read_parameters(); }

void Phold::read_parameters()
{

    amrex::Array<amrex::Real, AMREX_SPACEDIM> plo{AMREX_D_DECL(0.0, 0.0, 0.0)};
    amrex::Array<amrex::Real, AMREX_SPACEDIM> phi{AMREX_D_DECL(0.0, 0.0, 0.0)};
    amrex::Array<int, AMREX_SPACEDIM> periodic{AMREX_D_DECL(0, 0, 0)};
    {
        amrex::ParmParse pp("geometry");
        pp.get("prob_lo", plo);
        pp.get("prob_hi", phi);
        pp.get("is_periodic", periodic);
    }

    amrex::Array<int, AMREX_SPACEDIM> ncells{AMREX_D_DECL(0, 0, 0)};
    {
        amrex::ParmParse pp("amr");
        pp.get("n_cell", ncells);
    }

    amrex::RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
        real_box.setLo(n, plo[n]);
        real_box.setHi(n, phi[n]);
    }

    amrex::IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    amrex::IntVect domain_hi(
        AMREX_D_DECL(ncells[0] - 1, ncells[1] - 1, ncells[2] - 1));

    const amrex::Box domain(domain_lo, domain_hi);
    int coord = 0;
    m_geom.define(domain, &real_box, coord, periodic.data());

    m_process_op.m_dx = m_geom.CellSizeArray();
    m_process_op.m_dom = m_geom.Domain();
    m_process_op.m_plo = m_geom.ProbLoArray();

    amrex::Real lookahead = 1.0;
    amrex::Real lambda = 1.0;
    int entities_per_lp = 1;
    int messages_per_lp = 1;
    {
        amrex::ParmParse pp("spades");
        pp.query("lookahead", lookahead);
        pp.query("lambda", lambda);
        pp.query("entities_per_lp", entities_per_lp);
        pp.query("messages_per_lp", messages_per_lp);
    }

    m_process_op.m_lookahead = lookahead;
    m_process_op.m_lambda = lambda;
    m_process_op.m_entities_per_lp = entities_per_lp;
    m_init_entity_op.m_entities_per_lp = entities_per_lp;
    m_init_message_op.m_messages_per_lp = messages_per_lp;
    m_init_message_op.m_lambda = lambda;
    m_init_message_op.m_lookahead = lookahead;
}

} // namespace spades::models
