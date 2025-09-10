#include "Model.H"

namespace spades {

void Phold::initialize(const amrex::Geometry& geom)
{
    read_parameters();

    m_process_op.m_dx = geom.CellSizeArray();
    m_process_op.m_dom = geom.Domain();
    m_process_op.m_plo = geom.ProbLoArray();
}

void Phold::read_parameters()
{
    amrex::Real lookahead = 1.0;
    amrex::Real lambda = 1.0;
    int entities_per_lp = 1;
    {
        amrex::ParmParse pp("spades");
        pp.query("lookahead", lookahead);
        pp.query("lambda", lambda);
        pp.query("entities_per_lp", entities_per_lp);
    }

    m_process_op.m_lookahead = lookahead;
    m_process_op.m_lambda = lambda;
    m_process_op.m_entities_per_lp = entities_per_lp;
}

} // namespace spades
