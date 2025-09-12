#ifndef MODEL_H
#define MODEL_H
#include "Factory.H"

namespace spades {

class Phold
{
public:
    static std::string identifier() { return "phold"; }

    explicit Phold();

    struct ProcessOp
    {
        AMREX_GPU_DEVICE void operator()() const
        {
            amrex::Print() << "Success!!!!" << std::endl;
        }
    };

    using ProcessType = ProcessOp;

    ProcessType process_op() const { return m_process_op; }

    ProcessOp m_process_op;
};

} // namespace spades
#endif /* MODEL_H */
