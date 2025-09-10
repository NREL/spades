#ifndef MODEL_H
#define MODEL_H
#include "Factory.H"

namespace spades {
class Model : public spades::Factory<Model>
{
public:
    static std::string base_identifier() { return "Model"; }
    ~Model() override = default;

    //! String representation of this model for I/O
    virtual std::string model_name() const = 0;
};

class Phold : public Model::Register<Phold>
{
public:
    static std::string identifier() { return "phold"; }

    explicit Phold();

    std::string model_name() const override { return identifier(); };

    struct ProcessOp
    {};

    using ProcessType = ProcessOp;

    ProcessType process_instance() const { return m_process_op; }

    ProcessOp m_process_op;
};

} // namespace spades
#endif /* MODEL_H */
