#ifndef RUNTIME_H
#define RUNTIME_H
#include "SPADES.H"
#include "Constants.H"

class RunTime
{
public:
    virtual ~RunTime() = default;
    virtual void run() = 0;
};

template <typename Model>
class TemplatedRunTime : public RunTime
{
public:
    explicit TemplatedRunTime(const Model& model) : m_spades(model) {}

    void run() override
    {
        amrex::Real start_time = amrex::ParallelDescriptor::second();
        amrex::Print() << "Initializing SPADES ..." << std::endl;
        m_spades.initialize_data();

        amrex::Real init_time =
            amrex::ParallelDescriptor::second() - start_time;
        amrex::ParallelDescriptor::ReduceRealMax(
            init_time, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::Print() << "Initialization successful. Time elapsed = "
                       << init_time << std::endl;

        m_spades.evolve();

        amrex::Real end_time = amrex::ParallelDescriptor::second() - start_time;
        amrex::ParallelDescriptor::ReduceRealMax(
            end_time, amrex::ParallelDescriptor::IOProcessorNumber());

        amrex::Print() << "\n==============================================="
                          "================================================="
                       << std::endl
                       << std::endl;
        amrex::Print() << "Time spent in initialize_data():    " << init_time
                       << std::endl;
        amrex::Print() << "Time spent in evolve():       "
                       << end_time - init_time << std::endl;
    }

private:
    spades::SPADES<Model> m_spades;
};

void create_and_run_model(const std::string& model_name)
{
    std::unique_ptr<RunTime> runtime;

    if (model_name == "phold") {
        spades::models::Phold model;
        runtime =
            std::make_unique<TemplatedRunTime<spades::models::Phold>>(model);
        runtime->run();
    } else {
        amrex::Abort("Invalid model type" + model_name);
    }
}

#endif
