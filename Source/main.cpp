#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_FileSystem.H>
#include <typeinfo>
#include "Constants.H"
#include "SPADES.H"
#include "ConsoleIO.H"

/**
   @brief Main function
   @param argc [in] An integer argument count of the command line arguments
   @param argv [in] An argument vector of the command line arguments
   @return an integer 0 upon exit success
 **/
int main(int argc, char* argv[]) // NOLINT(bugprone-exception-escape)
{
#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#endif

    using namespace amrex::mpidatatypes;

    if (argc < 2) {
        // Print usage and exit with error code if no input file was provided.
        spades::io::print_usage(MPI_COMM_WORLD, std::cout);
        spades::io::print_error(
            MPI_COMM_WORLD, "No input file provided. Exiting!!");
        return 1;
    }

    // Look for "-h" or "--help" flag and print usage
    for (auto i = 1; i < argc; i++) {
        const std::string param(argv[i]);
        if ((param == "--help") || (param == "-h")) {
            spades::io::print_banner(MPI_COMM_WORLD, std::cout);
            spades::io::print_usage(MPI_COMM_WORLD, std::cout);
            return 0;
        }
    }

    if (!amrex::FileSystem::Exists(std::string(argv[1]))) {
        // Print usage and exit with error code if we cannot find the input file
        spades::io::print_usage(MPI_COMM_WORLD, std::cout);
        spades::io::print_error(
            MPI_COMM_WORLD, "Input file does not exist = " +
                                std::string(argv[1]) + ". Exiting!!");
        return 1;
    }

    spades::io::print_banner(MPI_COMM_WORLD, std::cout);

    amrex::Initialize(argc, argv);

    {
        BL_PROFILE("SPADES::main()");

        amrex::Real start_time = amrex::ParallelDescriptor::second();
        amrex::Print() << "Initializing SPADES ..." << std::endl;

        spades::SPADES spades_obj;
        spades_obj.init_data();

        amrex::Real init_time =
            amrex::ParallelDescriptor::second() - start_time;
        amrex::ParallelDescriptor::ReduceRealMax(
            init_time, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::Print() << "Initialization successful. Time elapsed = "
                       << init_time << std::endl;

        spades_obj.evolve();

        amrex::Real end_time = amrex::ParallelDescriptor::second() - start_time;
        amrex::ParallelDescriptor::ReduceRealMax(
            end_time, amrex::ParallelDescriptor::IOProcessorNumber());

        amrex::Print() << "\n==============================================="
                          "================================================="
                       << std::endl
                       << std::endl;
        amrex::Print() << "Time spent in init_data():    " << init_time
                       << std::endl;
        amrex::Print() << "Time spent in evolve():       "
                       << end_time - init_time << std::endl;
    }

    amrex::Finalize();
#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
