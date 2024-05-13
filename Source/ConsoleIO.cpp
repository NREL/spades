#include <chrono>
#include <ctime>
#include "AMReX_Array.H"
#include "ConsoleIO.H"
#include "Source/SpadesVersion.H"

namespace amrex {
const char* buildInfoGetBuildDate();
const char* buildInfoGetComp();
const char* buildInfoGetGitHash(int i);
const char* buildInfoGetCompVersion();
} // namespace amrex

namespace spades::io {

namespace {
const std::string dbl_line = std::string(78, '=') + "\n";
const std::string dash_line = "\n" + std::string(78, '-') + "\n";
} // namespace

void print_usage(MPI_Comm comm, std::ostream& out)
{
#ifdef AMREX_USE_MPI
    int irank = 0;
    int num_ranks = 1;
    MPI_Comm_size(comm, &num_ranks);
    MPI_Comm_rank(comm, &irank);

    // Only root process does the printing
    if (irank != 0) {
        return;
    }
#else
    amrex::ignore_unused(comm);
#endif

    out << R"doc(Usage:
    spades <input_file> [param=value] [param=value] ...

Required:
    input_file   : Input file with simulation settings

Optional:
    param=value  : Overrides for parameters during runtime
)doc" << std::endl;
}

void print_error(MPI_Comm comm, const std::string& msg)
{
#ifdef AMREX_USE_MPI
    int irank = 0;
    int num_ranks = 1;
    MPI_Comm_size(comm, &num_ranks);
    MPI_Comm_rank(comm, &irank);

    // Only root process does the printing
    if (irank != 0) {
        return;
    }
#else
    amrex::ignore_unused(comm);
#endif

    std::cout << "ERROR: " << msg << std::endl;
}

void print_banner(MPI_Comm comm, std::ostream& out)
{
#ifdef AMREX_USE_MPI
    int irank = 0;
    int num_ranks = 1;
    MPI_Comm_size(comm, &num_ranks);
    MPI_Comm_rank(comm, &irank);

    // Only root process does the printing
    if (irank != 0) {
        return;
    }
#else
    amrex::ignore_unused(comm);
#endif

    auto etime = std::chrono::system_clock::now();
    auto etimet = std::chrono::system_clock::to_time_t(etime);
    amrex::Array<char, 64> time_buf;
    ctime_r(&etimet, time_buf.begin());
    const std::string tstamp(time_buf.begin());

    const std::string dirty_tag = (version::spades_dirty_repo == "DIRTY")
                                      ? ("-" + version::spades_dirty_repo)
                                      : "";
    const std::string spades_version = version::spades_version + dirty_tag;
    const std::string spades_git_sha = version::spades_git_sha + dirty_tag;

    // clang-format off
    out << dbl_line
        << "                SPADES (https://github.com/NREL/spades)"
        << std::endl << std::endl
        << "  SPADES version :: " << spades_version << std::endl
        << "  SPADES Git SHA :: " << spades_git_sha << std::endl
        << "  AMReX version  :: " << amrex::Version() << std::endl << std::endl
        << "  Exec. time     :: " << tstamp
        << "  Build time     :: " << amrex::buildInfoGetBuildDate() << std::endl
        << "  C++ compiler   :: " << amrex::buildInfoGetComp()
        << " " << amrex::buildInfoGetCompVersion() << std::endl << std::endl
        << "  MPI            :: "
#ifdef AMREX_USE_MPI
        << "ON    (Num. ranks = " << num_ranks << ")" << std::endl
#else
        << "OFF " << std::endl
#endif
        << "  GPU            :: "
#ifdef AMREX_USE_GPU
        << "ON    "
#if defined(AMREX_USE_CUDA)
        << "(Backend: CUDA)"
#elif defined(AMREX_USE_HIP)
        << "(Backend: HIP)"
#elif defined(AMREX_USE_SYCL)
        << "(Backend: SYCL)"
#endif
        << std::endl
#else
        << "OFF" << std::endl
#endif
        << "  OpenMP         :: "
#ifdef AMREX_USE_OMP
        << "ON    (Num. threads = " << omp_get_max_threads() << ")" << std::endl
#else
        << "OFF" << std::endl
#endif
        << std::endl;

    out << " This software is released under the Apache 2.0 license.           "
        << std::endl
        << " See https://github.com/NREL/spades/blob/main/LICENSE for details. "
        << dash_line << std::endl;
    // clang-format on
}

} // namespace spades::io
