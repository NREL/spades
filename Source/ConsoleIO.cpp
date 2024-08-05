#include <chrono>
#include <ctime>
#include "AMReX_Array.H"
#include "ConsoleIO.H"
#include "Source/SpadesVersion.H"

namespace amrex {
const char* buildInfoGetBuildDate();    // NOLINT(readability-identifier-naming)
const char* buildInfoGetComp();         // NOLINT(readability-identifier-naming)
const char* buildInfoGetGitHash(int i); // NOLINT(readability-identifier-naming)
const char* buildInfoGetCompVersion();  // NOLINT(readability-identifier-naming)
} // namespace amrex

namespace spades::io {

namespace {
const std::string DBL_LINE = std::string(78, '=') + "\n";
const std::string DASH_LINE = "\n" + std::string(78, '-') + "\n";
constexpr int TIME_BUF_SIZE = 64;
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
    amrex::Array<char, TIME_BUF_SIZE> time_buf;
    ctime_r(&etimet, time_buf.begin());
    const std::string tstamp(time_buf.begin());

    const std::string dirty_tag = (version::SPADES_DIRTY_REPO == "DIRTY")
                                      ? ("-" + version::SPADES_DIRTY_REPO)
                                      : "";
    const std::string spades_version = version::SPADES_VERSION + dirty_tag;
    const std::string spades_git_sha = version::SPADES_GIT_SHA + dirty_tag;

    // clang-format off
    out << DBL_LINE
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
        << "ON    (max threads = " << amrex::OpenMP::get_max_threads()
        << ")" << std::endl
#else
        << "OFF" << std::endl
#endif
        << std::endl;

    out << "  This software is released under the Apache 2.0 license.           "
        << std::endl
        << "  See https://github.com/NREL/spades/blob/main/LICENSE for details. "
        << DASH_LINE << std::endl;
    // clang-format on
}

} // namespace spades::io
