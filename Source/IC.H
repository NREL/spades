#ifndef IC_H
#define IC_H
#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include "Utilities.H"

/**
   @brief SPADES initial conditions
 **/
namespace spades::ic {

//! Constant initial condition
struct Constant
{
    struct DeviceOp
    {

        AMREX_GPU_DEVICE
        void operator()(
            const amrex::IntVect& /*iv*/,
            amrex::GeometryData const& /*geom*/,
            amrex::Real& state) const
        {
            state = 0.0;
        }
    };

    using DeviceType = DeviceOp;

    explicit Constant();

    static std::string identifier() { return "ic_constant"; }

    [[nodiscard]] DeviceType device_instance() const { return m_op; }

    DeviceOp m_op;
};

//! Base class for initial conditions
class InitializerBase
{
public:
    InitializerBase() = default;

    virtual ~InitializerBase() = default;

    virtual void initialize(const int lev, amrex::GeometryData const& geom) = 0;
};

//! Initial condition class
template <typename ICOp>
class Initializer : public InitializerBase
{
public:
    using ICDeviceOp = typename ICOp::DeviceType;

    Initializer(const ICOp ic_op, amrex::Vector<amrex::MultiFab>& state)
        : m_op(ic_op), m_state(state)
    {}

    ~Initializer() override = default;

    /**
       @brief Initialize m_state using the initial condition operator
       @param lev [in] level
       @param geom [in] geometry information
     **/
    void initialize(const int lev, amrex::GeometryData const& geom) override
    {
        auto const& s_arrs = m_state[lev].arrays();
        ICDeviceOp ic(m_op.device_instance());
        amrex::ParallelFor(
            m_state[lev], m_state[lev].nGrowVect(), constants::N_STATES,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
                const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                amrex::Real state = 0.0;
                ic(iv, geom, state);

                s_arrs[nbx](i, j, k, n) = state;
            });
        amrex::Gpu::synchronize();
    }

    //! Initial condition operator
    const ICOp m_op;

    //! State that will be filled with the initial condition operator
    amrex::Vector<amrex::MultiFab>& m_state;
};

} // namespace spades::ic
#endif
