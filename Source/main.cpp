#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <typeinfo>
#include "Constants.H"
#include "SPADES.H"

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    BL_PROFILE("SPADES::main()");

    {
        spades::SPADES spades_obj;

        spades_obj.init_data();

        spades_obj.evolve();
    }

    amrex::Finalize();
}
