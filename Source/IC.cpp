#include "IC.H"
namespace spades::ic {

Constant::Constant() { amrex::ParmParse pp(identifier()); }

} // namespace spades::ic
