#ifndef CONSOLEIO_H
#define CONSOLEIO_H

#include <iostream>
#include <AMReX.H>

namespace spades::io {

void print_banner(MPI_Comm /*comm*/, std::ostream& /*out*/);

void print_usage(MPI_Comm /*comm*/, std::ostream& /*out*/);

void print_error(MPI_Comm /*comm*/, const std::string& msg);

void print_summary(std::ostream&);

} // namespace spades::io

#endif /* CONSOLEIO_H */
