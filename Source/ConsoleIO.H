#ifndef CONSOLEIO_H
#define CONSOLEIO_H

#include <iostream>
#include <AMReX.H>

/**
   @brief IO functionality for SPADES
 **/
namespace spades::io {

/**
   @brief Print SPADES banner
   @param comm [in] MPI communicator
   @param out [out] output stream
 **/
void print_banner(MPI_Comm comm, std::ostream& out);

/**
   @brief Print SPADES usage
   @param comm [in] MPI communicator
   @param out [out] output stream
 **/
void print_usage(MPI_Comm comm, std::ostream& out);

/**
   @brief Print SPADES error
   @param comm [in] MPI communicator
   @param msg [in] output message
 **/
void print_error(MPI_Comm comm, const std::string& msg);

} // namespace spades::io

#endif
