Input Files and Controls
------------------------

The input file specified on the command line is a free-format text file, one entry per row, that specifies input data processed by the AMReX ``ParmParse`` module.

This file needs to specified along with the executable as an `argv` option, for example::

  mpirun -np 64 ./spades inputs


Entries can be overwritten on the command line: ``./spades inputs amr.plot_int=10``.

Example input files can be found in the ``Tests/test_files`` directory of this repository.
