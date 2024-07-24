SPADES: Scalable PArallel Discrete Event Simulation
---------------------------------------------------

|CI Badge| |Documentation Badge| |License Badge| |AMReX Badge| |C++ Badge|

.. |CI Badge| image:: https://github.com/NREL/spades/workflows/SPADES-CI/badge.svg
   :target: https://github.com/NREL/spades/actions

.. |Documentation Badge| image:: https://github.com/NREL/spades/workflows/SPADES-Docs/badge.svg
   :target: https://https://nrel.github.io/spades

.. |License Badge| image:: https://img.shields.io/badge/License-Apache%20v2.0-blue.svg
   :target: https://www.apache.org/licenses/LICENSE-2.0

.. |AMReX Badge| image:: https://img.shields.io/static/v1?label=%22powered%20by%22&message=%22AMReX%22&color=%22blue%22
   :target: https://amrex-codes.github.io/amrex/

.. |C++ Badge| image:: https://img.shields.io/badge/language-C%2B%2B17-blue
   :target: https://isocpp.org/

SPADES (Scalable PArallel Discrete Event Simulation) is an open-source
parallel discrete event simulation (PDES) package built on the AMReX
library. Targeted at solving discrete event systems in parallel, this
software package aims to be performance portable and scalable on
heterogeneous computing architectures, e.g., graphic processing units
(GPU). SPADES implements optimistic synchronization with rollback
through an implementation of the Time Warp algorithm. An alternative
conservative synchronization approach is also implemented using the
Lower Bound on Incoming Time Stamp. In our implementation, logical
processes are represented as cells in a grid and event messages are
represented as particles. SPADES supports various parallel
decomposition strategies, including the use of the Message Passing
Interface (MPI) and OpenMP threading. All major GPU architectures
(e.g., Intel, AMD, NVIDIA) are supported through the use of
performance portability functionalities implemented in AMReX.

Getting Started
~~~~~~~~~~~~~~~

To compile and run `SPADES`, one needs a C++ compiler that supports the C++17 standard, and then execute ::

    $ git clone --recursive git@github.com:NREL/spades.git
    $ cd Build
    $ ./cmake.sh
    $ ./spades example.inp

Dependencies
~~~~~~~~~~~~

`SPADES` is built on the `AMReX <https://github.com/AMReX-Codes/amrex>`_ library.


Documentation
~~~~~~~~~~~~~

The full documentation for `SPADES` exists in the Docs directory; at present this is maintained inline using `Sphinx <https://www.sphinx-doc.org/>`_ and `Doxygen <https://www.doxygen.nl/index.html>`_. To build the documentation ::

    $ cd Build && cmake -B build-docs ../Docs && cmake --build build-docs
