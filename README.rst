SPADES: Solver for PArallel Discrete Event Simulation
-----------------------------------------------------

|CI Badge| |Documentation Badge| |OpenSSF Badge| |SWR| |License Badge| |AMReX Badge| |C++ Badge|

.. |CI Badge| image:: https://github.com/NREL/spades/actions/workflows/ci.yml/badge.svg
   :target: https://github.com/NREL/spades/actions

.. |Documentation Badge| image:: https://github.com/NREL/spades/actions/workflows/docs.yml/badge.svg
   :target: https://nrel.github.io/spades/

.. |OpenSSF Badge| image:: https://www.bestpractices.dev/projects/11128/badge
   :target: https://www.bestpractices.dev/projects/11128

.. |License Badge| image:: https://img.shields.io/badge/License-Apache%20v2.0-blue.svg
   :target: https://www.apache.org/licenses/LICENSE-2.0

.. |SWR| image:: https://img.shields.io/badge/SWR-10.11578/dc.20250905.4-blue.svg
   :target: https://doi.org/10.11578/dc.20250905.4

.. |AMReX Badge| image:: https://img.shields.io/static/v1?label=%22powered%20by%22&message=%22AMReX%22&color=%22blue%22
   :target: https://amrex-codes.github.io/amrex/

.. |C++ Badge| image:: https://img.shields.io/badge/language-C%2B%2B17-blue
   :target: https://isocpp.org/

SPADES (Solver for PArallel Discrete Event Simulation) is an
open-source parallel discrete event simulation (PDES) package built on
the AMReX library. Targeted at solving discrete event systems in
parallel, this software package aims to be performance portable and
scalable on heterogeneous computing architectures, e.g., graphic
processing units (GPU). SPADES implements optimistic synchronization
with rollback through an implementation of the Time Warp algorithm. An
alternative conservative synchronization approach is also implemented
using the Lower Bound on Incoming Time Stamp. In our implementation,
logical processes are represented as cells in a grid and event
messages are represented as particles. SPADES supports various
parallel decomposition strategies, including the use of the Message
Passing Interface (MPI) and OpenMP threading. All major GPU
architectures (e.g., Intel, AMD, NVIDIA) are supported through the use
of performance portability functionalities implemented in AMReX. The
SPADES software is released in Department of Energy Software Record
`SWR-24-99 “SPADES (Scalable Parallel Discrete Events Simulation)”
<https://doi.org/10.11578/dc.20250905.4>`_.


Getting Started
~~~~~~~~~~~~~~~

To compile and run `SPADES`, one needs a C++ compiler that supports the C++17 standard, and then execute:

.. code-block:: console

    $ git clone --recursive git@github.com:NREL/spades.git
    $ cd Build
    $ ./cmake.sh
    $ ./spades example.inp

Dependencies
~~~~~~~~~~~~

`SPADES` is built on the `AMReX <https://github.com/AMReX-Codes/amrex>`_ library.


Documentation
~~~~~~~~~~~~~

The full documentation for `SPADES` exists in the Docs directory; at present this is maintained inline using `Sphinx <https://www.sphinx-doc.org/>`_ and `Doxygen <https://www.doxygen.nl/index.html>`_. To build the documentation:

.. code-block:: console

    $ cd Build && cmake -DSPADES_ENABLE_DOCUMENTATION:BOOL=ON .. && cmake --build . -t docs

Contributing, reporting bugs, and requesting help
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To report issues or bugs please `create a new issue <https://github.com/NREL/spades/issues>`_ on GitHub.

We welcome contributions from the community in form of bug fixes,
feature enhancements, documentation updates, etc. Please refer to the
`contributing guidelines
<https://github.com/NREL/spades/blob/main/CONTRIBUTING.md>`_ for more
information. All contributions are processed through pull-requests on
GitHub. Please refer to the `style guide
<https://nrel.github.io/spades/StyleGuide.html>`_ as a reference for
the best practices currently used to develop SPADES.

Please acknowledge as a publication co-author any developer that has
significantly contributed to implementing or improving specific
capability that was used for that publication.

Funding
~~~~~~~

This work was authored by the National Renewable Energy Laboratory (NREL) for the U.S. Department of Energy (DOE) under Contract No. DE-AC36-08GO28308. This work was supported by the Laboratory Directed Research and Development (LDRD) Program at NREL. The views expressed in the article do not necessarily represent the views of the DOE or the U.S. Government. The U.S. Government retains and the publisher, by accepting the article for publication, acknowledges that the U.S. Government retains a nonexclusive, paid-up, irrevocable, worldwide license to publish or reproduce the published form of this work, or allow others to do so, for U.S. Government purposes. A portion of the research was performed using computational resources sponsored by the Department of Energy’s Office of Energy Efficiency and Renewable Energy and located at the National Renewable Energy Laboratory.
