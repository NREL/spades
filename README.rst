SPADES: Scalable PArallel Discrete Events Simulation
----------------------------------------------------

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

The full documentation for `SPADES` exists in the Docs directory; at present this is maintained inline using `Sphinx <http://www.sphinx-doc.org>`_. To build the documentation ::

    $ cd Docs && mkdir build && cd build && sphinx-build -M html ../sphinx .

Or, using cmake ::

    $ cd Build && cmake -B build-docs ../Docs && cmake --build build-docs
