.. _Building:

Building
--------

`SPADES` is built using CMake. This is an example CMake configure command executed in the build directory:

.. code-block:: console

    $ cd Build
    $ cmake -DCMAKE_BUILD_TYPE:STRING=Release \
            -DSPADES_ENABLE_MPI:BOOL=ON \
            -DCMAKE_CXX_COMPILER:STRING=mpicxx \
            -DCMAKE_C_COMPILER:STRING=mpicc \
            ..
    $ cmake --build . --parallel $(sysctl -n hw.ncpu)

CMake configuration reference
`````````````````````````````

.. cmakeval:: SPADES_DIM "2" CACHE STRING "Number of physical dimensions")
.. cmakeval:: SPADES_ENABLE_ALL_WARNINGS "Show most warnings for most compilers" ON)
.. cmakeval:: SPADES_ENABLE_CLANG_TIDY "Compile with clang-tidy static analysis" OFF)
.. cmakeval:: SPADES_ENABLE_CPPCHECK "Enable cppcheck static analysis target" OFF)
.. cmakeval:: SPADES_ENABLE_FCOMPARE "Enable building fcompare when not testing" OFF)
.. cmakeval:: SPADES_TEST_WITH_FCOMPARE "Check test plots against gold files" OFF)
.. cmakeval:: SPADES_SAVE_GOLDS "Provide a directory in which to save golds during testing" OFF)
.. cmakeval:: SPADES_ENABLE_FPE_TRAP_FOR_TESTS "Enable FPE trapping in tests" ON)
.. cmakeval:: SPADES_ENABLE_MPI "Enable MPI" OFF)
.. cmakeval:: SPADES_ENABLE_OPENMP "Enable OpenMP" OFF)
.. cmakeval:: SPADES_ENABLE_CUDA "Enable CUDA" OFF)
.. cmakeval:: SPADES_ENABLE_ROCM "Enable ROCm/HIP" OFF)
.. cmakeval:: SPADES_ENABLE_SYCL "Enable Intel OneAPI SyCL" OFF)
.. cmakeval:: SPADES_ENABLE_TINY_PROFILE "Enable AMReX TinyProfile support" ON)
.. cmakeval:: SPADES_ENABLE_DOCUMENTATION "Enable documentation target" OFF)
.. cmakeval:: SPADES_PRECISION "DOUBLE" CACHE STRING "Floating point precision SINGLE or DOUBLE")
.. cmakeval:: SPADES_ENABLE_HDF5 "Enable HDF5 library" OFF)
.. cmakeval:: SPADES_ENABLE_HDF5_ZFP "Enable ZFP compression in HDF5 library" OFF)

