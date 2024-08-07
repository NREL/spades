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

.. cmakeval:: SPADES_DIM

   Number of physical dimensions. Default: "2"

.. cmakeval:: SPADES_ENABLE_ALL_WARNINGS

   Show most warnings for most compilers. Default: ON
.. cmakeval:: SPADES_ENABLE_CLANG_TIDY

   Compile with clang-tidy static analysis. Default: OFF

.. cmakeval:: SPADES_ENABLE_CPPCHECK

   Enable cppcheck static analysis target. Default: OFF

.. cmakeval:: SPADES_ENABLE_FCOMPARE

   Enable building fcompare when not testing. Default: OFF

.. cmakeval:: SPADES_TEST_WITH_FCOMPARE

   Check test plots against gold files. Default: OFF

.. cmakeval:: SPADES_SAVE_GOLDS

   Provide a directory in which to save golds during testing. Default: OFF

.. cmakeval:: SPADES_ENABLE_FPE_TRAP_FOR_TESTS

   Enable FPE trapping in tests. Default: ON

.. cmakeval:: SPADES_ENABLE_MPI

   Enable MPI. Default: OFF

.. cmakeval:: SPADES_ENABLE_OPENMP

   Enable OpenMP. Default: OFF

.. cmakeval:: SPADES_ENABLE_CUDA

   Enable CUDA. Default: OFF

.. cmakeval:: SPADES_ENABLE_ROCM

   Enable ROCm/HIP. Default: OFF

.. cmakeval:: SPADES_ENABLE_SYCL

   Enable Intel OneAPI SyCL. Default: OFF

.. cmakeval:: SPADES_ENABLE_TINY_PROFILE

   Enable AMReX TinyProfile support. Default: ON

.. cmakeval:: SPADES_ENABLE_DOCUMENTATION

   Enable documentation target. Default: OFF

.. cmakeval:: SPADES_ENABLE_ASAN

   Enable AddressSanitizer. Default: OFF

.. cmakeval:: SPADES_ENABLE_LSAN

   Enable LeakSanitizer. Default: OFF

.. cmakeval:: SPADES_ENABLE_UBSAN

   Enable UndefinedBehaviorSanitizer. Default: OFF

.. cmakeval:: SPADES_ENABLE_TSAN

   Enable ThreadSanitizer. Default: OFF

.. cmakeval:: SPADES_PRECISION DOUBLE.

   Floating point precision SINGLE or DOUBLE. Default: "DOUBLE"
