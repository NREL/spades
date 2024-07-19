.. _Building:

Building
--------

`SPADES` is built using CMake. This is an example CMake configure command executed in the build directory::

  $ cd Build
  $ cmake -DCMAKE_BUILD_TYPE:STRING=Release \
          -DSPADES_ENABLE_MPI:BOOL=ON \
          -DCMAKE_CXX_COMPILER:STRING=mpicxx \
          -DCMAKE_C_COMPILER:STRING=mpicc \
          ..
  $ cmake --build . --parallel $(sysctl -n hw.ncpu)
