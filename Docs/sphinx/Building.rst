.. _Building:

Building
--------

`SPADES` is built using CMake.

CMake
~~~~~

To use the CMake option, one executes the CMake configure command in the build directory::

  cd Build
  cmake -DCMAKE_BUILD_TYPE:STRING=Release \
        -DSPADES_ENABLE_MPI:BOOL=ON \
        -DCMAKE_CXX_COMPILER:STRING=mpicxx \
        -DCMAKE_C_COMPILER:STRING=mpicc \
        .. && make
