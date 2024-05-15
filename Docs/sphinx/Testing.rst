#.. _Testing:

Testing
-------

Testing of `SPADES` can be performed using CTest, which is included in the CMake build system. With CMake, this can be enabled with the following options ::

  $ cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
          -DCMAKE_CXX_COMPILER:STRING=clang++ \
          -DCMAKE_C_COMPILER:STRING=clang \
          -DCMAKE_BUILD_TYPE:STRING=Release \
          -DSPADES_DIM:STRING=3 \
          -DSPADES_ENABLE_MPI:BOOL=OFF \
          -DSPADES_ENABLE_CPPCHECK:BOOL=OFF \
          -DSPADES_ENABLE_CLANG_TIDY:BOOL=OFF \
          ..

To perform the tests and compare to previously generated gold files, use the following additional options::

  -DSPADES_TEST_WITH_FCOMPARE:BOOL=ON \
  -DSPADES_SAVED_GOLDS_DIRECTORY:STRING=$(pwd)/golds/tmp \
  -DSPADES_REFERENCE_GOLDS_DIRECTORY:STRING=$(pwd)/golds/current \
