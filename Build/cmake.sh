#!/bin/bash

function clean(){
    rm -f CMakeCache.txt
    rm -f CTestCustom.cmake
    rm -f CTestTestfile.cmake
    rm -f DartConfiguration.tcl
    rm -f Makefile
    rm -f SpadesConfig.cmake
    rm -f spades
    rm -f cmake_install.cmake
    rm -f compile_commands.json
    rm -rf CMakeDoxyfile.in
    rm -rf CMakeDoxygenDefaults.cmake
    rm -rf CMakeFiles
    rm -rf Source
    rm -rf Submodules
    rm -rf Testing
    rm -rf Tests
    rm -rf Docs
    rm -rf spades_obj
    rm -rf clang-tidy-full-report.txt
    rm -rf clang-tidy-ci-report.txt
    rm -rf cppcheck-ci-report.txt
    rm -rf cppcheck
    rm -rf bin
    rm -rf ../Docs/sphinx/doxygen
}

CLEAN='false'
TEST='false'

while getopts :ct flag
do
    case "${flag}" in
        c)
            CLEAN='true'
            ;;
        t)
            TEST='true'
            ;;
        '?')
            echo "INVALID OPTION -- ${OPTARG}" >&2
            exit 1
            ;;
        ':')
            echo "MISSING ARGUMENT for option -- ${OPTARG}" >&2
            exit 1
            ;;
        *)
            echo "UNIMPLEMENTED OPTION -- ${flag}" >&2
            exit 1
            ;;
    esac
done

if ${CLEAN}; then
    echo "Cleaning build"
    clean
    exit
fi

cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DCMAKE_CXX_COMPILER:STRING=clang++ \
      -DCMAKE_C_COMPILER:STRING=clang \
      -DSPADES_ENABLE_MPI:BOOL=OFF \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DSPADES_DIM:STRING=2 \
      -DSPADES_TEST_WITH_FCOMPARE:BOOL=OFF \
      -DSPADES_TEST_WITH_PYTHON:BOOL=OFF \
      -DSPADES_SAVE_GOLDS:BOOL=OFF \
      -DSPADES_SAVED_GOLDS_DIRECTORY:STRING=$(pwd)/golds/tmp \
      -DSPADES_REFERENCE_GOLDS_DIRECTORY:STRING=$(pwd)/golds/current \
      -DSPADES_ENABLE_FPE_TRAP_FOR_TESTS:BOOL=OFF \
      -DSPADES_ENABLE_CPPCHECK:BOOL=OFF \
      -DSPADES_ENABLE_CLANG_TIDY:BOOL=OFF \
      -DSPADES_ENABLE_DOCUMENTATION:BOOL=OFF \
      -DSPADES_ENABLE_OPENMP:BOOL=OFF \
      -DSPADES_ENABLE_ASAN:BOOL=OFF \
      -DSPADES_ENABLE_TSAN:BOOL=ON \
      ..



nice cmake --build . --parallel $(sysctl -n hw.ncpu)

# # Docs
# cmake --build . -t docs

if ${TEST}; then
    ctest -j 8
fi
