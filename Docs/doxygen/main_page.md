# SPADES API documentation {#mainpage}

[SPADES](https://github.com/NREL/spades) is an open-source parallel
discrete events simulation (PDES) package built on the AMReX
library. Targeted at solving discrete events systems in parallel, this
software package aims to be performance portable and scalable on
heterogeneous computing architectures, e.g., graphic processing units
(GPU).

This document is intended for developers who want to understand the C++ code
structure and modify the codebase. It assumes that the reader is
familiar with the installation, compilation, and execution steps.

## How to use this API guide?

The API documentation is automatically generated from specially-formatted
comments in the source code using [Doxygen](https://www.doxygen.nl/index.html).

### Source code organization

After cloning the SPADES repository, the project directory (`spades`) contains the following:

- `Source` -- C++ source files
- `Build` -- Example build directory with a sample cmake script
- `CMake` -- Functions, utilities used during CMake configuration phase
- `Docs` -- Documentation (Sphinx-based and Doxygen)
- `Submodules` -- Third-party libraries and dependencies
- `Tests` -- Regression tests and associated infrastructure

### Building API documentation locally

To generate this documentation on a local machine, or to rebuild docs
during code development, `doxygen` and `graphviz` are required. 
```{shell}
$ cd Build && cmake -B build-docs ../Docs && cmake --build build-docs
```
The resulting documentation is in `build-docs/doxygen/html` directory.
