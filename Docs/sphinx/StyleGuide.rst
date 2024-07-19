.. _StyleGuide:

Style Guide
-----------

For consistency, we enforce some style guidelines with `ClangFormat
<https://clang.llvm.org/docs/ClangFormat.html>`_ and `Clang-Tidy
<https://clang.llvm.org/extra/clang-tidy/>`_. These are defined in the
:file:`.clang-format` and :file:`.clang-tidy` files.

Here are some conventions to follow:

#. All `SPADES` specific code must be within ``spades`` namespace.

#. Following AMReX convention, header files will use a ``.H`` extension and C++
   source files will use a ``.cpp`` extension.

#. Use `snake_case <https://en.wikipedia.org/wiki/Snake_case>`_ almost
   everywhere, i.e., for variables, namespaces, function and method names.
   Exceptions include when overriding AMReX class methods in inherited classes.

#. Use `CamelCase <https://en.wikipedia.org/wiki/Camel_case>`_ for class names.
   Capitalize the first letter of the first word in the compound name also.

#. Use ``m_`` prefix for class instance variables. No prefix for class methods.

#. Keep logic in functions short. Always use descriptive names for function
   names and variables that provide the reader with clues as to what the
   function does.

Formatting
``````````

This is a manual way to format all source files::

  $ cd Source
  $ find . \( -name "*.cpp" -o -name "*.H" \) -exec clang-format -i {} +


Linting
```````

Linting of `SPADES` can be performed using cmake options to configure
cppcheck and clang-tidy. For example, the continuous integration tools
through Github Actions uses the following for clang-tidy ::

  $ cmake -DCMAKE_BUILD_TYPE:STRING=Debug \
          -DCMAKE_CXX_COMPILER:STRING=clang++ \
          -DCMAKE_C_COMPILER:STRING=clang \
          -DSPADES_DIM:STRING=3 \
          -DSPADES_ENABLE_MPI:BOOL=OFF \
          -DSPADES_TEST_WITH_FCOMPARE:BOOL=OFF \
          -DSPADES_ENABLE_ALL_WARNINGS:BOOL=ON \
          -DSPADES_ENABLE_CPPCHECK:BOOL=OFF \
          -DSPADES_ENABLE_CLANG_TIDY:BOOL=ON \
          ..
  $ cmake --build . | tee -a clang-tidy-full-report.txt

A code spellchecker can be run by running `codespell
<https://github.com/codespell-project/codespell>`_ in the project
directory.
