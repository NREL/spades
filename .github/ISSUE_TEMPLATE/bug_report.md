---
name: Bug report
about: Report a bug to help us improve
title: 'Bug report'
labels: "bug:spades"
---

## Bug description
<!-- A clear and concise description of the bug. -->

## Steps to reproduce
<!-- Update the following list with your specific information. -->

Steps to reproduce the behavior:
1. Compiler used
   - [ ] GCC
   - [ ] LLVM
   - [ ] oneapi (Intel)
   - [ ] nvcc (NVIDIA)
   - [ ] rocm (AMD)
   - [ ] with MPI
   - [ ] other:
2. Operating system
   - [ ] Linux
   - [ ] OSX
   - [ ] Windows
   - [ ] other (do tell ;)):
3. Hardware:
   - [ ] CPU
   - [ ] GPU
4. Machine details ():
```
<!-- name, modules loaded, environment variables, etc -->
```
5. Input file attachments <!-- Please upload the input files in a zip or point to a public branch. -->
6. Error (paste or attach):
```
<!-- error output -->
```
7. If this is a segfault, a stack trace from a debug build (paste or attach):
```
<!-- stack trace -->
```

## Expected behavior
<!-- A clear and concise description of what is expected behavior. -->

## SPADES information
<!-- Please provide as much detail as possible including git commit. The best information is a snapshot of the SPADES header. -->

```
==============================================================================
                SPADES (https://github.com/NREL/spades)

  SPADES version :: 2306055-DIRTY
  SPADES Git SHA :: 2306055a77c48659db623b835a5817677be5fb47-DIRTY
  AMReX version  :: 24.04-6-g771c439170bf

  Exec. time     :: Mon May 13 14:46:55 2024
  Build time     :: May 13 2024 14:46:51
  C++ compiler   :: Clang 18.1.5

  MPI            :: OFF
  GPU            :: OFF
  OpenMP         :: OFF

 This software is released under the Apache 2.0 license.
 See https://github.com/NREL/spades/blob/main/LICENSE for details.
------------------------------------------------------------------------------
```

## Additional context
<!-- Screenshots, related issues, etc -->
