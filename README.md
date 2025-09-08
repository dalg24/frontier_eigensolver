# Solving Eigenvalue Problems on OLCF's Frontier System with Anasazi

This repository provides a small demonstration program showcasing how to use
the [Anasazi](https://trilinos.github.io/anasazi.html) package, part of the
Trilinos project, to solve eigenvalue problems on the Frontier supercomputer at
the Oak Ridge Leadership Computing Facility (OLCF).

Anasazi is a robust and flexible library for solving large-scale eigenvalue
problems. This example focuses on demonstrating its basic usage within a
GPU-accelerated environment on Frontier.

## Using Pre-installed Software on Frontier

This is the recommended approach for most users, leveraging the optimized
software stack provided on Frontier.

### Download the Source

First, clone this repository to your desired location:

```bash
git clone https://github.com/dalg24/frontier_eigensolver.git eigensolver
cd eigensolver
```

### Load Modules

Load the necessary modules for AMD GPUs, ROCm, Cray MPICH, and the
pre-installed GPU-enabled Trilinos. The specific versions below are confirmed
as of 2025/07/16, but always check for the latest recommended versions on
Frontier.

```bash
module load amd/6.3.1 rocm/6.3.1 cray-mpich/8.1.32
module load trilinos/16.0.0-gpu-mpi
```

### Configure and Build

Use CMake to configure and build the project. We specify hipcc as the C++
compiler to ensure GPU offloading is enabled.

```bash
cmake -B builddir -DCMAKE_CXX_COMPILER=hipcc
cmake --build builddir
```

### Run the Example

Execute the compiled program in help mode. This will print the commandline options

```bash
./builddir/eigensolver --help
```

Per default it computes the 3 lowest Evals of the matrix in matrix_export.mtx via a BlockKrylovSchur solver up to machine precision.

You should see output indicating the progress of the eigenvalue solver and the
computed eigenvalues.

## Building Trilinos from Source (Optional)
For users who require more control over the versions of Kokkos/Trilinos,
specific configuration options, or development builds, a helper script
`build_frontier.sh` is provided. This script automates the process of building
Trilinos from source with appropriate settings for Frontier.

To use it:

```
./build_frontier.sh
```

**Note:** Building Trilinos from source can be time-consuming and requires
significant disk space. Refer to the comments within `build_frontier.sh` for
specific details and potential modifications.
