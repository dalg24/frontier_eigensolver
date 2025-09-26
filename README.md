# Solving Eigenvalue Problems on OLCF's Frontier System with Anasazi

This repository provides a small demonstration program showcasing how to use
the [Anasazi](https://trilinos.github.io/anasazi.html) package, part of the
Trilinos project, to solve eigenvalue problems on the Frontier supercomputer at
the Oak Ridge Leadership Computing Facility (OLCF).

Anasazi is a robust and flexible library for solving large-scale eigenvalue
problems. This example focuses on demonstrating its basic usage within a
GPU-accelerated environment on Frontier.

## Building Trilinos from Source (Required for support of complex numbers, PREFERRED)

### Download the Source

First, clone this repository to your desired location:

```bash
git clone https://github.com/dalg24/frontier_eigensolver.git eigensolver
cd eigensolver
```
### Building Trilinos
The repository supplies a helper script
`build_frontier.sh`. This script automates the process of building
Trilinos from source with appropriate settings for Frontier including support for complex numbers.

To use it:

```
module load rocm/6 netcdf-c openblas git parmetis metis cmake
./build_frontier.sh
```

**Note:** Building Trilinos from source can be time-consuming and requires
significant disk space. Refer to the comments within `build_frontier.sh` for
specific details and potential modifications.

### Configure and Build

Use CMake to configure and build the project. We specify hipcc as the C++
compiler to ensure GPU offloading is enabled.

```bash
cmake -B builddir -DCMAKE_CXX_COMPILER=hipcc -DCMAKE_PREFIX_PATH=$HOME/eigensolver/trilinos_install
cmake --build builddir
```

## Building Trilinos using spack (Works for complex numbers, takes a lot of time and diskspace)

One simple way to install the dependencies is to use the spack package manager.
Since compiling Trilinos will take considerable disk space, go to a location that has at least >50GB

```bash
git clone --depth=2 --branch=releases/v1.0 https://github.com/spack/spack.git ./spack
cd spack
. share/spack/setup-env.sh
```

Next we will create a custom environment for the eigensolver (to not accidentially pollute the general environment)

```bash
spack env create frontier_eigensolver
spack env activate frontier_eigensolver
```

Once the environment is activated, we install and load trilinos with gpu and complex number support

```bash
spack install --add trilinos@16.1.0 +rocm amdgpu_target=gfx90a +complex
spack load hipcc
spack load trilinos
```

Then we can clone the eigensolver with the following steps
```bash
git clone https://github.com/dalg24/frontier_eigensolver.git eigensolver
cd eigensolver
```

And build it using CMake
```bash
cmake -B builddir -DCMAKE_CXX_COMPILER=hipcc
cmake --build builddir
```

## Using Pre-installed Software on Frontier (currently not working with complex valued matrices)

This is the recommended approach for most users, leveraging the optimized
software stack provided on Frontier. Unfortunately the current software stack does not include Trilinos with enabled complex support


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

## Run the Example

Execute the compiled program in help mode. This will print the commandline options

```bash
./builddir/eigensolver --help
```

Per default it computes the 3 lowest Evals of the matrix in `matrix_export.mtx` via a BlockKrylovSchur solver up to machine precision.

You should see output indicating the progress of the eigenvalue solver and the
computed eigenvalues.
Useful calculations will look similar to this:

```bash
./builddir/eigensolver --nev=10 --filename=mymatrixfile.mtx --which=LR --solver=BlockDavidson
```
