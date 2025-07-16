module load rocm/6 netcdf-c openblas git parmetis metis cmake

# Paths to repositories - fill in
export TRILINOS_SRC=$HOME/eigensolver/trilinos_src
mkdir -p $TRILINOS_SRC
cd $TRILINOS_SRC
git clone -b trilinos-release-16-1-0 https://github.com/trilinos/Trilinos.git

export TRILINOS_BUILD=$HOME/eigensolver/trilinos_build
export TRILINOS_INSTALL=$HOME/eigensolver/trilinos_install
export OMPI_CXX=$ROCM_PATH/bin/hipcc
mkdir -p $TRILINOS_BUILD
mkdir -p $TRILINOS_INSTALL

# Configure Trilinos using Kokkos and KokkosKernels as TPLs
cd $TRILINOS_BUILD
rm -rf CMake*
echo "Configuring Trilinos in directory: $PWD"
cmake \
  -B ${TRILINOS_BUILD} \
  -D CMAKE_BUILD_TYPE:STRING="RelWithDebInfo" \
  -D CMAKE_INSTALL_PREFIX:PATH=${TRILINOS_INSTALL} \
  -D CMAKE_EXE_LINKER_FLAGS="-L/opt/cray/pe/lib64/ -L${MPICH_DIR}/lib -lmpi ${CRAY_XPMEM_POST_LINK_OPTS} -lxpmem ${PE_MPICH_GTL_DIR_amd_gfx90a} ${PE_MPICH_GTL_LIBS_amd_gfx90a}" \
  -D CMAKE_SHARED_LINKER_FLAGS="-L/opt/cray/pe/lib64/ -L${MPICH_DIR}/lib -lmpi ${CRAY_XPMEM_POST_LINK_OPTS} -lxpmem ${PE_MPICH_GTL_DIR_amd_gfx90a} ${PE_MPICH_GTL_LIBS_amd_gfx90a}" \
  -D CMAKE_CXX_FLAGS="-I${MPICH_DIR}/include" \
  -D CMAKE_CXX_STANDARD=17 \
  -D CMAKE_CXX_COMPILER=hipcc \
  -D CMAKE_Fortran_FLAGS=-hnopattern \
  -D BUILD_SHARED_LIBS=OFF \
  -D TPL_ENABLE_MPI=ON \
  -D BLAS_LIBRARY_DIRS=${OLCF_OPENBLAS_ROOT}/lib \
  -D BLAS_LIBRARY_NAMES=openblas \
  -D Trilinos_ENABLE_Anasazi=ON \
  -D TPL_ENABLE_LAPACK=ON \
  -D LAPACK_LIBRARY_DIRS=${OLCF_OPENBLAS_ROOT}/lib \
  -D LAPACK_LIBRARY_NAMES=openblas \
  -D TPL_ENABLE_Boost=ON \
  -D TPL_ENABLE_BoostLib=OFF \
  -D TPL_ENABLE_Matio=OFF \
  -D TPL_ENABLE_METIS=ON \
  -D TPL_ENABLE_ParMETIS=OFF \
  -D TPL_ENABLE_PARDISO_MKL=OFF \
  -D TPL_ENABLE_ROCBLAS=ON \
  -D TPL_ENABLE_ROCSOLVER=ON \
  -D TPL_ENABLE_ROCSPARSE=ON \
  -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION=ON \
  -D Trilinos_ENABLE_ALL_PACKAGES=OFF \
  -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
  -D Trilinos_ENABLE_TESTS=OFF \
  -D Trilinos_ENABLE_EXAMPLES=OFF \
  -D Trilinos_ENABLE_OpenMP=ON \
  -D Trilinos_ENABLE_Kokkos=ON \
  -D Kokkos_ENABLE_SERIAL=ON \
  -D Kokkos_ENABLE_OPENMP=ON \
  -D Kokkos_ENABLE_HIP=ON \
  -D Kokkos_ARCH_AMD_GFX90A:BOOL=ON \
  -D Trilinos_ENABLE_Tpetra=ON \
  -D Tpetra_INST_COMPLEX_DOUBLE=OFF \
  -D Tpetra_INST_COMPLEX_FLOAT=OFF \
  -D Tpetra_INST_SERIAL=ON \
  -D Tpetra_INST_OPENMP=ON \
  -D Tpetra_INST_HIP=ON \
  -D Tpetra_ASSUME_GPU_AWARE_MPI=OFF \
  -D Trilinos_SHOW_DEPRECATED_WARNINGS=OFF \
  -D Panzer_ENABLE_EXAMPLES=OFF \
  -D Sacado_ENABLE_HIERARCHICAL_DFAD=ON \
  -D TPL_ENABLE_gtest=OFF \
  -D Trilinos_ENABLE_Gtest=OFF \
  -S ${TRILINOS_SRC}/Trilinos \

cmake --build ${TRILINOS_BUILD} -j 32 && cmake --install ${TRILINOS_BUILD}
export TPETRA_ASSUME_GPU_AWARE_MPI=0

cd ..
cmake -B build -DCMAKE_CXX_COMPILER=hipcc -DCMAKE_PREFIX_PATH=${TRILINOS_INSTALL}
cmake --build build
