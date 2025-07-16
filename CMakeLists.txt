cmake_minimum_required(VERSION 3.23)

project(Eigensolver CXX)

# Look for an installed Trilinos
find_package(Trilinos REQUIRED)
find_package(MPI REQUIRED)

#SET(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")

#INCLUDE_DIRECTORIES (${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})

#LINK_DIRECTORIES (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})


add_executable(eigensolver eigensolver.cpp)

target_include_directories(eigensolver PRIVATE ${MPI_INCLUDE_PATH})

target_link_libraries(eigensolver PRIVATE
                      ${Trilinos_LIBRARIES}
                      ${Trilinos_TPL_LIBRARIES}
                      ${MPI_CXX_LIBRARIES})
