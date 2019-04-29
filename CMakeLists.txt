cmake_minimum_required(VERSION 3.10)
project(vmc)

find_package(MPI REQUIRED)

# Set compiler
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

# Compiling flags used
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -o - DNDEBUG")
SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0")


# Include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/Basis
    ${PROJECT_SOURCE_DIR}/src/Hamiltonians
    ${PROJECT_SOURCE_DIR}/src/InitialStates
    ${PROJECT_SOURCE_DIR}/src/InitialWeights
    ${PROJECT_SOURCE_DIR}/src/Metropolis
    ${PROJECT_SOURCE_DIR}/src/Optimization
    ${PROJECT_SOURCE_DIR}/src/Resampling
    ${PROJECT_SOURCE_DIR}/src/RNG
    ${PROJECT_SOURCE_DIR}/src/WaveFunctions
    ${PROJECT_SOURCE_DIR}/src/Eigen
    ${MPI_INCLUDE_PATH}
)

# Specify which functions to compile
set(SOURCES
    src/main.cpp
    src/system.cpp
    src/Hamiltonians/hamiltonian.cpp
    src/Hamiltonians/harmonicoscillator.cpp
    src/WaveFunctions/wavefunction.cpp
    src/InitialStates/initialstate.cpp
    src/InitialStates/randomuniform.cpp
    src/sampler.cpp
    src/InitialStates/randomnormal.cpp
    src/Hamiltonians/atomicnucleus.cpp
    src/InitialWeights/initialweights.cpp
    src/InitialWeights/randomize.cpp
    src/Metropolis/metropolis.cpp
    src/Metropolis/bruteforce.cpp
    src/Metropolis/importancesampling.cpp
    src/InitialWeights/constant.cpp
    src/Optimization/optimization.cpp
    src/WaveFunctions/gaussian.cpp
    src/WaveFunctions/partlyrestricted.cpp
    src/WaveFunctions/slaterdeterminant.cpp
    src/Optimization/barzilaiborwein.cpp
    src/RNG/rng.cpp
    src/RNG/mersennetwister.cpp
    src/RNG/parkmiller.cpp
    src/Optimization/gradientdescent.cpp
    src/Basis/basis.cpp
    src/Basis/hermite.cpp
    src/Basis/hydrogenorbital.cpp
    src/WaveFunctions/hydrogenlike.cpp
    src/Optimization/sgd.cpp
    src/Optimization/asgd.cpp
    src/Optimization/adam.cpp
    src/Basis/none.cpp
    src/WaveFunctions/padejastrow.cpp
    src/Resampling/blocker.cpp
    src/WaveFunctions/simplejastrow.cpp
    src/WaveFunctions/samsethjastrow.cpp
    src/WaveFunctions/padejastrow2.cpp
    src/WaveFunctions/rbmjastrow5.cpp
    src/WaveFunctions/rbmgaussian.cpp
    src/WaveFunctions/rbmgaussian2.cpp
    src/WaveFunctions/rbmjastrow.cpp
    src/WaveFunctions/rbmjastrow2.cpp
    src/WaveFunctions/rbmjastrow3.cpp
)


# Add executable file vmc
add_executable(vmc ${SOURCES})
target_link_libraries(vmc ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(vmc PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(vmc PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

