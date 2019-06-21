cmake_minimum_required(VERSION 3.10)
project(DotNet)

find_package(MPI REQUIRED)
find_package(GTest REQUIRED)

# Set compiler
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

# Compilation flags used
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -o -DNDEBUG")
SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

# Set default build configuration
IF( NOT CMAKE_BUILD_TYPE )
   SET(CMAKE_BUILD_TYPE Release)
ENDIF()

# Print which build configuration used
message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

# Include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/Basis
    ${PROJECT_SOURCE_DIR}/src/Hamiltonians
    ${PROJECT_SOURCE_DIR}/src/InitialStates
    ${PROJECT_SOURCE_DIR}/src/InitialWeights
    ${PROJECT_SOURCE_DIR}/src/Metropolis
    ${PROJECT_SOURCE_DIR}/src/Optimization
    ${PROJECT_SOURCE_DIR}/src/RNG
    ${PROJECT_SOURCE_DIR}/src/WaveFunctions
    ${PROJECT_SOURCE_DIR}/src/Eigen
    ${PROJECT_SOURCE_DIR}/src/blocker/c++
    ${PROJECT_SOURCE_DIR}/src/tqdm
    ${PROJECT_SOURCE_DIR}/test
    ${GTEST_INCLUDE_DIRS}
    ${MPI_INCLUDE_PATH}
)

# Specify which functions to compile
set(SOURCES
    src/main.cpp
    src/system.cpp
    src/Hamiltonians/hamiltonian.cpp
    src/Hamiltonians/harmonicoscillator.cpp
    src/Hamiltonians/doublewell.cpp
    src/WaveFunctions/wavefunction.cpp
    src/InitialStates/initialstate.cpp
    src/InitialStates/randomuniform.cpp
    src/sampler.cpp
    src/InitialStates/randomnormal.cpp
    src/Hamiltonians/atomicnucleus.cpp
    src/InitialWeights/initialweights.cpp
    src/InitialWeights/randomize.cpp
    src/InitialWeights/automatize.cpp
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
    src/Optimization/gradientdescent.cpp
    src/Basis/basis.cpp
    src/Basis/hermite.cpp
    src/Basis/hartreefock.cpp
    src/Basis/hydrogenorbital.cpp
    src/Basis/hermiteexpansion.cpp
    src/WaveFunctions/hydrogenlike.cpp
    src/Optimization/sgd.cpp
    src/Optimization/asgd.cpp
    src/Optimization/adam.cpp
    src/Basis/none.cpp
    src/WaveFunctions/padejastrow.cpp
    src/WaveFunctions/simplejastrow.cpp
    src/WaveFunctions/rbmgaussian.cpp
    src/WaveFunctions/rbmjastrow.cpp
    src/block/c++/blocker.cpp
)

set(TESTSOURCES 
    tests/main.cpp
    tests/test_basis.cpp
)

set(TEST_EXECUTABLE dotnetTests)
set(PROJECT_EXECUTABLE dotnet)


# Add executable file dotnet
add_executable(${PROJECT_EXECUTABLE} ${SOURCES})
target_link_libraries(${PROJECT_EXECUTABLE} ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(${PROJECT_EXECUTABLE} PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${PROJECT_EXECUTABLE} PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()


# Add test executable file
add_executable(${TEST_EXECUTABLE} ${TESTSOURCES})
target_link_libraries(${TEST_EXECUTABLE} ${GTEST_LIBRARIES} ${MPI_LIBRARIES})
add_test(${TEST_EXECUTABLE} ${TEST_EXECUTABLE})

if(MPI_COMPILE_FLAGS)
  set_target_properties(${TEST_EXECUTABLE} PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${TEST_EXECUTABLE} PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

