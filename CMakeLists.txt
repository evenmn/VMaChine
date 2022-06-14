cmake_minimum_required(VERSION 3.10)
project(VMaChine)

find_package(MPI REQUIRED)
# find_package(GTest REQUIRED)

# Set compiler
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

# Compilation flags used
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -o -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")

# Set default build configuration
if( NOT CMAKE_BUILD_TYPE )
   SET(CMAKE_BUILD_TYPE Release)
endif()

# Include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/Activation
    ${PROJECT_SOURCE_DIR}/src/Basis
    ${PROJECT_SOURCE_DIR}/src/Hamiltonians
    ${PROJECT_SOURCE_DIR}/src/InitialStates
    ${PROJECT_SOURCE_DIR}/src/InitialWeights
    ${PROJECT_SOURCE_DIR}/src/Interaction
    ${PROJECT_SOURCE_DIR}/src/Layer
    ${PROJECT_SOURCE_DIR}/src/Metropolis
    ${PROJECT_SOURCE_DIR}/src/Optimization
    ${PROJECT_SOURCE_DIR}/src/RNG
    ${PROJECT_SOURCE_DIR}/src/WaveFunctions
#    ${PROJECT_SOURCE_DIR}/src/Eigen
    ${PROJECT_SOURCE_DIR}/include/Eigen
    ${PROJECT_SOURCE_DIR}/include/block
#    ${PROJECT_SOURCE_DIR}/src/blocker/c++
#    ${PROJECT_SOURCE_DIR}/test
#	   ${GTEST_INCLUDE_DIRS}
    ${MPI_INCLUDE_PATH}
)

# Specify which functions to compile
set(SOURCES
    src/main.cpp
    src/system.cpp
    src/Activation/activation.cpp
    src/Activation/elu.cpp
    src/Activation/leakyrelu.cpp
    src/Activation/relu.cpp
    src/Activation/purelinear.cpp
    src/Activation/sigmoid.cpp
    src/Hamiltonians/hamiltonian.cpp
    src/Hamiltonians/harmonicoscillator.cpp
    src/Hamiltonians/ellipticalharmonicoscillator.cpp
    src/Hamiltonians/doublewell.cpp
    src/WaveFunctions/wavefunction.cpp
    src/InitialStates/initialstate.cpp
    src/InitialStates/randomuniform.cpp
    src/sampler.cpp
    src/InitialStates/randomnormal.cpp
    src/Hamiltonians/atomicnucleus.cpp
    src/InitialWeights/initialweights.cpp
    src/InitialWeights/randomuniform.cpp
    src/InitialWeights/randomnormal.cpp
    src/InitialWeights/automatize.cpp
    src/InitialWeights/customized.cpp
    src/Interaction/interaction.cpp
    src/Interaction/nointeraction.cpp
    src/Interaction/coulomb.cpp
    src/Layer/layer.cpp
    src/Layer/dense.cpp
    src/Layer/input.cpp
    src/Layer/output.cpp
    src/Metropolis/metropolis.cpp
    src/Metropolis/bruteforce.cpp
    src/Metropolis/importancesampling.cpp
    src/InitialWeights/constant.cpp
    src/Optimization/optimization.cpp
    src/WaveFunctions/gaussian.cpp
    src/WaveFunctions/partlyrestricted.cpp
    src/WaveFunctions/slaterdeterminant.cpp
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
    src/WaveFunctions/hardcorejastrow.cpp
    src/WaveFunctions/rbmgaussian.cpp
    src/WaveFunctions/rbmproduct.cpp
    src/WaveFunctions/fnn.cpp
    src/block/c++/blocker.cpp
)

set(TESTSOURCES
    tests/main.cpp
    tests/test_basis.cpp
)

set(TEST_EXECUTABLE vmachinetest)
set(PROJECT_EXECUTABLE vmachine)


# Add executable file VMaChine
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

if(armadillo)
  target_link_libraries(${PROJECT_EXECUTABLE} ${ARMADILLO_LIBRARIES})
endif()


# Add test executable file
add_executable(${TEST_EXECUTABLE} ${TESTSOURCES})
# target_link_libraries(${TEST_EXECUTABLE} ${GTEST_LIBRARIES} ${MPI_LIBRARIES})
target_link_libraries(${TEST_EXECUTABLE} ${MPI_LIBRARIES})
add_test(${TEST_EXECUTABLE} ${TEST_EXECUTABLE})

if(MPI_COMPILE_FLAGS)
  set_target_properties(${TEST_EXECUTABLE} PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${TEST_EXECUTABLE} PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

# Add library to make source files visible for doxygen
# add_library(SourceFiles ${SOURCES})
# target_include_directories(SourceFiles PUBLIC .)
