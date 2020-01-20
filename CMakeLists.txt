cmake_minimum_required(VERSION 3.10)
project(VMaChine)

find_package(MPI REQUIRED)
find_package(GTest REQUIRED)

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

# Set Eigen as the default linear algebra library
OPTION(eigen "use Eigen library" ON)
OPTION(armadillo "use Armadillo library" OFF)

# Print which build configuration used
message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

# Add documentation
add_subdirectory ("docs")

if(armadillo)
    find_package(Armadillo REQUIRED)
    # Include directories
    include_directories(
        ${PROJECT_SOURCE_DIR}/src-arma
        ${PROJECT_SOURCE_DIR}/src-arma/Basis
        ${PROJECT_SOURCE_DIR}/src-arma/Hamiltonians
        ${PROJECT_SOURCE_DIR}/src-arma/InitialStates
        ${PROJECT_SOURCE_DIR}/src-arma/InitialWeights
        ${PROJECT_SOURCE_DIR}/src-arma/Metropolis
        ${PROJECT_SOURCE_DIR}/src-arma/Optimization
        ${PROJECT_SOURCE_DIR}/src-arma/RNG
        ${PROJECT_SOURCE_DIR}/src-arma/WaveFunctions
        ${PROJECT_SOURCE_DIR}/src-arma/Eigen
        ${PROJECT_SOURCE_DIR}/src-arma/blocker/c++
        ${PROJECT_SOURCE_DIR}/src-arma/tqdm
        ${PROJECT_SOURCE_DIR}/test
        ${GTEST_INCLUDE_DIRS}
        ${MPI_INCLUDE_PATH}
        ${ARMADILLO_INCLUDE_DIRS}
    )

    # Specify which functions to compile
    set(SOURCES
        src-arma/main.cpp
        src-arma/system.cpp
        src-arma/Hamiltonians/hamiltonian.cpp
        src-arma/Hamiltonians/harmonicoscillator.cpp
        src-arma/Hamiltonians/doublewell.cpp
        src-arma/WaveFunctions/wavefunction.cpp
        src-arma/InitialStates/initialstate.cpp
        src-arma/InitialStates/randomuniform.cpp
        src-arma/sampler.cpp
        src-arma/InitialStates/randomnormal.cpp
        src-arma/Hamiltonians/atomicnucleus.cpp
        src-arma/InitialWeights/initialweights.cpp
        src-arma/InitialWeights/randomuniform.cpp
        src-arma/InitialWeights/automatize.cpp
        src-arma/Metropolis/metropolis.cpp
        src-arma/Metropolis/bruteforce.cpp
        src-arma/Metropolis/importancesampling.cpp
        src-arma/InitialWeights/constant.cpp
        src-arma/Optimization/optimization.cpp
        src-arma/WaveFunctions/gaussian.cpp
        src-arma/WaveFunctions/partlyrestricted.cpp
        src-arma/WaveFunctions/slaterdeterminant.cpp
        src-arma/RNG/rng.cpp
        src-arma/RNG/mersennetwister.cpp
        src-arma/Optimization/gradientdescent.cpp
        src-arma/Basis/basis.cpp
        src-arma/Basis/hermite.cpp
        src-arma/Basis/hartreefock.cpp
        src-arma/Basis/hydrogenorbital.cpp
        src-arma/Basis/hermiteexpansion.cpp
        src-arma/WaveFunctions/hydrogenlike.cpp
        src-arma/Optimization/sgd.cpp
        src-arma/Optimization/asgd.cpp
        src-arma/Optimization/adam.cpp
        src-arma/Basis/none.cpp
        src-arma/WaveFunctions/padejastrow.cpp
        src-arma/WaveFunctions/simplejastrow.cpp
        src-arma/WaveFunctions/rbmgaussian.cpp
        src-arma/WaveFunctions/rbmproduct.cpp
        src-arma/block/c++/blocker.cpp
    )
else()
    # Include directories
    include_directories(
        ${PROJECT_SOURCE_DIR}/src-eigen
        ${PROJECT_SOURCE_DIR}/src-eigen/Activation
        ${PROJECT_SOURCE_DIR}/src-eigen/Basis
        ${PROJECT_SOURCE_DIR}/src-eigen/Hamiltonians
        ${PROJECT_SOURCE_DIR}/src-eigen/InitialStates
        ${PROJECT_SOURCE_DIR}/src-eigen/InitialWeights
        ${PROJECT_SOURCE_DIR}/src-eigen/Layer
        ${PROJECT_SOURCE_DIR}/src-eigen/Metropolis
        ${PROJECT_SOURCE_DIR}/src-eigen/Optimization
        ${PROJECT_SOURCE_DIR}/src-eigen/RNG
        ${PROJECT_SOURCE_DIR}/src-eigen/WaveFunctions
        ${PROJECT_SOURCE_DIR}/src-eigen/Eigen
        ${PROJECT_SOURCE_DIR}/src-eigen/blocker/c++
        ${PROJECT_SOURCE_DIR}/test
        ${GTEST_INCLUDE_DIRS}
        ${MPI_INCLUDE_PATH}
    )

    # Specify which functions to compile
    set(SOURCES
        src-eigen/main.cpp
        src-eigen/system.cpp
        src-eigen/Activation/activation.cpp
        src-eigen/Activation/elu.cpp
        src-eigen/Activation/leakyrelu.cpp
        src-eigen/Activation/relu.cpp
        src-eigen/Activation/purelinear.cpp
        src-eigen/Activation/sigmoid.cpp
        src-eigen/Hamiltonians/hamiltonian.cpp
        src-eigen/Hamiltonians/harmonicoscillator.cpp
        src-eigen/Hamiltonians/doublewell.cpp
        src-eigen/WaveFunctions/wavefunction.cpp
        src-eigen/InitialStates/initialstate.cpp
        src-eigen/InitialStates/randomuniform.cpp
        src-eigen/sampler.cpp
        src-eigen/InitialStates/randomnormal.cpp
        src-eigen/Hamiltonians/atomicnucleus.cpp
        src-eigen/InitialWeights/initialweights.cpp
        src-eigen/InitialWeights/randomuniform.cpp
        src-eigen/InitialWeights/randomnormal.cpp
        src-eigen/InitialWeights/automatize.cpp
        src-eigen/Layer/layer.cpp
        src-eigen/Layer/dense.cpp
        src-eigen/Layer/input.cpp
        src-eigen/Layer/output.cpp
        src-eigen/Metropolis/metropolis.cpp
        src-eigen/Metropolis/bruteforce.cpp
        src-eigen/Metropolis/importancesampling.cpp
        src-eigen/InitialWeights/constant.cpp
        src-eigen/Optimization/optimization.cpp
        src-eigen/WaveFunctions/gaussian.cpp
        src-eigen/WaveFunctions/partlyrestricted.cpp
        src-eigen/WaveFunctions/slaterdeterminant.cpp
        src-eigen/RNG/rng.cpp
        src-eigen/RNG/mersennetwister.cpp
        src-eigen/Optimization/gradientdescent.cpp
        src-eigen/Basis/basis.cpp
        src-eigen/Basis/hermite.cpp
        src-eigen/Basis/hartreefock.cpp
        src-eigen/Basis/hydrogenorbital.cpp
        src-eigen/Basis/hermiteexpansion.cpp
        src-eigen/WaveFunctions/hydrogenlike.cpp
        src-eigen/Optimization/sgd.cpp
        src-eigen/Optimization/asgd.cpp
        src-eigen/Optimization/adam.cpp
        src-eigen/Basis/none.cpp
        src-eigen/WaveFunctions/padejastrow.cpp
        src-eigen/WaveFunctions/simplejastrow.cpp
        src-eigen/WaveFunctions/rbmgaussian.cpp
        src-eigen/WaveFunctions/rbmproduct.cpp
        src-eigen/WaveFunctions/fnn.cpp
        src-eigen/block/c++/blocker.cpp
    )
endif()

set(TESTSOURCES 
    tests/main.cpp
    tests/test_basis.cpp
)

set(TEST_EXECUTABLE VMaChineTests)
set(PROJECT_EXECUTABLE VMaChine)


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

if(armadillo)
  target_link_libraries(${PROJECT_EXECUTABLE} ${ARMADILLO_LIBRARIES})
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

# Add library to make source files visible for doxygen
# add_library(SourceFiles ${SOURCES})
# target_include_directories(SourceFiles PUBLIC .)
