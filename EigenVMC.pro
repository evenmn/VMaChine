TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -llapack -lblas

QT += widgets
QT += core
QT += charts

# Remove possible other optimization flags
QMAKE_CXXFLAGS_RELEASE -= -O
QMAKE_CXXFLAGS_RELEASE -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2

# Add the desired -O3 if not present
QMAKE_CXXFLAGS_RELEASE *= -O3

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE       = $$QMAKE_CXX
QMAKE_CXX_DEBUG         = $$QMAKE_CXX
QMAKE_LINK              = $$QMAKE_CXX
QMAKE_CC                = mpicc

QMAKE_CFLAGS           += $$system(mpicc  --showme:compile)
QMAKE_LFLAGS           += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS         += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK

# Load internal source files
SOURCES += src/main.cpp \
    src/WaveFunctions/drbmproduct.cpp \
    src/WaveFunctions/rbmproduct.cpp \
    src/system.cpp \
    src/sampler.cpp \
    src/Hamiltonians/hamiltonian.cpp \
    src/Hamiltonians/harmonicoscillator.cpp \
    src/Hamiltonians/atomicnucleus.cpp \
    src/Hamiltonians/doublewell.cpp \
    src/InitialStates/initialstate.cpp \
    src/InitialStates/randomuniform.cpp \
    src/InitialStates/randomnormal.cpp \
    src/InitialWeights/initialweights.cpp \
    src/InitialWeights/randomize.cpp \
    src/InitialWeights/constant.cpp \
    src/InitialWeights/automatize.cpp \
    src/Metropolis/metropolis.cpp \
    src/Metropolis/bruteforce.cpp \
    src/Metropolis/importancesampling.cpp \
    src/WaveFunctions/wavefunction.cpp \
    src/WaveFunctions/gaussian.cpp \
    src/WaveFunctions/slaterdeterminant.cpp \
    src/WaveFunctions/hydrogenlike.cpp \
    src/WaveFunctions/padejastrow.cpp \
    src/WaveFunctions/simplejastrow.cpp \
    src/WaveFunctions/rbmgaussian.cpp \
    src/WaveFunctions/partlyrestricted.cpp \
    src/Optimization/optimization.cpp \
    src/Optimization/gradientdescent.cpp \
    src/Optimization/sgd.cpp \
    src/Optimization/asgd.cpp \
    src/Optimization/adam.cpp \
    src/RNG/rng.cpp \
    src/RNG/mersennetwister.cpp \
    src/Basis/basis.cpp \
    src/Basis/hermite.cpp \
    src/Basis/hydrogenorbital.cpp \
    src/Basis/none.cpp \
    src/Basis/hartreefock.cpp \
    src/Basis/hermiteexpansion.cpp

# Load external source files
SOURCES += \
    src/block/c++/blocker.cpp

# Load internal header files
HEADERS += \
    src/WaveFunctions/drbmproduct.h \
    src/WaveFunctions/rbmproduct.h \
    src/allheaders.h \
    src/system.h \
    src/sampler.h \
    src/Hamiltonians/hamiltonian.h \
    src/Hamiltonians/harmonicoscillator.h \
    src/Hamiltonians/atomicnucleus.h \
    src/Hamiltonians/doublewell.h \
    src/InitialStates/initialstate.h \
    src/InitialStates/randomuniform.h \
    src/InitialStates/randomnormal.h \
    src/InitialWeights/initialweights.h \
    src/InitialWeights/randomize.h \
    src/InitialWeights/constant.h \
    src/InitialWeights/automatize.h \
    src/Metropolis/metropolis.h \
    src/Metropolis/bruteforce.h \
    src/Metropolis/importancesampling.h \
    src/WaveFunctions/wavefunction.h \
    src/WaveFunctions/gaussian.h \
    src/WaveFunctions/slaterdeterminant.h \
    src/WaveFunctions/hydrogenlike.h \
    src/WaveFunctions/padejastrow.h \
    src/WaveFunctions/simplejastrow.h \
    src/WaveFunctions/rbmgaussian.h \
    src/WaveFunctions/partlyrestricted.h \
    src/Optimization/optimization.h \
    src/Optimization/gradientdescent.h \
    src/Optimization/sgd.h \
    src/Optimization/asgd.h \
    src/Optimization/adam.h \
    src/RNG/rng.h \
    src/RNG/mersennetwister.h \
    src/Basis/basis.h \
    src/Basis/hermite.h \
    src/Basis/hydrogenorbital.h \
    src/Basis/none.h \
    src/Basis/hartreefock.h \
    src/Basis/hermiteexpansion.h

# Load external header files
HEADERS += \
    src/block/c++/blocker.h
