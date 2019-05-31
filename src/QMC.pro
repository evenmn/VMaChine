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
SOURCES += main.cpp \
    system.cpp \
    sampler.cpp \
    Hamiltonians/hamiltonian.cpp \
    Hamiltonians/harmonicoscillator.cpp \
    Hamiltonians/atomicnucleus.cpp \
    Hamiltonians/doublewell.cpp \
    InitialStates/initialstate.cpp \
    InitialStates/randomuniform.cpp \
    InitialStates/randomnormal.cpp \
    InitialWeights/initialweights.cpp \
    InitialWeights/randomize.cpp \
    InitialWeights/constant.cpp \
    InitialWeights/automatize.cpp \
    Metropolis/metropolis.cpp \
    Metropolis/bruteforce.cpp \
    Metropolis/importancesampling.cpp \
    WaveFunctions/wavefunction.cpp \
    WaveFunctions/gaussian.cpp \
    WaveFunctions/slaterdeterminant.cpp \
    WaveFunctions/hydrogenlike.cpp \
    WaveFunctions/padejastrow.cpp \
    WaveFunctions/simplejastrow.cpp \
    WaveFunctions/rbmgaussian.cpp \
    WaveFunctions/rbmjastrow.cpp \
    WaveFunctions/drbmjastrow.cpp \
    WaveFunctions/partlyrestricted.cpp \
    Optimization/optimization.cpp \
    Optimization/barzilaiborwein.cpp \
    Optimization/gradientdescent.cpp \
    Optimization/sgd.cpp \
    Optimization/asgd.cpp \
    Optimization/adam.cpp \
    RNG/rng.cpp \
    RNG/mersennetwister.cpp \
    Basis/basis.cpp \
    Basis/hermite.cpp \
    Basis/hydrogenorbital.cpp \
    Basis/none.cpp \
    Basis/hartreefock.cpp \
    Basis/hermiteexpansion.cpp

# Load external source files
SOURCES += \
    block/c++/blocker.cpp

# Load internal header files
HEADERS += \
    system.h \
    sampler.h \
    Hamiltonians/hamiltonian.h \
    Hamiltonians/harmonicoscillator.h \
    Hamiltonians/atomicnucleus.h \
    Hamiltonians/doublewell.h \
    InitialStates/initialstate.h \
    InitialStates/randomuniform.h \
    InitialStates/randomnormal.h \
    InitialWeights/initialweights.h \
    InitialWeights/randomize.h \
    InitialWeights/constant.h \
    InitialWeights/automatize.h \
    Metropolis/metropolis.h \
    Metropolis/bruteforce.h \
    Metropolis/importancesampling.h \
    WaveFunctions/wavefunction.h \
    WaveFunctions/gaussian.h \
    WaveFunctions/slaterdeterminant.h \
    WaveFunctions/hydrogenlike.h \
    WaveFunctions/padejastrow.h \
    WaveFunctions/simplejastrow.h \
    WaveFunctions/rbmgaussian.h \
    WaveFunctions/rbmjastrow.h \
    WaveFunctions/drbmjastrow.h \
    WaveFunctions/partlyrestricted.h \
    Optimization/optimization.h \
    Optimization/barzilaiborwein.h \
    Optimization/gradientdescent.h \
    Optimization/sgd.h \
    Optimization/asgd.h \
    Optimization/adam.h \
    RNG/rng.h \
    RNG/mersennetwister.h \
    Basis/basis.h \
    Basis/hermite.h \
    Basis/hydrogenorbital.h \
    Basis/none.h \
    Basis/hartreefock.h \
    Basis/hermiteexpansion.h

# Load external header files
HEADERS += \
    block/c++/blocker.h \
    tqdm/tqdm.h
