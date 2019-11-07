TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -llapack -lblas -larmadillo

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

# Load arma::uwordernal source files
SOURCES += src-arma/main.cpp \
    src-arma/WaveFunctions/rbmproduct.cpp \
    src-arma/system.cpp \
    src-arma/sampler.cpp \
    src-arma/Hamiltonians/hamiltonian.cpp \
    src-arma/Hamiltonians/harmonicoscillator.cpp \
    src-arma/Hamiltonians/atomicnucleus.cpp \
    src-arma/Hamiltonians/doublewell.cpp \
    src-arma/InitialStates/initialstate.cpp \
    src-arma/InitialStates/randomuniform.cpp \
    src-arma/InitialStates/randomnormal.cpp \
    src-arma/InitialWeights/initialweights.cpp \
    src-arma/InitialWeights/randomize.cpp \
    src-arma/InitialWeights/constant.cpp \
    src-arma/InitialWeights/automatize.cpp \
    src-arma/Metropolis/metropolis.cpp \
    src-arma/Metropolis/bruteforce.cpp \
    src-arma/Metropolis/importancesampling.cpp \
    src-arma/WaveFunctions/wavefunction.cpp \
    src-arma/WaveFunctions/gaussian.cpp \
    src-arma/WaveFunctions/slaterdeterminant.cpp \
    src-arma/WaveFunctions/hydrogenlike.cpp \
    src-arma/WaveFunctions/padejastrow.cpp \
    src-arma/WaveFunctions/simplejastrow.cpp \
    src-arma/WaveFunctions/rbmgaussian.cpp \
    src-arma/WaveFunctions/partlyrestricted.cpp \
    src-arma/Optimization/optimization.cpp \
    src-arma/Optimization/gradientdescent.cpp \
    src-arma/Optimization/sgd.cpp \
    src-arma/Optimization/asgd.cpp \
    src-arma/Optimization/adam.cpp \
    src-arma/RNG/rng.cpp \
    src-arma/RNG/mersennetwister.cpp \
    src-arma/Basis/basis.cpp \
    src-arma/Basis/hermite.cpp \
    src-arma/Basis/hydrogenorbital.cpp \
    src-arma/Basis/none.cpp \
    src-arma/Basis/hartreefock.cpp \
    src-arma/Basis/hermiteexpansion.cpp

# Load external source files
SOURCES += \
    src-arma/block/c++/blocker.cpp

# Load arma::uwordernal header files
HEADERS += \
    src-arma/WaveFunctions/rbmproduct.h \
    src-arma/main.h \
    src-arma/system.h \
    src-arma/sampler.h \
    src-arma/Hamiltonians/hamiltonian.h \
    src-arma/Hamiltonians/harmonicoscillator.h \
    src-arma/Hamiltonians/atomicnucleus.h \
    src-arma/Hamiltonians/doublewell.h \
    src-arma/InitialStates/initialstate.h \
    src-arma/InitialStates/randomuniform.h \
    src-arma/InitialStates/randomnormal.h \
    src-arma/InitialWeights/initialweights.h \
    src-arma/InitialWeights/randomize.h \
    src-arma/InitialWeights/constant.h \
    src-arma/InitialWeights/automatize.h \
    src-arma/Metropolis/metropolis.h \
    src-arma/Metropolis/bruteforce.h \
    src-arma/Metropolis/importancesampling.h \
    src-arma/WaveFunctions/wavefunction.h \
    src-arma/WaveFunctions/gaussian.h \
    src-arma/WaveFunctions/slaterdeterminant.h \
    src-arma/WaveFunctions/hydrogenlike.h \
    src-arma/WaveFunctions/padejastrow.h \
    src-arma/WaveFunctions/simplejastrow.h \
    src-arma/WaveFunctions/rbmgaussian.h \
    src-arma/WaveFunctions/partlyrestricted.h \
    src-arma/Optimization/optimization.h \
    src-arma/Optimization/gradientdescent.h \
    src-arma/Optimization/sgd.h \
    src-arma/Optimization/asgd.h \
    src-arma/Optimization/adam.h \
    src-arma/RNG/rng.h \
    src-arma/RNG/mersennetwister.h \
    src-arma/Basis/basis.h \
    src-arma/Basis/hermite.h \
    src-arma/Basis/hydrogenorbital.h \
    src-arma/Basis/none.h \
    src-arma/Basis/hartreefock.h \
    src-arma/Basis/hermiteexpansion.h

# Load external header files
HEADERS += \
    src-arma/block/c++/blocker.h
