TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG += qt
LIBS += -llapack -lblas

QT += core
QT += charts
# greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

# target.path = $$/home/evenmn/VMaChine/src-eigen/

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
SOURCES += src-eigen/main.cpp \
    src-eigen/Activation/activation.cpp \
    src-eigen/Activation/elu.cpp \
    src-eigen/Activation/leakyrelu.cpp \
    src-eigen/Activation/purelinear.cpp \
    src-eigen/Activation/relu.cpp \
    src-eigen/Activation/sigmoid.cpp \
    src-eigen/Layer/dense.cpp \
    src-eigen/Layer/input.cpp \
    src-eigen/Layer/layer.cpp \
    src-eigen/Plotter/plotter.cpp \
    src-eigen/WaveFunctions/doubleproduct.cpp \
    src-eigen/WaveFunctions/drbmproduct.cpp \
    src-eigen/WaveFunctions/fnn.cpp \
    src-eigen/WaveFunctions/rbmproduct.cpp \
    src-eigen/system.cpp \
    src-eigen/sampler.cpp \
    src-eigen/Hamiltonians/hamiltonian.cpp \
    src-eigen/Hamiltonians/harmonicoscillator.cpp \
    src-eigen/Hamiltonians/atomicnucleus.cpp \
    src-eigen/Hamiltonians/doublewell.cpp \
    src-eigen/InitialStates/initialstate.cpp \
    src-eigen/InitialStates/randomuniform.cpp \
    src-eigen/InitialStates/randomnormal.cpp \
    src-eigen/InitialWeights/initialweights.cpp \
    src-eigen/InitialWeights/randomize.cpp \
    src-eigen/InitialWeights/constant.cpp \
    src-eigen/InitialWeights/automatize.cpp \
    src-eigen/Metropolis/metropolis.cpp \
    src-eigen/Metropolis/bruteforce.cpp \
    src-eigen/Metropolis/importancesampling.cpp \
    src-eigen/WaveFunctions/wavefunction.cpp \
    src-eigen/WaveFunctions/gaussian.cpp \
    src-eigen/WaveFunctions/slaterdeterminant.cpp \
    src-eigen/WaveFunctions/hydrogenlike.cpp \
    src-eigen/WaveFunctions/padejastrow.cpp \
    src-eigen/WaveFunctions/simplejastrow.cpp \
    src-eigen/WaveFunctions/rbmgaussian.cpp \
    src-eigen/WaveFunctions/partlyrestricted.cpp \
    src-eigen/Optimization/optimization.cpp \
    src-eigen/Optimization/gradientdescent.cpp \
    src-eigen/Optimization/sgd.cpp \
    src-eigen/Optimization/asgd.cpp \
    src-eigen/Optimization/adam.cpp \
    src-eigen/RNG/rng.cpp \
    src-eigen/RNG/mersennetwister.cpp \
    src-eigen/Basis/basis.cpp \
    src-eigen/Basis/hermite.cpp \
    src-eigen/Basis/hydrogenorbital.cpp \
    src-eigen/Basis/none.cpp \
    src-eigen/Basis/hartreefock.cpp \
    src-eigen/Basis/hermiteexpansion.cpp

# Load external source files
SOURCES += \
    src-eigen/block/c++/blocker.cpp

# Load internal header files
HEADERS += \
    src-eigen/Activation/activation.h \
    src-eigen/Activation/elu.h \
    src-eigen/Activation/leakyrelu.h \
    src-eigen/Activation/purelinear.h \
    src-eigen/Activation/relu.h \
    src-eigen/Activation/sigmoid.h \
    src-eigen/Layer/dense.h \
    src-eigen/Layer/input.h \
    src-eigen/Layer/layer.h \
    src-eigen/WaveFunctions/doubleproduct.h \
    src-eigen/WaveFunctions/drbmproduct.h \
    src-eigen/WaveFunctions/fnn.h \
    src-eigen/WaveFunctions/rbmproduct.h \
    src-eigen/main.h \
    src-eigen/system.h \
    src-eigen/sampler.h \
    src-eigen/Hamiltonians/hamiltonian.h \
    src-eigen/Hamiltonians/harmonicoscillator.h \
    src-eigen/Hamiltonians/atomicnucleus.h \
    src-eigen/Hamiltonians/doublewell.h \
    src-eigen/InitialStates/initialstate.h \
    src-eigen/InitialStates/randomuniform.h \
    src-eigen/InitialStates/randomnormal.h \
    src-eigen/InitialWeights/initialweights.h \
    src-eigen/InitialWeights/randomize.h \
    src-eigen/InitialWeights/constant.h \
    src-eigen/InitialWeights/automatize.h \
    src-eigen/Metropolis/metropolis.h \
    src-eigen/Metropolis/bruteforce.h \
    src-eigen/Metropolis/importancesampling.h \
    src-eigen/WaveFunctions/wavefunction.h \
    src-eigen/WaveFunctions/gaussian.h \
    src-eigen/WaveFunctions/slaterdeterminant.h \
    src-eigen/WaveFunctions/hydrogenlike.h \
    src-eigen/WaveFunctions/padejastrow.h \
    src-eigen/WaveFunctions/simplejastrow.h \
    src-eigen/WaveFunctions/rbmgaussian.h \
    src-eigen/WaveFunctions/partlyrestricted.h \
    src-eigen/Optimization/optimization.h \
    src-eigen/Optimization/gradientdescent.h \
    src-eigen/Optimization/sgd.h \
    src-eigen/Optimization/asgd.h \
    src-eigen/Optimization/adam.h \
    src-eigen/RNG/rng.h \
    src-eigen/RNG/mersennetwister.h \
    src-eigen/Basis/basis.h \
    src-eigen/Basis/hermite.h \
    src-eigen/Basis/hydrogenorbital.h \
    src-eigen/Basis/none.h \
    src-eigen/Basis/hartreefock.h \
    src-eigen/Basis/hermiteexpansion.h

# Load external header files
HEADERS += \
    src-eigen/block/c++/blocker.h
