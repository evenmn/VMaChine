QT += widgets
QT += core
QT += charts

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK

SOURCES += main.cpp \
    system.cpp \
    sampler.cpp \
    Hamiltonians/hamiltonian.cpp \
    Hamiltonians/harmonicoscillator.cpp \
    Hamiltonians/atomicnucleus.cpp \
    InitialStates/initialstate.cpp \
    InitialStates/randomuniform.cpp \
    InitialStates/randomnormal.cpp \
    InitialWeights/initialweights.cpp \
    InitialWeights/randomize.cpp \
    InitialWeights/constant.cpp \
    Metropolis/metropolis.cpp \
    Metropolis/bruteforce.cpp \
    Metropolis/importancesampling.cpp \
    WaveFunctions/wavefunction.cpp \
    WaveFunctions/gaussian.cpp \
    WaveFunctions/partlyrestricted.cpp \
    WaveFunctions/slaterdeterminant.cpp \
    WaveFunctions/hydrogenlike.cpp \
    WaveFunctions/padejastrow.cpp \
    WaveFunctions/padejastrow2.cpp \
    WaveFunctions/simplejastrow.cpp \
    WaveFunctions/samsethjastrow.cpp \
    WaveFunctions/rbmgaussian.cpp \
    WaveFunctions/rbmgaussian2.cpp \
    WaveFunctions/rbmjastrow.cpp \
    WaveFunctions/rbmjastrow2.cpp \
    WaveFunctions/rbmjastrow3.cpp \
    Optimization/optimization.cpp \
    Optimization/barzilaiborwein.cpp \
    Optimization/gradientdescent.cpp \
    Optimization/sgd.cpp \
    Optimization/asgd.cpp \
    Optimization/adam.cpp \
    RNG/rng.cpp \
    RNG/mersennetwister.cpp \
    RNG/parkmiller.cpp \
    Basis/basis.cpp \
    Basis/hermite.cpp \
    Basis/hydrogenorbital.cpp \
    Basis/none.cpp \
    Resampling/blocker.cpp

HEADERS += \
    system.h \
    sampler.h \
    Hamiltonians/hamiltonian.h \
    Hamiltonians/harmonicoscillator.h \
    Hamiltonians/atomicnucleus.h \
    InitialStates/initialstate.h \
    InitialStates/randomuniform.h \
    InitialStates/randomnormal.h \
    InitialWeights/initialweights.h \
    InitialWeights/randomize.h \
    InitialWeights/constant.h \
    Metropolis/metropolis.h \
    Metropolis/bruteforce.h \
    Metropolis/importancesampling.h \
    WaveFunctions/wavefunction.h \
    WaveFunctions/gaussian.h \
    WaveFunctions/partlyrestricted.h \
    WaveFunctions/slaterdeterminant.h \
    WaveFunctions/hydrogenlike.h \
    WaveFunctions/padejastrow.h \
    WaveFunctions/padejastrow2.h \
    WaveFunctions/simplejastrow.h \
    WaveFunctions/samsethjastrow.h \
    WaveFunctions/rbmgaussian.h \
    WaveFunctions/rbmgaussian2.h \
    WaveFunctions/rbmjastrow.h \
    WaveFunctions/rbmjastrow2.h \
    WaveFunctions/rbmjastrow3.h \
    Optimization/optimization.h \
    Optimization/barzilaiborwein.h \
    Optimization/gradientdescent.h \
    Optimization/sgd.h \
    Optimization/asgd.h \
    Optimization/adam.h \
    RNG/rng.h \
    RNG/mersennetwister.h \
    RNG/parkmiller.h \
    Basis/basis.h \
    Basis/hermite.h \
    Basis/hydrogenorbital.h \
    Basis/none.h \
    Resampling/blocker.h
