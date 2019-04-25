QT += widgets
QT += core
QT += charts

SOURCES += main.cpp \
    system.cpp \
    Hamiltonians/hamiltonian.cpp \
    Hamiltonians/harmonicoscillator.cpp \
    WaveFunctions/wavefunction.cpp \
    InitialStates/initialstate.cpp \
    InitialStates/randomuniform.cpp \
    sampler.cpp \
    InitialStates/randomnormal.cpp \
    Hamiltonians/atomicnucleus.cpp \
    InitialWeights/initialweights.cpp \
    InitialWeights/randomize.cpp \
    Metropolis/metropolis.cpp \
    Metropolis/bruteforce.cpp \
    Metropolis/importancesampling.cpp \
    WaveFunctions/nqsjastrow.cpp \
    InitialWeights/constant.cpp \
    Optimization/optimization.cpp \
    WaveFunctions/gaussian.cpp \
    WaveFunctions/partlyrestricted.cpp \
    WaveFunctions/slaterdeterminant.cpp \
    Optimization/barzilaiborwein.cpp \
    RNG/rng.cpp \
    RNG/mersennetwister.cpp \
    RNG/parkmiller.cpp \
    Optimization/gradientdescent.cpp \
    Plotter/plotter.cpp \
    Basis/basis.cpp \
    Basis/hermite.cpp \
    Basis/hydrogenorbital.cpp \
    WaveFunctions/hydrogenlike.cpp \
    Optimization/sgd.cpp \
    Optimization/asgd.cpp \
    Optimization/adam.cpp \
    Basis/none.cpp \
    WaveFunctions/padejastrow.cpp \
    WaveFunctions/nqsgaussian.cpp \
    Resampling/AutoBlocking/blocker.cpp \
    WaveFunctions/simplejastrow.cpp \
    WaveFunctions/nqsjastrow2.cpp \
    WaveFunctions/samsethjastrow.cpp \
    WaveFunctions/nqsjastrow3.cpp \
    WaveFunctions/padejastrow2.cpp \
    WaveFunctions/nqsgaussian2.cpp

HEADERS += \
    system.h \
    Hamiltonians/hamiltonian.h \
    Hamiltonians/harmonicoscillator.h \
    WaveFunctions/wavefunction.h \
    InitialStates/initialstate.h \
    InitialStates/randomuniform.h \
    sampler.h \
    InitialStates/randomnormal.h \
    Hamiltonians/atomicnucleus.h \
    InitialWeights/initialweights.h \
    InitialWeights/randomize.h \
    Metropolis/metropolis.h \
    Metropolis/bruteforce.h \
    Metropolis/importancesampling.h \
    WaveFunctions/nqsjastrow.h \
    InitialWeights/constant.h \
    Optimization/optimization.h \
    WaveFunctions/gaussian.h \
    WaveFunctions/partlyrestricted.h \
    WaveFunctions/slaterdeterminant.h \
    Optimization/barzilaiborwein.h \
    RNG/rng.h \
    RNG/mersennetwister.h \
    RNG/parkmiller.h \
    Optimization/gradientdescent.h \
    Plotter/plotter.h \
    Basis/basis.h \
    Basis/hermite.h \
    Basis/hydrogenorbital.h \
    WaveFunctions/hydrogenlike.h \
    Optimization/sgd.h \
    Optimization/asgd.h \
    Optimization/adam.h \
    Basis/none.h \
    WaveFunctions/padejastrow.h \
    WaveFunctions/nqsgaussian.h \
    Resampling/AutoBlocking/blocker.h \
    WaveFunctions/simplejastrow.h \
    WaveFunctions/nqsjastrow2.h \
    WaveFunctions/samsethjastrow.h \
    WaveFunctions/nqsjastrow3.h \
    WaveFunctions/padejastrow2.h \
    WaveFunctions/nqsgaussian2.h

