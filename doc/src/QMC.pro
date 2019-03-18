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
    WaveFunctions/mlgaussian.cpp \
    Metropolis/metropolis.cpp \
    Metropolis/bruteforce.cpp \
    Metropolis/importancesampling.cpp \
    WaveFunctions/nqsjastrow.cpp \
    InitialWeights/constant.cpp \
    Optimization/optimization.cpp \
    WaveFunctions/gaussian.cpp \
    WaveFunctions/padejastrow.cpp \
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
    WaveFunctions/padejastrownew.cpp \
    WaveFunctions/nqsjastrowold.cpp

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
    WaveFunctions/mlgaussian.h \
    Metropolis/metropolis.h \
    Metropolis/bruteforce.h \
    Metropolis/importancesampling.h \
    WaveFunctions/nqsjastrow.h \
    InitialWeights/constant.h \
    Optimization/optimization.h \
    WaveFunctions/gaussian.h \
    WaveFunctions/padejastrow.h \
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
    WaveFunctions/padejastrownew.h \
    WaveFunctions/nqsjastrowold.h

