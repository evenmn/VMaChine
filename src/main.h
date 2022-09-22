#include "sampler.h"
#include "system.h"

#include "WaveFunctions/doubleproduct.h"
#include "WaveFunctions/gaussian.h"
#include "WaveFunctions/hydrogenlike.h"
#include "WaveFunctions/padejastrow.h"
#include "WaveFunctions/partlyrestricted.h"
#include "WaveFunctions/rbmgaussian.h"
#include "WaveFunctions/rbmproduct.h"
#include "WaveFunctions/simplejastrow.h"
#include "WaveFunctions/slaterdeterminant.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/fnn.h"
#include "WaveFunctions/hardcorejastrow.h"

#include "Activation/activation.h"
#include "Activation/elu.h"
#include "Activation/leakyrelu.h"
#include "Activation/purelinear.h"
#include "Activation/relu.h"
#include "Activation/sigmoid.h"

#include "Layer/layer.h"
#include "Layer/dense.h"
#include "Layer/input.h"
#include "Layer/output.h"

#include "Hamiltonians/atomicnucleus.h"
#include "Hamiltonians/doublewell.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/ellipticalharmonicoscillator.h"

#include "Basis/basis.h"
//#include "Basis/hartreefock.h"
#include "Basis/hermite.h"
//#include "Basis/hermiteexpansion.h"
//#include "Basis/hermitespin.h"
#include "Basis/hydrogenorbital.h"
#include "Basis/none.h"

#include "InitialStates/initialstate.h"
#include "InitialStates/randomnormal.h"
#include "InitialStates/randomuniform.h"

#include "InitialWeights/automatize.h"
#include "InitialWeights/constant.h"
#include "InitialWeights/initialweights.h"
#include "InitialWeights/randomuniform.h"
#include "InitialWeights/randomnormal.h"
#include "InitialWeights/customized.h"
#include "InitialWeights/xavier.h"
#include "InitialWeights/fromfile.h"

#include "Interaction/interaction.h"
#include "Interaction/nointeraction.h"
#include "Interaction/coulomb.h"

#include "Metropolis/bruteforce.h"
#include "Metropolis/importancesampling.h"
#include "Metropolis/metropolis.h"

#include "Optimization/adam.h"
#include "Optimization/gradientdescent.h"
#include "Optimization/optimization.h"
#include "Optimization/sgd.h"

#include "RNG/mersennetwister.h"
#include "RNG/rng.h"
