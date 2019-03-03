#include "system.h"

#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/gaussian.h"
#include "WaveFunctions/mlgaussian.h"
#include "WaveFunctions/hydrogenlike.h"
#include "WaveFunctions/padejastrow.h"
#include "WaveFunctions/nqsjastrow.h"
#include "WaveFunctions/nqsjastrowreal.h"
#include "WaveFunctions/partlyrestricted.h"
#include "WaveFunctions/slaterdeterminant.h"

#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/atomicnucleus.h"

#include "Basis/basis.h"
#include "Basis/hermite.h"
#include "Basis/hydrogenorbital.h"

#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "InitialStates/randomnormal.h"

#include "InitialWeights/initialweights.h"
#include "InitialWeights/constant.h"
#include "InitialWeights/randomize.h"

#include "Metropolis/metropolis.h"
#include "Metropolis/bruteforce.h"
#include "Metropolis/importancesampling.h"

#include "Optimization/optimization.h"
#include "Optimization/gradientdescent.h"
#include "Optimization/barzilaiborwein.h"

#include "RNG/rng.h"
#include "RNG/mersennetwister.h"
#include "RNG/parkmiller.h"

#include "Plotter/plotter.cpp"

int main(int argc, char *argv[]) {
    int     numberOfDimensions  = 2;
    int     numberOfParticles   = 2;
    int     numberOfHiddenNodes = 2;
    int     numberOfSteps       = int(pow(2,20));
    int     numberOfIterations  = 100;
    double  eta                 = 0.05;         // Learning rate
    double  omega               = 1.0;          // Oscillator frequency
    double  sigma               = 1.0;          // Width of probability distribution
    double  stepLength          = 0.1;          // Metropolis step length
    double  equilibration       = 0.0;          // Amount of the total steps used
    bool    interaction         = true;
    int     maxNumberOfParametersPerElement = numberOfParticles*numberOfDimensions*numberOfParticles*numberOfDimensions;

    System* system = new System();
    system->setEquilibrationFraction    (equilibration);
    system->setInteraction              (interaction);
    system->setStepLength               (stepLength);
    system->setFrequency                (omega);
    system->setWidth                    (sigma);
    system->setLearningRate             (eta);
    system->setNumberOfParticles        (numberOfParticles);
    system->setNumberOfDimensions       (numberOfDimensions);
    system->setNumberOfHiddenNodes      (numberOfHiddenNodes);
    system->setNumberOfMetropolisSteps  (numberOfSteps);
    system->setMaxNumberOfParametersPerElement (maxNumberOfParametersPerElement);
    system->setTotalNumberOfSteps       ();
    system->setNumberOfFreeDimensions   ();

    system->setBasis                    (new Hermite(system));
    std::vector<class WaveFunction*> WaveFunctionElements;
    //WaveFunctionElements.push_back      (new class HydrogenLike         (system));
    WaveFunctionElements.push_back      (new class Gaussian             (system));
    //WaveFunctionElements.push_back      (new class MLGaussian           (system));
    //WaveFunctionElements.push_back      (new class NQSJastrow           (system));
    //WaveFunctionElements.push_back      (new class PartlyRestricted     (system));
    //WaveFunctionElements.push_back      (new class SlaterDeterminant    (system));
    WaveFunctionElements.push_back      (new class PadeJastrow          (system));

    system->setNumberOfWaveFunctionElements(int(WaveFunctionElements.size()));
    system->setWaveFunction             (WaveFunctionElements);
    system->setRandomNumberGenerator    (new MersenneTwister());
    system->setInitialWeights           (new Randomize(system, 0.1));
    system->setInitialState             (new RandomNormal(system));
    system->setHamiltonian              (new HarmonicOscillator(system));
    system->setMetropolis               (new ImportanceSampling(system));
    system->setOptimization             (new GradientDescent(system, 0.0));
    system->setGradients                ();
    system->runMetropolisSteps          (numberOfIterations);

    //plotEnergy(argc, argv);

    return 0;
}
