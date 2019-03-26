#include "system.h"

#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/gaussian.h"
#include "WaveFunctions/padejastrow.h"
#include "WaveFunctions/slaterdeterminant.h"
#include "WaveFunctions/nqsgaussian.h"
#include "WaveFunctions/nqsjastrow.h"
#include "WaveFunctions/simplejastrow.h"
#include "WaveFunctions/partlyrestricted.h"
#include "WaveFunctions/hydrogenlike.h"

#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/atomicnucleus.h"

#include "Basis/basis.h"
#include "Basis/none.h"
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
#include "Optimization/sgd.h"
#include "Optimization/adam.h"

#include "RNG/rng.h"
#include "RNG/mersennetwister.h"
#include "RNG/parkmiller.h"

#include "Plotter/plotter.h"

int main(int argc, char *argv[]) {
    int     numberOfDimensions  = 3;
    int     numberOfParticles   = 8;
    int     numberOfHiddenNodes = numberOfParticles;
    int     numberOfSteps       = int(pow(2,22));
    int     numberOfIterations  = 1000;
    double  eta                 = 0.01;                      // Learning rate
    double  omega               = 0.1;                      // Oscillator frequency
    int     Z                   = numberOfParticles;        // Atomic number (nucleus charge)
    double  sigma               = 1/sqrt(omega);            // Width of probability distribution
    double  stepLength          = 0.1;                      // Metropolis step length
    double  equilibration       = 0.2;                      // Amount of the total steps used
    bool    interaction         = true;
    int     maxNumberOfParametersPerElement = numberOfParticles*numberOfDimensions*numberOfParticles*numberOfDimensions;

    System* system = new System();
    system->setEquilibrationFraction    (equilibration);
    system->setInteraction              (interaction);
    system->setStepLength               (stepLength);
    system->setFrequency                (omega);
    system->setAtomicNumber             (Z);
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
    WaveFunctionElements.push_back      (new class Gaussian             (system));
    //WaveFunctionElements.push_back      (new class NQSGaussian          (system));
    //WaveFunctionElements.push_back      (new class NQSJastrow           (system));
    //WaveFunctionElements.push_back      (new class SimpleJastrow        (system));
    WaveFunctionElements.push_back      (new class SlaterDeterminant    (system));
    WaveFunctionElements.push_back      (new class PadeJastrow          (system));

    system->setNumberOfWaveFunctionElements(int(WaveFunctionElements.size()));
    system->setWaveFunction             (WaveFunctionElements);
    system->setRandomNumberGenerator    (new MersenneTwister());
    system->setInitialWeights           (new Constant(system, 1.0));
    system->setInitialState             (new RandomNormal(system));
    system->setHamiltonian              (new HarmonicOscillator(system));
    system->setMetropolis               (new ImportanceSampling(system));
    system->setOptimization             (new SGD(system,0.0,0.0));
    system->setGradients                ();
    system->runMetropolisSteps          (numberOfIterations);

    //class Plotter* plots = new Plotter(system);

    //plots->plotEnergy(argc, argv);
    //plots->plotOneBodyDensity(argc, argv);

    return 0;
}
