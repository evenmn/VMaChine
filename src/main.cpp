#include <mpi.h>
#include "system.h"
#include <iostream>
#include <string>

#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/gaussian.h"
#include "WaveFunctions/padejastrow.h"
#include "WaveFunctions/slaterdeterminant.h"
#include "WaveFunctions/rbmgaussian.h"
#include "WaveFunctions/rbmjastrow.h"
#include "WaveFunctions/rbmjastrow2.h"
#include "WaveFunctions/rbmjastrow3.h"
#include "WaveFunctions/drbmjastrow.h"
#include "WaveFunctions/simplejastrow.h"
#include "WaveFunctions/partlyrestricted.h"
#include "WaveFunctions/hydrogenlike.h"

#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "Hamiltonians/atomicnucleus.h"
#include "Hamiltonians/doublewell.h"

#include "Basis/basis.h"
#include "Basis/none.h"
#include "Basis/hermite.h"
#include "Basis/hydrogenorbital.h"
#include "Basis/hartreefockhermite.h"

#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "InitialStates/randomnormal.h"

#include "InitialWeights/initialweights.h"
#include "InitialWeights/constant.h"
#include "InitialWeights/randomize.h"
#include "InitialWeights/automatize.h"

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

//#include "Plotter/plotter.h"

int main(int argc, char *argv[]) {
    // MPI initializations
    int numberOfProcesses, myRank;
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numberOfProcesses);
    MPI_Comm_rank (MPI_COMM_WORLD, &myRank);

    // --- SYSTEM SETTINGS ---
    // Parameters
    int     numberOfDimensions  = 2;
    int     numberOfParticles   = 2;
    int     numberOfHiddenNodes = numberOfParticles;
    int     numberOfSteps       = int(pow(2,18));
    int     numberOfIterations  = 1000;
    double  learningRate        = 0.5;
    double  omega               = 1.0;                      // Oscillator frequency
    int     Z                   = numberOfParticles;        // Atomic number (nucleus charge)
    double  sigma               = 1/sqrt(omega);            // Width of probability distribution
    double  stepLength          = 0.1;                      // Metropolis step length
    double  equilibration       = 0.001;                      // Amount of the total steps used

    // Switches
    bool    interaction             = true;                     // Repulsive interaction on or off
    bool    checkConvergence        = false;                    // Stops the program after it has converged
    bool    applyAdaptiveSteps      = true;                     // Increase the number of MC-cycles for the last iterations
    bool    computeOneBodyDensity   = true;                     // Compute one-body density and print to file
    bool    computeTwoBodyDensity   = true;
    bool    printEnergyFile         = true;                     // Print energy for every iteration to file
    bool    printParametersToFile   = true;
    bool    doResampling            = true;                     // Print blocking file for the last iteration and do blocking


    // --- ADVANCED SETTINGS ---
    // Path to data files
    std::string path = "data/";

    // Convergence tools
    int     numberOfEnergies                = 5;            // Check this number of energies for convergence
    double  tolerance                       = 1e-7;         // Convergence tolerance

    // Dynamic step tools
    int     rangeOfDynamicSteps             = 10;           // For how many iterations should we increase # MC-cycles?
    int     additionalSteps                 = 4;            // How much should we increase it? (as a power of 2)
    int     additionalStepsLastIteration    = 8;            // How much should we increase the very last? (as a power of 2)

    // Density tools
    double  maxRadius                       = 10;          // Max radius of one-body density plots
    int     numberOfBins                    = 3000;        // 100 bins per radius unit


    // --- SET PARAMETERS ---
    System* system = new System();

    system->setPath                     (path);
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->setFrequency                (omega);
    system->setAtomicNumber             (Z);
    system->setWidth                    (sigma);
    system->setLearningRate             (learningRate);
    system->setNumberOfParticles        (numberOfParticles);
    system->setNumberOfDimensions       (numberOfDimensions);
    system->setNumberOfHiddenNodes      (numberOfHiddenNodes);
    system->setMPITools                 (myRank, numberOfProcesses);
    system->setNumberOfMetropolisSteps  (numberOfSteps);

    system->setInteraction              (interaction);
    system->setParameterPrintingTools   (printParametersToFile);
    system->setConvergenceTools         (checkConvergence, numberOfEnergies, tolerance);
    system->setDynamicStepTools         (applyAdaptiveSteps, rangeOfDynamicSteps, additionalSteps, additionalStepsLastIteration);
    system->setDensityTools             (computeOneBodyDensity, computeTwoBodyDensity, numberOfBins, maxRadius);
    system->setEnergyPrintingTools      (printEnergyFile, doResampling);

    if(argc == 2) system->parser        (argv[1], numberOfIterations);

    //system->setBasis                    (new HartreeFockHermite(system, new Hermite(system)));
    system->setBasis                    (new Hermite(system));
    std::vector<class WaveFunction*> WaveFunctionElements;
    //WaveFunctionElements.push_back      (new class HydrogenLike      (system));
    WaveFunctionElements.push_back      (new class Gaussian          (system));
    //WaveFunctionElements.push_back      (new class RBMGaussian       (system));
    //WaveFunctionElements.push_back      (new class RBMJastrow        (system));
    //WaveFunctionElements.push_back      (new class SimpleJastrow     (system));
    //WaveFunctionElements.push_back      (new class DRBMJastrow       (system, 3));
    WaveFunctionElements.push_back      (new class SlaterDeterminant (system));
    //WaveFunctionElements.push_back      (new class PartlyRestricted  (system));
    WaveFunctionElements.push_back      (new class PadeJastrow       (system));

    system->setWaveFunctionElements     (WaveFunctionElements);
    system->setNumberOfElements         (WaveFunctionElements.size());
    system->setRandomNumberGenerator    (new MersenneTwister());
    system->setInitialWeights           (new Constant(system, 1.0));
    system->setInitialState             (new RandomNormal(system));
    system->setHamiltonian              (new HarmonicOscillator(system));
    //system->setHamiltonian              (new DoubleWell(system));
    //system->setHamiltonian              (new AtomicNucleus(system));
    system->setGlobalArraysToCalculate  ();
    system->setMetropolis               (new ImportanceSampling(system));
    system->setOptimization             (new ADAM(system));
    system->setGradients                ();
    system->runIterations               (numberOfIterations);

    //class Plotter* plots = new Plotter(system);

    //plots->plotEnergy(argc, argv);
    //plots->plotOneBodyDensity(argc, argv);

    return 0;
}
