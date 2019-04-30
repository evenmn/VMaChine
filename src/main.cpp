#include <mpi.h>
#include "system.h"
#include <iostream>
#include <string>

#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/gaussian.h"
#include "WaveFunctions/padejastrow.h"
#include "WaveFunctions/padejastrow2.h"
#include "WaveFunctions/slaterdeterminant.h"
#include "WaveFunctions/rbmgaussian.h"
#include "WaveFunctions/rbmgaussian2.h"
#include "WaveFunctions/rbmjastrow.h"
#include "WaveFunctions/rbmjastrow2.h"
#include "WaveFunctions/rbmjastrow3.h"
//#include "WaveFunctions/rbmjastrow4.h"
#include "WaveFunctions/rbmjastrow5.h"
#include "WaveFunctions/simplejastrow.h"
#include "WaveFunctions/samsethjastrow.h"
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
    int     numberOfIterations  = 100;
    double  eta                 = 0.5;                      // Learning rate
    double  omega               = 1.0;                      // Oscillator frequency
    int     Z                   = numberOfParticles;        // Atomic number (nucleus charge)
    double  sigma               = 1/sqrt(omega);            // Width of probability distribution
    double  stepLength          = 0.1;                      // Metropolis step length
    double  equilibration       = 0.2;                      // Amount of the total steps used

    // Switches
    bool    interaction         = true;                     // Repulsive interaction on or off
    bool    checkConvergence    = false;                    // Stops the program after it has converged
    bool    applyDynamicSteps   = false;                    // Increase the number of MC-cycles for the last iterations
    bool    computeDensity      = true;                     // Compute one-body density and print to file
    bool    printEnergyFile     = true;                     // Print energy for every iteration to file
    bool    doBlocking          = true;                     // Print blocking file for the last iteration and do blocking


    // --- ADVANCED SETTINGS ---
    // Path to data file
    std::string path = "/home/evenmn/VMC/data/";

    // Convergence tools
    int     numberOfEnergies                = 5;            // Check this number of energies for convergence
    double  tolerance                       = 1e-7;         // Convergence tolerance

    // Dynamic step tools
    int     rangeOfDynamicSteps             = 10;           // For how many iterations should we increase # MC-cycles?
    int     additionalSteps                 = 4;            // How much should we increase it? (as a power of 2)
    int     additionalStepsLastIteration    = 8;            // How much should we increase the very last? (as a power of 2)

    // Density tools
    double  maxRadius                       = 15;                       // Max radius of one-body density plots
    int     numberOfBins                    = int(100 * maxRadius);     // Number of bins


    // --- SET PARAMETERS ---
    System* system = new System();
    system->setPath                     (path);
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->setFrequency                (omega);
    system->setAtomicNumber             (Z);
    system->setWidth                    (sigma);
    system->setLearningRate             (eta);
    system->setNumberOfParticles        (numberOfParticles);
    system->setNumberOfDimensions       (numberOfDimensions);
    system->setNumberOfHiddenNodes      (numberOfHiddenNodes);
    system->setMPITools                 (myRank, numberOfProcesses);
    system->setNumberOfMetropolisSteps  (numberOfSteps);
    system->setTotalNumberOfSteps       ();
    system->setNumberOfFreeDimensions   ();

    system->setInteraction              (interaction);
    system->setConvergenceTools         (checkConvergence, numberOfEnergies, tolerance);
    system->setDynamicStepTools         (applyDynamicSteps, rangeOfDynamicSteps, additionalSteps, additionalStepsLastIteration);
    system->setDensityTools             (computeDensity, numberOfBins, maxRadius);
    system->setEnergyPrintingTools      (printEnergyFile, doBlocking);

    system->setBasis                    (new Hermite(system));
    std::vector<class WaveFunction*> WaveFunctionElements;
    //WaveFunctionElements.push_back      (new class HydrogenLike         (system));
    WaveFunctionElements.push_back      (new class Gaussian             (system));
    //WaveFunctionElements.push_back      (new class RBMGaussian          (system));
    //WaveFunctionElements.push_back      (new class RBMGaussian2         (system));
    //WaveFunctionElements.push_back      (new class RBMJastrow           (system));
    //WaveFunctionElements.push_back      (new class SimpleJastrow        (system));
    //WaveFunctionElements.push_back      (new class RBMJastrow2          (system));
    //WaveFunctionElements.push_back      (new class RBMJastrow5          (system));
    //WaveFunctionElements.push_back      (new class SlaterDeterminant    (system));
    //WaveFunctionElements.push_back      (new class PartlyRestricted     (system));
    WaveFunctionElements.push_back      (new class PadeJastrow          (system));
    //WaveFunctionElements.push_back      (new class PadeJastrow2         (system));
    //WaveFunctionElements.push_back      (new class SamsethJastrow       (system));

    system->setNumberOfWaveFunctionElements(int(WaveFunctionElements.size()));
    system->setWaveFunctionElements     (WaveFunctionElements);
    system->setMaxNumberOfParametersPerElement ();
    system->setRandomNumberGenerator    (new MersenneTwister());
    system->setInitialWeights           (new Constant(system, 1.0));
    system->setInitialState             (new RandomNormal(system));
    system->setHamiltonian              (new HarmonicOscillator(system));
    system->setGlobalArraysToCalculate  ();
    system->setMetropolis               (new ImportanceSampling(system));
    system->setOptimization             (new SGD(system,0.0,0.0));
    system->setGradients                ();
    system->runIterations               (numberOfIterations);

    //class Plotter* plots = new Plotter(system);

    //plots->plotEnergy(argc, argv);
    //plots->plotOneBodyDensity(argc, argv);

    MPI_Finalize();
    return 0;
}
