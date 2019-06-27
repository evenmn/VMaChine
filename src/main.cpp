#include <mpi.h>
#include <string>
#include "allheaders.h"

int main(int argc, char *argv[]) {
    // MPI initializations
    int numberOfProcesses, myRank;
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numberOfProcesses);
    MPI_Comm_rank (MPI_COMM_WORLD, &myRank);

    // --- SYSTEM SETTINGS ---
    // Parameters
    int     numberOfDimensions  = 2;
    int     numberOfParticles   = 30;
    int     numberOfHiddenNodes = numberOfParticles;
    int     numberOfSteps       = int(pow(2,15));
    int     numberOfIterations  = 10000;
    double  totalSpin           = 0;                    // totalSpin is half-integer
    double  learningRate        = 0.001;
    double  omega               = 1.0;                 // Oscillator frequency
    int     Z                   = numberOfParticles;    // Atomic number (nucleus charge)
    double  sigma               = 1/sqrt(omega);        // Width of probability distribution
    double  stepLength          = 0.1;                  // Metropolis step length
    double  equilibration       = 0.00;                // Amount of the total steps used

    // Switches
    bool    interaction             = true;     // Repulsive interaction on or off
    bool    checkConvergence        = false;    // Stops the program after it has converged
    bool    applyAdaptiveSteps      = true;     // Increase the number of MC-cycles for the last iterations
    bool    computeOneBodyDensity   = true;     // Compute one-body density and print to file
    bool    computeTwoBodyDensity   = true;
    bool    printEnergyFile         = true;     // Print energy for every iteration to file
    bool    printParametersToFile   = true;
    bool    doResampling            = true;     // Print blocking file for the last iteration and do blocking


    // --- ADVANCED SETTINGS ---
    // Path to data files
    std::string path = "data/";

    // Convergence tools
    int     numberOfEnergies        = 5;        // Check this number of energies for convergence
    double  tolerance               = 1e-7;     // Convergence tolerance

    // Dynamic step toolssimply
    int     rangeOfAdaptiveSteps    = 10;       // For how many iterations should we increase # MC-cycles?
    int     additionalSteps         = 4;        // How much should we increase it? (as a power of 2)
    int     additionalStepsLastIter = 8;        // How much should we increase the very last? (as a power of 2)

    // Density tools
    double  maxRadius               = 15;       // Max radius of one-body density plots
    int     numberOfBins            = 3000;     // 100 bins per radius unit

    // Screening tools
    double  screeningStrength       = 1;        // Screening parameter
    double  dsl                     = 100;      // Debye Screening length


    // --- SET PARAMETERS ---
    System* system = new System();

    system->setPath                     (path);
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->setFrequency                (omega);
    system->setAtomicNumber             (Z);
    system->setWidth                    (sigma);
    system->setTotalSpin                (totalSpin);
    system->setLearningRate             (learningRate);
    system->setNumberOfParticles        (numberOfParticles);
    system->setNumberOfDimensions       (numberOfDimensions);
    system->setNumberOfHiddenNodes      (numberOfHiddenNodes);
    system->setMPITools                 (myRank, numberOfProcesses);
    system->setNumberOfMetropolisSteps  (numberOfSteps);

    system->setInteraction              (interaction);
    system->setParameterPrintingTools   (printParametersToFile);
    system->setConvergenceTools         (checkConvergence, numberOfEnergies, tolerance);
    system->setAdaptiveStepTools        (applyAdaptiveSteps, rangeOfAdaptiveSteps, additionalSteps, additionalStepsLastIter);
    system->setDensityTools             (computeOneBodyDensity, computeTwoBodyDensity, numberOfBins, maxRadius);
    system->setEnergyPrintingTools      (printEnergyFile, doResampling);
    system->setScreeningTools           (screeningStrength, dsl);

    if(argc == 2) system->parserConstants(argv[1], numberOfIterations);

    system->setBasis                    (new Hermite(system));
    std::vector<class WaveFunction*> waveFunctionElements;
    //waveFunctionElements.push_back      (new class Gaussian          (system));
    waveFunctionElements.push_back      (new class SlaterDeterminant (system));
    waveFunctionElements.push_back      (new class RBMGaussian       (system));
    waveFunctionElements.push_back      (new class RBMJastrow        (system));
    //waveFunctionElements.push_back      (new class SimpleJastrow     (system));
    waveFunctionElements.push_back      (new class PadeJastrow       (system));
    //waveFunctionElements.push_back      (new class PartlyRestricted  (system));
    //waveFunctionElements.push_back      (new HydrogenLike      (system));

    system->setWaveFunctionElements     (waveFunctionElements);
    system->setRandomNumberGenerator    (new MersenneTwister());
    system->setOptimization             (new ADAM(system));
    system->setInitialWeights           (new Randomize(system, 0.1));
    system->setInitialState             (new RandomNormal(system));
    system->setHamiltonian              (new HarmonicOscillator(system));
    system->setMetropolis               (new ImportanceSampling(system));

    if(argc == 2) system->parserObjects (argv[1]);

    system->runIterations               (numberOfIterations);

    return 0;
}
