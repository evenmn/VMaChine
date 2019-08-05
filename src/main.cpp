#include "allheaders.h"
#include <mpi.h>
#include <string>

int main(int argc, char *argv[])
{
    // MPI initializations
    int numberOfProcesses, myRank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    // --- SYSTEM SETTINGS ---
    // Parameters
    int numberOfDimensions = 2;
    int numberOfParticles = 56;
    int numberOfHiddenNodes = numberOfParticles;
    int numberOfSteps = int(pow(2, 18));
    int numberOfIterations = 3000;
    double totalSpin = 0;           // TotalSpin is half-integer
    double learningRate = 0.0001;   // Learning rate
    double omega = 0.1;             // Oscillator frequency
    int Z = numberOfParticles;      // Atomic number (nucleus charge)
    double sigma = 1 / sqrt(omega); // Width of probability distribution
    double stepLength = 0.1;        // Metropolis step length
    double equilibration = 0.0;     // Amount of the total steps used

    // Switches
    bool interaction = true;           // Repulsive interaction on or off
    bool checkConvergence = false;     // Stops the program after it has converged
    bool applyAdaptiveSteps = true;    // Increase the number of MC-cycles for the last iterations
    bool computeOneBodyDensity = true; // Compute one-body density and print to file
    bool computeTwoBodyDensity = true; // Compute one-body density and print to file
    bool printEnergyFile = true;       // Print energy for every iteration to file
    bool printParametersToFile = true; // Print parameter matrix to file
    bool doResampling = true;          // Print blocking file for the last iteration and do blocking
    bool screening = false;

    // --- ADVANCED SETTINGS ---
    // Path to data files
    std::string path = "data/";

    // Convergence tools
    int numberOfEnergies = 5; // Check this number of energies for convergence
    double tolerance = 1e-7;  // Convergence tolerance

    // Dynamic step tools
    int rangeOfAdaptiveSteps = 10;   // For how many iterations should we increase # MC-cycles?
    int additionalSteps = 4;         // How much should we increase it? (as a power of 2)
    int additionalStepsLastIter = 8; // How much should we increase the very last? (as a power of 2)

    // Density tools
    double maxRadius = 75;   // Max radius of one-body density plots
    int numberOfBins = 3000; // 100 bins per radius unit

    // Screening tools
    double screeningStrength = 1; // Screening parameter
    double dsl = 100;             // Debye Screening length

    // --- SET PARAMETERS ---
    System *system = new System();

    system->setPath(path);
    system->setEquilibrationFraction(equilibration);
    system->setStepLength(stepLength);
    system->setFrequency(omega);
    system->setAtomicNumber(Z);
    system->setWidth(sigma);
    system->setTotalSpin(totalSpin);
    system->setLearningRate(learningRate);
    system->setNumberOfParticles(numberOfParticles);
    system->setNumberOfDimensions(numberOfDimensions);
    system->setNumberOfHiddenNodes(numberOfHiddenNodes);
    system->setMPITools(myRank, numberOfProcesses);
    system->setNumberOfMetropolisSteps(numberOfSteps);

    system->setInteraction(interaction);
    system->setParameterPrintingTools(printParametersToFile);
    system->setConvergenceTools(checkConvergence, numberOfEnergies, tolerance);
    system->setAdaptiveStepTools(applyAdaptiveSteps,
                                 rangeOfAdaptiveSteps,
                                 additionalSteps,
                                 additionalStepsLastIter);
    system->setDensityTools(computeOneBodyDensity, computeTwoBodyDensity, numberOfBins, maxRadius);
    system->setEnergyPrintingTools(printEnergyFile, doResampling);
    system->setScreeningTools(screening, screeningStrength, dsl);

    if (argc == 2)
        system->parserConstants(argv[1], numberOfIterations);

    system->setBasis(new Hermite(system));
    //system->setWaveFunctionElement(new Gaussian(system));
    system->setWaveFunctionElement(new SlaterDeterminant(system));
    system->setWaveFunctionElement(new RBMGaussian(system));
    system->setWaveFunctionElement(new RBMProduct(system));
    //system->setWaveFunctionElement(new PadeJastrow(system));

    system->setRandomNumberGenerator(new MersenneTwister());
    system->setOptimization(new ADAM(system));
    system->setInitialWeights(new Automatize(system));
    system->setInitialState(new RandomNormal(system));
    system->setHamiltonian(new HarmonicOscillator(system));
    system->setMetropolis(new ImportanceSampling(system));

    if (argc == 2)
        system->parserObjects(argv[1]);

    system->runIterations(numberOfIterations);

    return 0;
}
