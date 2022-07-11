#include <iostream>
#include <fstream>
#include <cassert>
#include <ctime>
#include <chrono>
#include <iterator>
#include <mpi.h>
#include <string>
#include <vector>
#include <Eigen/Dense>

#include "system.h"

/* ----------------------------------------------------------------------------
  System class constructor. Setting default classes
---------------------------------------------------------------------------- */

System::System()
{
    m_basis = new None(this);
    m_initialState = new RandomNormal(this);
    m_initialWeights = new Automatize(this);
    m_metropolis = new ImportanceSampling(this);
    m_optimization = new ADAM(this);
    m_randomNumberGenerator = new MersenneTwister();
    m_interactionStyle = new NoInteraction(this);
}


/* ----------------------------------------------------------------------------
  Print logo to terminal
---------------------------------------------------------------------------- */

void System::printLogo()
{
    std::cout << "__      ____  __          _____ _    _ _____ _   _ ______ " << std::endl;
    std::cout << "\\ \\    / /  \\/  |   /\\   / ____| |  | |_   _| \\ | |  ____|" << std::endl;
    std::cout << " \\ \\  / /| \\  / |  /  \\ | |    | |__| | | | |  \\| | |__   " << std::endl;
    std::cout << "  \\ \\/ / | |\\/| | / /\\ \\| |    |  __  | | | | . ` |  __|  " << std::endl;
    std::cout << "   \\  /  | |  | |/ ____ \\ |____| |  | |_| |_| |\\  | |____ " << std::endl;
    std::cout << "    \\/   |_|  |_/_/    \\_\\_____|_|  |_|_____|_| \\_|______|" << std::endl;
    std::cout << "==========================================================" << std::endl;
    std::cout << std::endl;
}


/* ----------------------------------------------------------------------------
  Initialize MPI based on command line arguments. Command line arguments are
  taken automatically from main
---------------------------------------------------------------------------- */

void System::initializeMPI()
{
    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized)
        MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &m_numberOfProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
}


/* ----------------------------------------------------------------------------
  Initialize system. The various objects need to be initialized in the correct
  order, which is the task of this function. This allows for the configuration
  keys to be in arbitrary order in main/configuration file
---------------------------------------------------------------------------- */

void System::initializeSystem()
{
    initializeMPI();
    parser(m_configFile);
    setInputLayer(m_degreesOfFreedom);      // Add input layer
    //setOutputLayer(new Sigmoid(this));      // Add output layer
    m_hamiltonian->initialize();
    m_basis->initialize();
    setAllConstants();
    setMaxParameters();
    setGradients();
    m_optimization->initialize();
    m_initialWeights->setupInitialWeights();
    m_parameters = m_initialWeights->getParameters();
    m_initialState->setupInitialState();
    m_positions = m_initialState->getParticles();
    m_distanceMatrix = m_initialState->getDistanceMatrix();
    m_radialVector = m_initialState->getRadialVector();
    m_metropolis->initialize();
    m_interactionStyle->initialize();

    m_sampler = new Sampler(this); //sampler has to be initialized lastly 
    m_sampler->openOutputFiles();
}


/* ----------------------------------------------------------------------------
  Print initial information to terminal
---------------------------------------------------------------------------- */

void System::printInitialInformation()
{
    m_start = std::chrono::system_clock::now();
    std::time_t start_time = std::chrono::system_clock::to_time_t(m_start);
    std::cout << "Started computation at " << std::ctime(&start_time)
              << "Running on " << m_numberOfProcesses << " CPU threads using Open MPI" << std::endl;
    std::cout << "Simulation is run from the directory: " << m_path << std::endl;
}


/* ----------------------------------------------------------------------------
  Print system information to terminal
---------------------------------------------------------------------------- */

void System::printSystemInformation()
{
    std::cout << std::fixed;
    std::cout << std::boolalpha;
    std::cout << std::setprecision(6);
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "             SYSTEM INFORMATION" << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Number of particles:      " << m_numberOfParticles << std::endl;
    std::cout << "Number of dimensions:     " << m_numberOfDimensions << std::endl;
    std::cout << "Interaction style:        " << m_interactionStyle->getLabel() << std::endl;
    std::cout << "Hamiltonian:              " << m_hamiltonian->getLabel() << std::endl;
    std::cout << "Initial state:            " << m_initialState->getLabel() << std::endl;
    if (m_hamiltonian->getLabel() == "harmonic oscillator") {
        std::cout << "Oscillator frequency:     " << m_omega << std::endl;
    } else if (m_hamiltonian->getLabel() == "double well"){
        std::cout << "Oscillator frequency:     " << m_omega << std::endl;
    } else if (m_hamiltonian->getLabel() == "atom") {
        std::cout << "Atomic number:            " << m_Z << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "           WAVE FUNCTION INFORMATION" << std::endl;
    std::cout << "==============================================" << std::endl;
    for (int i = 0; i < m_numberOfElements; i++) {
        std::cout << "Element " << i << ":                " << m_waveFunctionElements[unsigned(i)]->getLabel() << std::endl;
    }
    std::cout << "Basis:                    " << m_basis->getLabel() << std::endl;
    std::cout << "Initial parameters:       " << m_initialWeights->getLabel() << std::endl;
    std::cout << "Number of parameters:     " << this->getTotalNumberOfParameters() << std::endl;
    std::cout << "Number of hidden nodes:   " << m_numberOfHiddenUnits << std::endl;
    std::cout << "Print parameters to file: " << m_printParametersToFile << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "            SIMULATION INFORMATION" << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Max number of iterations: " << m_numberOfIterations << std::endl;
    std::cout << "Number of MC cycles:      " << m_totalStepsWOEqui << std::endl;
    std::cout << "Number of burn-in cycles: " << m_equilibriationSteps << std::endl;
    std::cout << "Learning rate:            " << m_eta << std::endl;
    std::cout << "Step length:              " << m_stepLength << std::endl;
    std::cout << "Equilibration fraction:   " << m_equilibrationFraction << std::endl;
    std::cout << "Sampling:                 " << m_metropolis->getLabel() << std::endl;
    std::cout << "Optimization:             " << m_optimization->getLabel() << std::endl;
    std::cout << "Random number generator:  " << m_randomNumberGenerator->getLabel() << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "               PARTICLE DENSITY" << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Compute radial one-body density:  " << m_computeOneBodyDensity << std::endl;
    std::cout << "Compute spatial one-body density: " << m_computeOneBodyDensity2 << std::endl;
    std::cout << "Compute radial two-body density:  " << m_computeTwoBodyDensity << std::endl;
    std::cout << "Max radius of electron density:   " << m_maxRadius << std::endl;
    std::cout << "Number of bins in each direction: " << m_numberOfBins << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "                  RESAMPLING" << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Do resampling:            " << m_doResampling << std::endl;
    std::cout << "Print energies to file:   " << m_printEnergyToFile << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "                  CONVERGENCE" << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Check convergence:        " << m_checkConvergence << std::endl;
    std::cout << "Tolerance:                " << m_tolerance << std::endl;
    std::cout << "Number of energies:       " << m_numberOfEnergies << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "                 ADAPTIVE STEPS" << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Apply adaptive steps:       " << m_applyAdaptiveSteps << std::endl;
    std::cout << "Range of adaptive steps:    " << m_rangeOfAdaptiveSteps << std::endl;
    std::cout << "Additional steps:           " << m_additionalSteps << std::endl;
    std::cout << "Additional steps last iter: " << m_additionalStepsLastIter << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
}


/* ----------------------------------------------------------------------------
  Print header line to terminal, containing information about the various 
  columns
---------------------------------------------------------------------------- */

void System::printHeaderLine()
{
    std::cout << "Step  " << "Energy  " << "Energy_STD  " << "Kinetic  "
              << "Kinetic_STD  " << "External  " << "External_STD  "
              << "Interaction  " << "Interaction_STD  " << "Acceptence  "
              << "CPU_time" << std::endl;
}


/* ----------------------------------------------------------------------------
  Run simulation specified in main/configuration file. This the main loop in
  VMC, where we update parameters
---------------------------------------------------------------------------- */

void System::runSimulation()
{
    if (m_rank == 0) {
        printLogo();
        printInitialInformation();
    }
    initializeSystem();
    if (m_rank == 0) {
        printSystemInformation();
        printHeaderLine();
    }
    m_numberOfNormalIterations = m_numberOfIterations - m_rangeOfAdaptiveSteps - 1;
    for (m_iter = 0; m_iter < m_numberOfIterations; m_iter++) {
        if (m_applyAdaptiveSteps) {
            m_stepsWOEqui = m_initialStepsWOEqui * adaptiveSteps();
            m_stepsWEqui = m_stepsWOEqui + m_equilibriationSteps;
            m_totalStepsWOEqui = m_initialTotalStepsWOEqui * adaptiveSteps();
            m_totalStepsWEqui = m_totalStepsWOEqui + m_totalEquilibriationSteps;
        }
        m_sampler->setNumberOfSteps(m_stepsWOEqui, m_totalStepsWOEqui, m_totalStepsWEqui);
        double startTime = MPI_Wtime();
        runMetropolisCycles();
        double endTime = MPI_Wtime();
        double time = endTime - startTime;

        MPI_Reduce(&time, &m_totalTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (m_iter < m_numberOfNormalIterations) {
            m_globalTime += m_totalTime;
        }
        m_sampler->computeTotals();

        if (m_rank == 0) {
            m_sampler->computeAverages();
            m_parameters -= m_optimization->updateParameters();
        }
        m_sampler->printParametersToFile();
        m_sampler->printEnergyToFile();
        if (m_iter == m_numberOfNormalIterations + m_rangeOfAdaptiveSteps) {
            m_sampler->printOneBodyDensityToFile();
            m_sampler->printOneBodyDensity2ToFile();
            m_sampler->printTwoBodyDensityToFile();
        }
        printToTerminal();

        if (m_checkConvergence && m_rank == 0) {
            checkingConvergence();
        }
        MPI_Bcast(&m_numberOfNormalIterations, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(m_parameters.data(),
                  int(m_numberOfElements * m_maxParameters),
                  MPI_DOUBLE,
                  0,
                  MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        updateAllParameters(m_parameters);
    }
}


/* ----------------------------------------------------------------------------
  This is the sampling loop, where we sample the particles with Markov chains
  and compute expectation values
---------------------------------------------------------------------------- */

void System::runMetropolisCycles()
{
    for (int i = 0; i < m_stepsWEqui; i++) {
        bool acceptedStep = m_metropolis->acceptMove();
        m_positions = m_metropolis->updatePositions();
        m_distanceMatrix = m_metropolis->updateDistanceMatrix();
        m_radialVector = m_metropolis->updateRadialVector();
        if (i >= m_equilibriationSteps) {
            m_sampler->sample(acceptedStep, i);
            if (m_iter == m_numberOfNormalIterations + m_rangeOfAdaptiveSteps) {
                m_sampler->printInstantValuesToFile();
                m_sampler->computeOneBodyDensity(m_radialVector);
                m_sampler->computeTwoBodyDensity(m_radialVector);
                m_sampler->computeOneBodyDensity2(m_positions);
            }
        }
    }
}


/* ----------------------------------------------------------------------------
  Call this function to print to terminal
---------------------------------------------------------------------------- */

void System::printToTerminal()
{
    if (m_iter == m_numberOfNormalIterations + m_rangeOfAdaptiveSteps) {
        m_sampler->closeOutputFiles();
        if (m_rank == 0) {
            m_sampler->doResampling();
            m_sampler->printFinalOutputToTerminal(m_start);
            std::cout << std::endl;
            std::cout << "Average CPU time before convergence: " << m_globalTime / m_numberOfNormalIterations
                      << std::endl;
        }

        MPI_Finalize();
        exit(0);
    } else {
        if (m_rank == 0) {
            m_sampler->printOutputToTerminal(m_numberOfIterations, m_totalTime);
        }
    }
}


/* ----------------------------------------------------------------------------
  This function checks if the trial wave function has converged
---------------------------------------------------------------------------- */

void System::checkingConvergence()
{
    m_energies.head(m_numberOfEnergies - 1) = m_energies.tail(m_numberOfEnergies - 1);
    m_energies(m_numberOfEnergies - 1) = m_sampler->getAverageEnergy();
    if (fabs(m_energies(0) - m_energies(m_numberOfEnergies - 1)) < m_tolerance) {
        std::cout << "The system has converged! Let's run one more cycle to collect data"
                  << std::endl;
        m_numberOfNormalIterations = m_iter + 1;
        m_checkConvergence = false;
    }
}


/* ----------------------------------------------------------------------------
  Is responsible for the adaptive steps
---------------------------------------------------------------------------- */

int System::adaptiveSteps()
{
    int stepRatio = 1;
    if (m_iter == m_numberOfNormalIterations + m_rangeOfAdaptiveSteps) {
        stepRatio = int(pow(2, m_additionalStepsLastIter));

    } else if (m_iter >= m_numberOfNormalIterations) {
        stepRatio = int(pow(2, m_additionalSteps));
    }
    return stepRatio;
}


// OPERATIONS ON WAVE FUNCTION ELEMENTS

/* ----------------------------------------------------------------------------
  Initialize the wave function elements with essential variables specified by
  used
---------------------------------------------------------------------------- */

void System::setAllConstants()
{
    for (int i = 0; i < m_numberOfElements; i++) {
        m_waveFunctionElements[unsigned(i)]->setConstants(i);
    }
}


/* ----------------------------------------------------------------------------
  Initialize the wave function elements with position
---------------------------------------------------------------------------- */

void System::initializeAllArrays(const Eigen::VectorXd positions,
                                 const Eigen::VectorXd radialVector,
                                 const Eigen::MatrixXd distanceMatrix)
{
    for (auto &i : m_waveFunctionElements) {
        i->initializeArrays(positions, radialVector, distanceMatrix);
        i->setArrays();
    }
}


/* ----------------------------------------------------------------------------
  Update positions in all wave function elements when a particle is moved
---------------------------------------------------------------------------- */

void System::updateAllArrays(const Eigen::VectorXd positions,
                             const Eigen::VectorXd radialVector,
                             const Eigen::MatrixXd distanceMatrix,
                             const int changedCoord)
{
    for (auto &i : m_waveFunctionElements) {
        i->setArrays();
        i->updateArrays(positions, radialVector, distanceMatrix, changedCoord);
    }
}


/* ----------------------------------------------------------------------------
  Reset positions in all ave function elements when a particle move is rejected
---------------------------------------------------------------------------- */

void System::resetAllArrays()
{
    for (auto &i : m_waveFunctionElements) {
        i->resetArrays();
    }
}


/* ----------------------------------------------------------------------------
  Update the parameters/weights in all the wave function element when the
  parameters are updated
---------------------------------------------------------------------------- */

void System::updateAllParameters(const Eigen::MatrixXd parameters)
{
    for (auto &i : m_waveFunctionElements) {
        i->updateParameters(parameters);
    }
}


/* ----------------------------------------------------------------------------
  Evaluate the collective probability ratio
---------------------------------------------------------------------------- */

double System::evaluateProbabilityRatio()
{
    double ratio = 1;
    for (auto &i : m_waveFunctionElements) {
        ratio *= i->evaluateRatio();
    }
    return ratio;
}


/* ----------------------------------------------------------------------------
  Obtain the total kinetic energy of the system based on the gradients and 
  Laplacians of all elements
---------------------------------------------------------------------------- */

double System::getKineticEnergy()
{
    double kineticEnergy = 0;
    for (auto &i : m_waveFunctionElements) {
        kineticEnergy += i->computeLaplacian();
    }
    for (int k = 0; k < m_degreesOfFreedom; k++) {
        double nablaLnPsi = 0;
        for (auto &i : m_waveFunctionElements) {
            nablaLnPsi += i->computeGradient(k);
        }
        kineticEnergy += nablaLnPsi * nablaLnPsi;
    }
    return -0.5 * kineticEnergy;
}


/* ----------------------------------------------------------------------------
  Get the gradient of all the wave function elements with respect to all the 
  parameters. To be used in the parameter update
---------------------------------------------------------------------------- */

Eigen::MatrixXd System::getAllParameterGradients()
{
    for (int i = 0; i < m_numberOfElements; i++) {
        m_gradients.row(i) = m_waveFunctionElements[unsigned(i)]->computeParameterGradient();
    }
    return m_gradients;
}


/* ----------------------------------------------------------------------------
  Check if the elements need distance matrix and/or radial distance vector
---------------------------------------------------------------------------- */

void System::setGlobalArraysToCalculate()
{
    // Check if the elements need distance matrix and/or radial distance vector
    for (auto &p : m_waveFunctionElements) {
        int need = p->getGlobalArrayNeed();
        if (need == 1) {
            m_calculateDistanceMatrix = true;
        } else if (need == 2) {
            m_calculateRadialVector = true;
        } else if (need == 3) {
            m_calculateDistanceMatrix = true;
            m_calculateRadialVector = true;
        }
    }
    // Check if the Hamiltonian needs distance matrix and/or radial distance vector
    int need = m_hamiltonian->getGlobalArrayNeed();
    if (need == 1) {
        m_calculateDistanceMatrix = true;
    } else if (need == 2) {
        m_calculateRadialVector = true;
    } else if (need == 3) {
        m_calculateDistanceMatrix = true;
        m_calculateRadialVector = true;
    }
    need = m_interactionStyle->getGlobalArrayNeed();
    if (need == 1) {
        m_calculateDistanceMatrix = true;
    } else if (need == 2) {
        m_calculateRadialVector = true;
    } else if (need == 3) {
        m_calculateDistanceMatrix = true;
        m_calculateRadialVector = true;
    }
}


/* ----------------------------------------------------------------------------
  Set the maximum number of parameters found in a wave function element. This
  will be used to create the parameter matrix, with dim
  (maxNumberOfParameter x numberOfElements). Determines also the total number
  of particles. Called automatically when needed
---------------------------------------------------------------------------- */

void System::setMaxParameters()
{
    int maxNumberOfElements = 0;
    int counter = 0;
    for (auto &i : m_waveFunctionElements) {
        int numberOfParameters = i->getNumberOfParameters();
        if (numberOfParameters > maxNumberOfElements) {
            maxNumberOfElements = numberOfParameters;
        }
        counter += numberOfParameters;
    }
    m_maxParameters = maxNumberOfElements;
    m_totalNumberOfParameters = counter;
}

// VARIABLES SPECIFIED IN MAIN

/* ----------------------------------------------------------------------------
  Specify the number of particles, N, to be used
---------------------------------------------------------------------------- */

void System::setNumberOfParticles(const int numberOfParticles)
{
    assert(numberOfParticles > 0);
    m_numberOfParticles = numberOfParticles;
    initializeMPI();
}


/* ----------------------------------------------------------------------------
  Specify the number of dimensions, d, of the system. Supported 2 and 3
  dimensions
---------------------------------------------------------------------------- */

void System::setNumberOfDimensions(const int numberOfDimensions)
{
    assert(numberOfDimensions > 1);
    assert(numberOfDimensions < 4);
    m_numberOfDimensions = numberOfDimensions;
    setNumberOfFreeDimensions();
}


/* ----------------------------------------------------------------------------
  Set the number of free dimensions, F=Nd. This is called automatically after
  the number of particles and number of dimensions are specified
---------------------------------------------------------------------------- */

void System::setNumberOfFreeDimensions()
{
    m_degreesOfFreedom = m_numberOfParticles * m_numberOfDimensions;
}


/* ----------------------------------------------------------------------------
  Set the number of hidden units used in RBM trial wave function
---------------------------------------------------------------------------- */

void System::setNumberOfHiddenUnits(const int numberOfHiddenUnits)
{
    assert(numberOfHiddenUnits > 0);
    m_numberOfHiddenUnits = numberOfHiddenUnits;
}


/* ----------------------------------------------------------------------------
  Determine the number of steps used by each process when the equilibriation
  (burn-in period) is excluded (a power of 2) and when it is included
---------------------------------------------------------------------------- */

void System::setNumberOfMetropolisCycles(const int steps)
{
    m_totalStepsWOEqui = steps;
    if (m_rank == 0) {
        m_stepsWOEqui = steps / m_numberOfProcesses + steps % m_numberOfProcesses;
    } else {
        m_stepsWOEqui = steps / m_numberOfProcesses;
    }

    // Store the initial steps in case adaptive step is chosen
    m_initialStepsWOEqui = m_stepsWOEqui;
    m_initialTotalStepsWOEqui = m_totalStepsWOEqui;

    // Calculate the number of equilibriation steps (needs to be unaffected by the number of processes)
    m_equilibriationSteps = int(m_totalStepsWOEqui * m_equilibrationFraction);
    m_totalEquilibriationSteps = int(m_totalStepsWOEqui * m_equilibrationFraction
                                     * m_numberOfProcesses);

    // Calculate the number of steps included equilibriation
    m_totalStepsWEqui = m_totalStepsWOEqui + m_totalEquilibriationSteps;
    m_stepsWEqui = m_stepsWOEqui + m_equilibriationSteps;
}


/* ----------------------------------------------------------------------------
  Set the number of wave function elements. This is called automatically when
  needed
---------------------------------------------------------------------------- */

void System::setNumberOfElements(const unsigned long numberOfElements)
{
    m_numberOfElements = static_cast<int>(numberOfElements);
    collectAllLabels();
}


/* ----------------------------------------------------------------------------
  Set maximum number of iterations used in the Monte Carlo integration. The
  actual number of iterations used can be overruled by the stop criterion
---------------------------------------------------------------------------- */

void System::setNumberOfIterations(const int numberOfIterations)
{
    assert(numberOfIterations > 0);
    m_numberOfIterations = numberOfIterations;
}


/* ----------------------------------------------------------------------------
  Set step length used in sampling
---------------------------------------------------------------------------- */

void System::setStepLength(const double stepLength)
{
    assert(stepLength > 0);
    m_stepLength = stepLength;
}


/* ----------------------------------------------------------------------------
  Set equilibriation fraction (burn-in period)
---------------------------------------------------------------------------- */

void System::setEquilibrationFraction(const double equilibrationFraction)
{
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
    setNumberOfMetropolisCycles(m_totalStepsWOEqui);
}


/* ----------------------------------------------------------------------------
  Set frequency of harmonic oscillator potential
---------------------------------------------------------------------------- */

void System::setFrequency(const double omega)
{
    assert(omega > 0);
    m_omega = omega;
}


/* ----------------------------------------------------------------------------
  By calling this function, the screening is enabled. The screening strength
  and Debye screening length (DSL) need to be specified
---------------------------------------------------------------------------- */

void System::setScreeningTools(const double screeningStrength,
                               const double dsl)
{
    assert(screeningStrength >= 1);
    assert(dsl > 0);
    m_screening = true;
    m_screeningStrength = screeningStrength;
    m_dsl = dsl;
}


/* ----------------------------------------------------------------------------
  Set atomic number, Z, when simulating atoms
---------------------------------------------------------------------------- */

void System::setAtomicNumber(const int Z)
{
    assert(Z > 0);
    m_Z = Z;
}


/* ----------------------------------------------------------------------------
  Set learning rate of simulation
---------------------------------------------------------------------------- */

void System::setLearningRate(const double eta)
{
    assert(eta > 0);
    m_eta = eta;
}


/* ----------------------------------------------------------------------------
  Set distribution width of Gaussian distribution when using a Gaussian-binary
  restricted Boltzmann machine as the trial wave function guess
---------------------------------------------------------------------------- */

void System::setWidth(const double sigma)
{
    assert(sigma > 0);
    m_sigma = sigma;
}


/* ----------------------------------------------------------------------------
  Specify convergence. Convergence is enabled when this function is called
---------------------------------------------------------------------------- */

void System::setConvergenceTools(int numberOfEnergies, double tolerance)
{
    m_checkConvergence = true;
    m_tolerance = tolerance;
    m_numberOfEnergies = numberOfEnergies;
    m_energies = Eigen::VectorXd::Zero(numberOfEnergies);
}


/* ----------------------------------------------------------------------------
  Set adaptive steps. Adaptive steps are enabled when this function is called.
  Default is 4 times as many cycles when 10 iterations are remaining and 8
  times as many cycles for the last iteration
---------------------------------------------------------------------------- */

void System::setAdaptiveStepTools(int rangeOfAdaptiveSteps,
                                  int additionalSteps,
                                  int additionalStepsLastIteration)
{
    m_applyAdaptiveSteps = true;
    m_rangeOfAdaptiveSteps = rangeOfAdaptiveSteps;
    m_additionalSteps = additionalSteps;
    m_additionalStepsLastIter = additionalStepsLastIteration;
}


/* ----------------------------------------------------------------------------
  Compute the radial one-body density. Default number of bins is 1000, default
  max radius is 50 given in natural units
---------------------------------------------------------------------------- */

void System::computeRadialOneBodyDensity(int numberOfBins, double maxRadius)
{
    m_computeOneBodyDensity = true;
    m_numberOfBins = numberOfBins;
    m_maxRadius = maxRadius;
}


/* ----------------------------------------------------------------------------
  Compute the spatial one-body density. Default number of bins is 1000 x 1000,
  default max radius is 50 given in natural units
---------------------------------------------------------------------------- */

void System::computeSpatialOneBodyDensity(int numberOfBins, double maxRadius)
{
    m_computeOneBodyDensity2 = true;
    m_numberOfBins = numberOfBins;
    m_maxRadius = maxRadius;
}


/* ----------------------------------------------------------------------------
  Compute the radial two-body density. Default number of bins is 1000, default
  max radius is 50 given in natural units
---------------------------------------------------------------------------- */

void System::computeTwoBodyDensity(int numberOfBins, double maxRadius)
{
    m_computeTwoBodyDensity = true;
    m_numberOfBins = numberOfBins;
    m_maxRadius = maxRadius;
}


/* ----------------------------------------------------------------------------
  Set total spin of the system. To be implemented
---------------------------------------------------------------------------- */

void System::setTotalSpin(const double totalSpin)
{
    double intpart;
    assert(std::modf(m_numberOfParticles / 2 - abs(totalSpin), &intpart) == 0.0);
    m_totalSpin = totalSpin;
}


/* ----------------------------------------------------------------------------
  Dump the energy expectation value after each iteration to fule. Default 
  setting is to print energy file
---------------------------------------------------------------------------- */

void::System::dumpEnergyToFile(bool printEnergyFile) {
    m_printEnergyToFile = printEnergyFile;
}


/* ----------------------------------------------------------------------------
  Dump the instant energy to file (energies from all the cycles). This well be
  used by resampling. Default setting is to print instant file for the last
  iteration
---------------------------------------------------------------------------- */

void::System::doBlocking(bool printInstantEnergyFile) {
    m_doResampling = printInstantEnergyFile;
}


/* ----------------------------------------------------------------------------
  Print parameters to file after each iteration. Makes it possible to start
  where the previous simulation ended. Default setting is to print parameters
  to file
---------------------------------------------------------------------------- */

void System::dumpParametersToFile(bool printParametersToFile)
{
    m_printParametersToFile = printParametersToFile;
}


/* ----------------------------------------------------------------------------
  Set path to where to save the files
---------------------------------------------------------------------------- */

void System::setPath(const std::string path)
{
    m_path = path;
}


// CLASSES SPECIFIED IN MAIN


/* ----------------------------------------------------------------------------
  Specify which Hamiltonian to be used
---------------------------------------------------------------------------- */

void System::setHamiltonian(Hamiltonian *hamiltonian)
{
    m_hamiltonian = hamiltonian;
}


/* ----------------------------------------------------------------------------
  Specify the basis to be used in the Slater determinant
---------------------------------------------------------------------------- */

void System::setBasis(Basis *basis)
{
    m_basis = basis;
}


/* ----------------------------------------------------------------------------
  Set all wave function elements, based on an std::vector
---------------------------------------------------------------------------- */

void System::setWaveFunctionElements(std::vector<class WaveFunction *> waveFunctionElements)
{
    m_waveFunctionElements = waveFunctionElements;
    setNumberOfElements(waveFunctionElements.size());
    setAllConstants();
}


/* ----------------------------------------------------------------------------
  Add a wave function element
---------------------------------------------------------------------------- */

void System::setWaveFunctionElement(WaveFunction *waveFunction)
{
    m_waveFunctionElements.push_back(waveFunction);
    setNumberOfElements(m_waveFunctionElements.size());
}


/* ----------------------------------------------------------------------------
  Set the input layer then a feed-forward neural network (FNN) is used as a
  trial wave function element
---------------------------------------------------------------------------- */

void System::setInputLayer(int numberOfUnits)
{
    m_layers.insert(m_layers.begin(), new Input(this, numberOfUnits));
}


/* ----------------------------------------------------------------------------
  Add dense layer to a feed-forward neural network
---------------------------------------------------------------------------- */

void System::addDenseLayer(int numberOfUnits, Activation *activation)
{
    m_layers.push_back(new Dense(this, numberOfUnits, activation));
    setNumberOfElements(m_waveFunctionElements.size());
}


/* ----------------------------------------------------------------------------
  Add dense layer to a feed-forward neural network
---------------------------------------------------------------------------- */

void System::setOutputLayer(Activation *activation)
{
    m_layers.push_back(new Output(this, activation));
    setNumberOfElements(m_waveFunctionElements.size());
}


/* ----------------------------------------------------------------------------
  Initialize particle positions. Possible choices are:
   - Random normal (RandomNormal, default)
   - Random uniform (RandomUniform)
---------------------------------------------------------------------------- */

void System::setInitialState(InitialState *initialState)
{
    m_initialState = initialState;
}


/* ----------------------------------------------------------------------------
  Initialize parameters. Possible choices are:
    - Automatize (Automatize, default)
    - Constant (ConstantParameters)
    - Random uniform (RandomUniform)
---------------------------------------------------------------------------- */

void System::setInitialWeights(InitialWeights *initialWeights)
{
    m_initialWeights = initialWeights;
}


/* ----------------------------------------------------------------------------
  Set interaction style
---------------------------------------------------------------------------- */

void System::setInteractionStyle(Interaction *interaction)
{
    m_interactionStyle = interaction;
}


/* ----------------------------------------------------------------------------
  Specify sampling algorithm. Possible choices:
    - Brute force (BruteForce)
    - Importance sampling (ImportanceSampling)
---------------------------------------------------------------------------- */
  
void System::setMetropolis(Metropolis *metropolis)
{
    m_metropolis = metropolis;
}


/* ----------------------------------------------------------------------------
  Specify optimization tool when updating parameters. Possible choices:
    - Gradient descent (GD)
    - Stochastic gradient descent (SGD)
    - Adam (ADAM)
---------------------------------------------------------------------------- */
  
void System::setOptimization(Optimization *optimization)
{
    m_optimization = optimization;
}


/* ----------------------------------------------------------------------------
  Specify which randum number generator to be used. Currently, only Mersenne
  Twister 19337 is available
---------------------------------------------------------------------------- */

void System::setRandomNumberGenerator(RandomNumberGenerator *randomNumberGenerator)
{
    m_randomNumberGenerator = randomNumberGenerator;
}


/* ----------------------------------------------------------------------------
  Initialize gradient matrix
---------------------------------------------------------------------------- */

void System::setGradients()
{
    m_gradients = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxParameters);
}

/* ----------------------------------------------------------------------------
  Initialize variables from configuration file. Will overwrite main and
  default settings
---------------------------------------------------------------------------- */

void System::initializeFromConfig(int argc, char** argv) {
    if (argc >= 2) {
        m_configFile = argv[1];
        m_args = argc;
    }
}

/* ----------------------------------------------------------------------------
  Parse input script
---------------------------------------------------------------------------- */

void System::parser(const std::string configFile)
{
    std::ifstream infile;
    infile.open(configFile.c_str());
    if (!infile.is_open() && m_args >= 2) {
        std::cout << std::endl;
        std::cerr << "File: '" << configFile << "'" << std::endl;
        perror("File not found");
        MPI_Abort(MPI_COMM_WORLD, 143);
    }
    else {
        std::string line;
        while (std::getline(infile, line)) {
            std::istringstream is_line(line);
            std::string key;
            if (line.empty()) {
                /* Continue if line is blank */
                continue;
            } else if (line.rfind("#", 0) == 0) {
                /* Continue if line starts with '#' */
                continue;
            } else if (std::getline(is_line, key, ':')) {
                key = trim(key);
                std::string value;
                if (std::getline(is_line, value)) {
                    std::vector<std::string> splitted = split(value);
                    if (key == "numParticles") {
                        m_numberOfParticles = std::stoi(splitted.at(0));
                        m_numberOfHiddenUnits = m_numberOfParticles;
                        m_Z = m_numberOfParticles;
                    } else if (key == "numDimensions") {
                        m_numberOfDimensions = std::stoi(splitted.at(0));
                    } else if (key == "omega") {
                        m_omega = std::stod(splitted.at(0));
                        m_stepLength = 0.1 / sqrt(m_omega);
                        m_sigma = 1.0 / sqrt(m_omega);
                    } else if (key == "atomicNumber") {
                        m_Z = std::stoi(splitted.at(0));
                    } else if (key == "learningRate") {
                        m_eta = std::stod(splitted.at(0));
                    } else if (key == "maxRadius") {
                        m_maxRadius = std::stod(splitted.at(0));
                    } else if (key == "numIterations") {
                        m_numberOfIterations = std::stoi(splitted.at(0));
                    } else if (key == "equilibration") {
                        setEquilibrationFraction(std::stod(splitted.at(0)));
                    } else if (key == "numSteps") {
                        setNumberOfMetropolisCycles(std::stoi(splitted.at(0)));
                    } else if (key == "numHiddenNodes") {
                        m_numberOfHiddenUnits = std::stoi(splitted.at(0));
                    } else if (key == "totalSpin") {
                        m_totalSpin = std::stod(splitted.at(0));
                    } else if (key == "stepLength") {
                        m_stepLength = std::stod(splitted.at(0));
                    } else if (key == "checkConvergence") {
                        std::istringstream(splitted.at(0)) >> std::boolalpha >> m_checkConvergence;
                    } else if (key == "applyAdaptiveSteps") {
                        std::istringstream(splitted.at(0)) >> std::boolalpha >> m_applyAdaptiveSteps;
                    } else if (key == "computeOneBodyDensity") {
                        std::istringstream(splitted.at(0)) >> std::boolalpha >> m_computeOneBodyDensity;
                    } else if (key == "computeOneBodyDensity2") {
                        std::istringstream(splitted.at(0)) >> std::boolalpha >> m_computeOneBodyDensity2;
                    } else if (key == "computeTwoBodyDensity") {
                        std::istringstream(splitted.at(0)) >> std::boolalpha >> m_computeTwoBodyDensity;
                    } else if (key == "printEnergyToFile") {
                        std::istringstream(splitted.at(0)) >> std::boolalpha >> m_printEnergyToFile;
                    } else if (key == "printParametersToFile") {
                        std::istringstream(splitted.at(0)) >> std::boolalpha >> m_printParametersToFile;
                    } else if (key == "doResampling") {
                        std::istringstream(splitted.at(0)) >> std::boolalpha >> m_doResampling;
                    } else if (key == "path") {
                        m_path = splitted.at(0);
                    } else if (key == "numberOfEnergies") {
                        m_numberOfEnergies = std::stoi(splitted.at(0));
                    } else if (key == "tolerance") {
                        m_tolerance = std::stod(splitted.at(0));
                    } else if (key == "rangeOfAdaptiveSteps") {
                        m_rangeOfAdaptiveSteps = std::stoi(splitted.at(0));
                    } else if (key == "additionalSteps") {
                        m_additionalSteps = std::stoi(splitted.at(0));
                    } else if (key == "additionalStepsLastIter") {
                        m_additionalStepsLastIter = std::stoi(splitted.at(0));
                    } else if (key == "numberOfBins") {
                        m_numberOfBins = std::stoi(splitted.at(0));
                    } else if (key == "basis") {
                        delete m_basis;
                        if (splitted.at(0) == "hermite") {
                            setBasis(new Hermite(this));
                        //} else if (splitted.at(0) == "hermiteExpansion") {
                        //    setBasis(new HermiteExpansion(this));
                        } else if (splitted.at(0) == "hydrogenOrbital") {
                            setBasis(new HydrogenOrbital(this));
                        } else {
                            std::cout << std::endl;
                            std::cerr << "Basis '"
                                      << splitted.at(0)
                                      << "' is not implemented!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                    } else if (key == "hamiltonian") {
                        delete m_hamiltonian;
                        if (splitted.at(0) == "harmonicOscillator") {
                            setHamiltonian(new HarmonicOscillator(this));
                        } else if (splitted.at(0) == "doubleWell") {
                            setHamiltonian(new DoubleWell(this, std::stod(splitted.at(1))));
                        } else if (splitted.at(0) == "atomicNucleus") {
                            setHamiltonian(new AtomicNucleus(this));
                        } else if (splitted.at(0) == "ellipticalHarmonicOscillator") {
                            setHamiltonian(new EllipticalHarmonicOscillator(this, std::stod(splitted.at(1))));
                        } else {
                            std::cout << std::endl;
                            std::cerr << "The Hamiltonian '"
                                      << splitted.at(0)
                                      << "' does not exist!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                    } else if (key == "optimization") {
                        delete m_optimization;
                        if (splitted.at(0) == "adam") {
                            setOptimization(new ADAM(this));
                        } else if (splitted.at(0) == "gd") {
                            if (splitted.size() >= 3) {
                                setOptimization(new GradientDescent(this, std::stod(splitted.at(1)), std::stod(splitted.at(2))));
                            } else {
                                setOptimization(new GradientDescent(this, 0.0, 0.0));
                            }
                        } else if (splitted.at(0) == "sgd") {
                          if (splitted.size() >= 3) {
                              setOptimization(new SGD(this, std::stod(splitted.at(1)), std::stod(splitted.at(2))));
                          } else {
                              setOptimization(new SGD(this, 0.0, 0.0));
                          }
                        } else {
                            std::cout << std::endl;
                            std::cerr << "Optimization method '"
                                      << splitted.at(0)
                                      << "' does not exist!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                    } else if (key == "initialWeights") {
                        delete m_initialWeights;
                        if (splitted.at(0) == "automatize") {
                            setInitialWeights(new Automatize(this));
                        } else if (splitted.at(0) == "randomuniform") {
                            if (splitted.size() >= 2) {
                                setInitialWeights(new RandomUniformWeights(this, std::stod(splitted.at(1))));
                            } else {
                                setInitialWeights(new RandomUniformWeights(this));
                            }
                        } else if (splitted.at(0) == "constant") {
                            if (splitted.size() >= 2) {
                                setInitialWeights(new Constant(this, std::stod(splitted.at(1))));
                            } else {
                                setInitialWeights(new Constant(this, 1.0));
                            }
                        } else {
                            std::cout << std::endl;
                            std::cerr << "Initial parameter configuration '"
                                      << splitted.at(0)
                                      << "' does not exist!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                    } else if (key == "initialState") {
                        delete m_initialState;
                        if (splitted.at(0) == "randomNormal") {
                            setInitialState(new RandomNormal(this));
                        } else if (splitted.at(0) == "randomUniform") {
                            setInitialState(new RandomUniform(this));
                        } else {
                            std::cout << std::endl;
                            std::cerr << "Initial state configuration '"
                                      << splitted.at(0)
                                      << "' does not exist!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                    } else if (key == "sampling") {
                        delete m_metropolis;
                        if (splitted.at(0) == "importanceSampling") {
                            setMetropolis(new ImportanceSampling(this));
                        } else if (splitted.at(0) == "bruteForce") {
                            setMetropolis(new BruteForce(this));
                        } else {
                            std::cout << std::endl;
                            std::cerr << "Sampling method '"
                                      << splitted.at(0)
                                      << "' does not exist!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                    } else if (key == "waveFunction") {
                        for (int i=m_numberOfElements; i--;)
                        {
                            delete m_waveFunctionElements[i];
                        }
                        m_waveFunctionElements.clear();

                        std::vector<WaveFunction *> waveFunctionElements;
                        if (splitted.at(0) == "VMC") {
                            m_waveFunctionElements.push_back(new class Gaussian(this));
                            //waveFunctionElements.push_back(new class Gaussian(this));
                            waveFunctionElements.push_back(new class SlaterDeterminant(this));
                            waveFunctionElements.push_back(new class PadeJastrow(this));
                        } else if (splitted.at(0) == "RBM") {
                            waveFunctionElements.push_back(new class SlaterDeterminant(this));
                            waveFunctionElements.push_back(new class RBMGaussian(this));
                            waveFunctionElements.push_back(new class RBMProduct(this));
                        } else if (splitted.at(0) == "RBMPJ") {
                            waveFunctionElements.push_back(new class SlaterDeterminant(this));
                            waveFunctionElements.push_back(new class RBMGaussian(this));
                            waveFunctionElements.push_back(new class RBMProduct(this));
                            waveFunctionElements.push_back(new class PadeJastrow(this));
                        } else if (splitted.at(0) == "RBMSJ") {
                            waveFunctionElements.push_back(new class SlaterDeterminant(this));
                            waveFunctionElements.push_back(new class RBMGaussian(this));
                            waveFunctionElements.push_back(new class RBMProduct(this));
                            waveFunctionElements.push_back(new class SimpleJastrow(this));
                        } else if (splitted.at(0) == "PRBM") {
                            waveFunctionElements.push_back(new class SlaterDeterminant(this));
                            waveFunctionElements.push_back(new class RBMGaussian(this));
                            waveFunctionElements.push_back(new class RBMProduct(this));
                            waveFunctionElements.push_back(new class PartlyRestricted(this));
                        } else {
                            std::cout << std::endl;
                            std::cerr << "Wave function '"
                                      << splitted.at(0)
                                      << "' does not exist!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                        setWaveFunctionElements(waveFunctionElements);
                    } else if (key == "waveFunctionElement") {
                        if (splitted.at(0) == "gaussian") {
                            setWaveFunctionElement(new class Gaussian(this));
                        } else if (splitted.at(0) == "slaterDeterminant") {
                            setWaveFunctionElement(new class SlaterDeterminant(this));
                        } else if (splitted.at(0) == "padeJastrow") {
                            setWaveFunctionElement(new class PadeJastrow(this));
                        } else if (splitted.at(0) == "simpleJastrow") {
                            setWaveFunctionElement(new class SimpleJastrow(this));
                        } else if (splitted.at(0) == "RBMGaussian") {
                            setWaveFunctionElement(new class RBMGaussian(this));
                        } else if (splitted.at(0) == "RBMProduct") {
                            setWaveFunctionElement(new class RBMProduct(this));
                        } else if (splitted.at(0) == "hydrogenLike") {
                            setWaveFunctionElement(new class HydrogenLike(this));
                        } else if (splitted.at(0) == "hardCoreJastrow") {
                            setWaveFunctionElement(new class HardCoreJastrow(this));
                        } else {
                            std::cout << std::endl;
                            std::cerr << "Wave function element '"
                                      << splitted.at(0)
                                      << "' does not exist!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                    } else if (key == "interactionStyle") {
                        delete m_interactionStyle;
                        if (splitted.at(0) == "noInteraction") {
                            setInteractionStyle(new class NoInteraction(this));
                        } else if (splitted.at(0) == "coulomb") {
                            setInteractionStyle(new class Coulomb(this));
                        } else {
                            std::cout << std::endl;
                            std::cerr << "Interaction style '"
                                      << splitted.at(0)
                                      << "' does not exist!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                    } else if (key == "rng") {
                        delete m_randomNumberGenerator;
                        if (splitted.at(0) == "MersenneTwister") {
                            setRandomNumberGenerator(new class MersenneTwister());
                        } else {
                            std::cout << std::endl;
                            std::cerr << "Random number generator '"
                                      << splitted.at(0)
                                      << "' does not exist!" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                    } else {
                        std::cout << std::endl;
                        std::cerr << "Invalid key '"
                                  << key
                                  << "' is passed to configuration file!" << std::endl;
                        MPI_Abort(MPI_COMM_WORLD, 143);
                    }
                }
            } else {
                std::cout << std::endl;
                std::cerr << "Invalid object detected in configuration file!" << std::endl;
                std::cerr << "Error raised when tried to read: "
                          << line
                          << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 143);
            }
        }
    }
    m_degreesOfFreedom = m_numberOfParticles * m_numberOfDimensions;
}


/* ----------------------------------------------------------------------------
  Map a wave function abbreviation to the actual wave function elements
---------------------------------------------------------------------------- */

void System::searchShortning(const std::vector<std::string> labels,
                             const std::string newLabel,
                             std::string &allLabels)
{
    if (labels.size() == 1) {
        if (allLabels == "_" + labels.at(0)) {
            allLabels = newLabel;
        }
    }
    if (labels.size() == 2) {
        for (auto &i : labels) {
            for (auto &j : labels) {
                if (i != j) {
                    if (allLabels == "_" + i + "_" + j) {
                        allLabels = newLabel;
                        break;
                    }
                }
            }
        }
    } else if (labels.size() == 3) {
        for (auto &i : labels) {
            for (auto &j : labels) {
                for (auto &k : labels) {
                    if (i != j || i != k || j != k) {
                        if (allLabels == "_" + i + "_" + j + "_" + k) {
                            allLabels = newLabel;
                            break;
                        }
                    }
                }
            }
        }
    } else if (labels.size() == 4) {
        for (auto &i : labels) {
            for (auto &j : labels) {
                for (auto &k : labels) {
                    for (auto &l : labels) {
                        if (i != j || i != k || i != l || j != k || j != l || k != l) {
                            if (allLabels == "_" + i + "_" + j + "_" + k + "_" + l) {
                                allLabels = newLabel;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
}


/* ----------------------------------------------------------------------------
  An overview of the possible trial wave function abbreviations
---------------------------------------------------------------------------- */

void System::collectAllLabels()
{
    m_trialWaveFunction = "";
    for (auto &i : m_waveFunctionElements) {
        m_trialWaveFunction += "_";
        m_trialWaveFunction += i->getLabel();
    }

    std::vector<std::string> testVMC1;
    testVMC1.push_back("Gaussian");
    searchShortning(testVMC1, "VMC", m_trialWaveFunction);

    std::vector<std::string> testVMC2;
    testVMC2.push_back("Gaussian");
    testVMC2.push_back("Padé-Jastrow");
    searchShortning(testVMC2, "VMC", m_trialWaveFunction);

    std::vector<std::string> testVMC3;
    testVMC3.push_back("Gaussian");
    testVMC3.push_back("Padé-Jastrow");
    testVMC3.push_back("Slater determinant");
    searchShortning(testVMC3, "VMC", m_trialWaveFunction);

    std::vector<std::string> testRBM1;
    testRBM1.push_back("RBM-Gaussian");
    testRBM1.push_back("RBM-product");
    searchShortning(testRBM1, "RBM", m_trialWaveFunction);

    std::vector<std::string> testRBM2;
    testRBM2.push_back("RBM-Gaussian");
    testRBM2.push_back("RBM-product");
    testRBM2.push_back("Slater determinant");
    searchShortning(testRBM2, "RBM", m_trialWaveFunction);

    std::vector<std::string> testRBMPJ1;
    testRBMPJ1.push_back("RBM-Gaussian");
    testRBMPJ1.push_back("RBM-product");
    testRBMPJ1.push_back("Padé-Jastrow");
    searchShortning(testRBMPJ1, "RBMPJ", m_trialWaveFunction);

    std::vector<std::string> testRBMPJ2;
    testRBMPJ2.push_back("RBM-Gaussian");
    testRBMPJ2.push_back("RBM-product");
    testRBMPJ2.push_back("Padé-Jastrow");
    testRBMPJ2.push_back("Slater determinant");
    searchShortning(testRBMPJ2, "RBMPJ", m_trialWaveFunction);

    std::vector<std::string> testPRBM1;
    testPRBM1.push_back("RBM-Gaussian");
    testPRBM1.push_back("RBM-product");
    testPRBM1.push_back("partlyrestricted");
    searchShortning(testPRBM1, "PRBM", m_trialWaveFunction);

    std::vector<std::string> testPRBM2;
    testPRBM2.push_back("RBM-Gaussian");
    testPRBM2.push_back("RBM-product");
    testPRBM2.push_back("partlyrestricted");
    testPRBM2.push_back("Slater determinant");
    searchShortning(testPRBM2, "PRBM", m_trialWaveFunction);

    std::vector<std::string> testDRBM1;
    testDRBM1.push_back("RBM-Gaussian");
    testDRBM1.push_back("drbmproduct");
    searchShortning(testDRBM1, "DRBM", m_trialWaveFunction);

    std::vector<std::string> testDRBM2;
    testDRBM2.push_back("RBM-Gaussian");
    testDRBM2.push_back("drbmproduct");
    testDRBM2.push_back("Slater determinant");
    searchShortning(testDRBM2, "DRBM", m_trialWaveFunction);

    std::vector<std::string> testRBMSJ1;
    testRBMSJ1.push_back("RBM-Gaussian");
    testRBMSJ1.push_back("RBM-product");
    testRBMSJ1.push_back("simple Jastrow");
    searchShortning(testRBMSJ1, "RBMSJ", m_trialWaveFunction);

    std::vector<std::string> testRBMSJ2;
    testRBMSJ2.push_back("RBM-Gaussian");
    testRBMSJ2.push_back("RBM-product");
    testRBMSJ2.push_back("simple Jastrow");
    testRBMSJ2.push_back("Slater determinant");
    searchShortning(testRBMSJ2, "RBMSJ", m_trialWaveFunction);

    std::vector<std::string> testVMC4;
    testVMC4.push_back("hydrogenlike");
    searchShortning(testVMC4, "VMC", m_trialWaveFunction);

    std::vector<std::string> testVMC5;
    testVMC5.push_back("Slater determinant");
    searchShortning(testVMC5, "VMC", m_trialWaveFunction);

    std::vector<std::string> testVMC6;
    testVMC6.push_back("hydrogenlike");
    testVMC6.push_back("Padé-Jastrow");
    searchShortning(testVMC6, "VMC", m_trialWaveFunction);

    std::vector<std::string> testVMC7;
    testVMC7.push_back("Slater determinant");
    testVMC7.push_back("Padé-Jastrow");
    searchShortning(testVMC7, "VMC", m_trialWaveFunction);

    std::vector<std::string> testVMC8;
    testVMC8.push_back("Slater determinant");
    testVMC8.push_back("Gaussian");
    searchShortning(testVMC8, "VMC", m_trialWaveFunction);

    std::vector<std::string> testSSJ1;
    testSSJ1.push_back("Slater determinant");
    testSSJ1.push_back("Gaussian");
    testSSJ1.push_back("simple Jastrow");
    searchShortning(testSSJ1, "SSJ", m_trialWaveFunction);

    std::vector<std::string> testSSJ2;
    testSSJ2.push_back("Gaussian");
    testSSJ2.push_back("simple Jastrow");
    searchShortning(testSSJ2, "SSJ", m_trialWaveFunction);

    std::vector<std::string> testFNN;
    testFNN.push_back("fnn");
    searchShortning(testFNN, "FNN", m_trialWaveFunction);

    std::vector<std::string> testBVMC;
    testBVMC.push_back("Gaussian");
    testBVMC.push_back("hard-core Jastrow");
    searchShortning(testBVMC, "BVMC", m_trialWaveFunction);
}


/* ----------------------------------------------------------------------------
  System destructor, free memory
---------------------------------------------------------------------------- */

System::~System()
{
    delete m_basis;
    delete m_initialState;
    delete m_initialWeights;
    delete m_sampler;
    delete m_optimization;
    delete m_randomNumberGenerator;
    delete m_interactionStyle;
    for (int i=m_numberOfElements; i--;)
    {
        delete m_waveFunctionElements[i];
    }
    for (unsigned int i=m_layers.size(); i--;)
    {
        delete m_layers[i];
    }
}
