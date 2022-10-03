#include <iostream>
#include <fstream>
#include <cassert>
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
            m_stepsWEqui = m_stepsWOEqui + m_burnInSteps;
            m_totalStepsWOEqui = m_initialTotalStepsWOEqui * adaptiveSteps();
            m_totalStepsWEqui = m_totalStepsWOEqui + m_totalBurnInSteps;
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
        if (m_iter % m_checkpointFreq == 0) { // checkpointing
            m_sampler->printParametersToFile("weights_" + std::to_string(m_iter) + ".dat");
        }
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
        if (i >= m_burnInSteps) {
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
