#include <iostream>
#include <Eigen/Dense>

#include "system.h"


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

