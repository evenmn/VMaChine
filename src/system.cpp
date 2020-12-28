#include "system.h"

void System::printLogo()
{
    /* Print system information to terminal */
    std::cout << "__      ____  __          _____ _    _ _____ _   _ ______ " << std::endl;
    std::cout << "\\ \\    / /  \\/  |   /\\   / ____| |  | |_   _| \\ | |  ____|" << std::endl;
    std::cout << " \\ \\  / /| \\  / |  /  \\ | |    | |__| | | | |  \\| | |__   " << std::endl;
    std::cout << "  \\ \\/ / | |\\/| | / /\\ \\| |    |  __  | | | | . ` |  __|  " << std::endl;
    std::cout << "   \\  /  | |  | |/ ____ \\ |____| |  | |_| |_| |\\  | |____ " << std::endl;
    std::cout << "    \\/   |_|  |_/_/    \\_\\_____|_|  |_|_____|_| \\_|______|" << std::endl;
    std::cout << "==========================================================" << std::endl;
    std::cout << std::endl;
}

void System::initializeMPI()
{
    /* Initialize MPI based on command line arguments.
     * Command line arguments are taken automatically
     * from main. */
    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized)
        MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &m_numberOfProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
}

void System::initializeSystem()
{
    /* Initialize system
     * The various objects need to be initialized in the correct order,
     * which is the task of this function. This lets the calls be in arbitrary
     * order in main/configuration file */
    initializeMPI();
    parser(m_configFile);
    setInputLayer(m_degreesOfFreedom);      // Add input layer
    setOutputLayer(new Sigmoid(this));    // Add output layer
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
    //parser(m_configFile);

    m_sampler = new Sampler(this);
    m_sampler->openOutputFiles();
}

void System::printSystemInformation()
{
    /* Print system information to terminal */
    std::cout << std::fixed;
    std::cout << std::setprecision(6);
    std::cout << std::endl;
    std::cout << "#############################################" << std::endl;
    std::cout << "### =====    SYSTEM INFORMATION     ===== ###" << std::endl;
    std::cout << "#############################################" << std::endl;
    std::cout << "Number of particles:      " << m_numberOfParticles << std::endl;
    std::cout << "Number of dimensions:     " << m_numberOfDimensions << std::endl;
    std::cout << "Interaction:              " << m_interaction << std::endl;
    std::cout << "Hamiltonian:              " << m_hamiltonian << std::endl;
    std::cout << "Initial state:            " << m_initialState << std::endl;
    std::cout << "Oscillator frequency:     " << m_omega << std::endl;
    std::cout << std::endl;
    std::cout << "#############################################" << std::endl;
    std::cout << "### ===== WAVE FUNCTION INFORMATION ===== ###" << std::endl;
    std::cout << "#############################################" << std::endl;
    for (int i = 0; i < m_numberOfElements; i++) {
        std::cout << "Element " << i << ":                " << m_waveFunctionElements[unsigned(i)] << std::endl;
    }
    std::cout << "Basis:                    " << m_basis << std::endl;
    std::cout << "Initial parameters:       " << m_initialWeights << std::endl;
    std::cout << "Number of parameters:     " << this->getTotalNumberOfParameters() << std::endl;
    std::cout << std::endl;
    std::cout << "#############################################" << std::endl;
    std::cout << "### =====  SIMULATION INFORMATION   ===== ###" << std::endl;
    std::cout << "#############################################" << std::endl;
    std::cout << "Max number of iterations: " << m_numberOfIterations << std::endl;
    std::cout << "Number of MPI threads:    " << m_numberOfProcesses << std::endl;
    std::cout << " # Monte Carlo Cycles    : " << m_totalStepsWEqui << " (" << m_totalStepsWOEqui << " after burn-in)" << std::endl;
    std::cout << std::endl;

}

void System::printHeaderLine()
{
    std::cout << "Step  " << "Energy  " << "Energy_STD  " << "Kinetic  "
              << "Kinetic_STD  " << "External  " << "External_STD  "
              << "Interaction  " << "Interaction_STD  " << "Acceptence  "
              << "CPU_time" << std::endl;
}

void System::runSimulation()
{
    /* Run simulation specified in main/configuration file.
     * This is the main loop in VMC, where we update the parameters. */
    printLogo();
    auto start = std::chrono::system_clock::now();
    std::time_t start_time = std::chrono::system_clock::to_time_t(start);
    std::cout << "Started computation at " << std::ctime(&start_time) << std::endl;
    initializeSystem();
    printSystemInformation();
    printHeaderLine();
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
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "Finished computation at " << std::ctime(&end_time)
              << "Elapsed time: " << elapsed_seconds.count() << "s\n";
}

void System::runMetropolisCycles()
{
    /* This is the sampling loop, where we move the particles
     * and calculate expectation values */
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

void System::printToTerminal()
{
    /* Call this function to print to terminal */
    if (m_iter == m_numberOfNormalIterations + m_rangeOfAdaptiveSteps) {
        m_sampler->closeOutputFiles();
        if (m_rank == 0) {
            m_sampler->doResampling();
            m_sampler->printFinalOutputToTerminal();
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

void System::checkingConvergence()
{
    /* This functions checks if the trial wave function has converged */
    m_energies.head(m_numberOfEnergies - 1) = m_energies.tail(m_numberOfEnergies - 1);
    m_energies(m_numberOfEnergies - 1) = m_sampler->getAverageEnergy();
    if (fabs(m_energies(0) - m_energies(m_numberOfEnergies - 1)) < m_tolerance) {
        std::cout << "The system has converged! Let's run one more cycle to collect data"
                  << std::endl;
        m_numberOfNormalIterations = m_iter + 1;
        m_checkConvergence = false;
    }
}

int System::adaptiveSteps()
{
    /* Is responsible for the adaptive steps */
    int stepRatio = 1;
    if (m_iter == m_numberOfNormalIterations + m_rangeOfAdaptiveSteps) {
        stepRatio = int(pow(2, m_additionalStepsLastIter));

    } else if (m_iter >= m_numberOfNormalIterations) {
        stepRatio = int(pow(2, m_additionalSteps));
    }
    return stepRatio;
}


// OPERATIONS ON WAVE FUNCTION ELEMENTS

void System::setAllConstants()
{
    /* Initialize the wave function elements with essential variables
     * specified by used. */
    for (int i = 0; i < m_numberOfElements; i++) {
        m_waveFunctionElements[unsigned(i)]->setConstants(i);
    }
}

void System::initializeAllArrays(const Eigen::VectorXd positions,
                                 const Eigen::VectorXd radialVector,
                                 const Eigen::MatrixXd distanceMatrix)
{
    /* Initialize the wave function elements with positions. */
    for (auto &i : m_waveFunctionElements) {
        i->initializeArrays(positions, radialVector, distanceMatrix);
        i->setArrays();
    }
}

void System::updateAllArrays(const Eigen::VectorXd positions,
                             const Eigen::VectorXd radialVector,
                             const Eigen::MatrixXd distanceMatrix,
                             const int changedCoord)
{
    /* Update positions in all wave function elements when a
     * particle is moved. */
    for (auto &i : m_waveFunctionElements) {
        i->setArrays();
        i->updateArrays(positions, radialVector, distanceMatrix, changedCoord);
    }
}

void System::resetAllArrays()
{
    /* Reset positions in all ave function elements when a
     * particle move is rejected. */
    for (auto &i : m_waveFunctionElements) {
        i->resetArrays();
    }
}

void System::updateAllParameters(const Eigen::MatrixXd parameters)
{
    /* Update the parameters/weights in all the wave function
     * elements when the parameters are updated. */
    for (auto &i : m_waveFunctionElements) {
        i->updateParameters(parameters);
    }
}

double System::evaluateProbabilityRatio()
{
    /* Evaluate the collective probability ratio. */
    double ratio = 1;
    for (auto &i : m_waveFunctionElements) {
        ratio *= i->evaluateRatio();
    }
    return ratio;
}

double System::getKineticEnergy()
{
    /* Obtain the total kinetic energy of the system
     * based on the gradients and Laplacians of all
     * elements. */
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

Eigen::MatrixXd System::getAllParameterGradients()
{
    /* Get the gradient of all the wave function elements with respect to
     * all the parameters. To be used in the parameter update. */
    for (int i = 0; i < m_numberOfElements; i++) {
        m_gradients.row(i) = m_waveFunctionElements[unsigned(i)]->computeParameterGradient();
    }
    return m_gradients;
}

void System::setGlobalArraysToCalculate()
{
    /* Check if the elements need distance matrix and/or radial distance vector */
    for (auto &p : m_waveFunctionElements) {
        int need = p->getGlobalArrayNeed();
        if (need == 1) {
            m_calculateDistanceMatrix = true;
        }
        if (need == 2) {
            m_calculateRadialVector = true;
        }
        if (need == 3) {
            m_calculateDistanceMatrix = true;
            m_calculateRadialVector = true;
        }
    }
    /* Check if the Hamiltonian needs distance matrix and/or radial distance vector */
    int need = m_hamiltonian->getGlobalArrayNeed();
    if (need == 1) {
        m_calculateDistanceMatrix = true;
    }
    if (need == 2) {
        m_calculateRadialVector = true;
    }
    if (need == 3) {
        m_calculateDistanceMatrix = true;
        m_calculateRadialVector = true;
    }
}

void System::setMaxParameters()
{
    /* Set the maximum number of parameters found in a
     * wave function element. This will be used to create
     * the parameter matrix, with dim
     * (maxNumberOfParameter x numberOfElements). Determines
     * also the total number of particles. Called
     * automatically when needed. */
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

void System::setNumberOfParticles(const int numberOfParticles)
{
    /* Specify the number of particles, N, to be used. */
    assert(numberOfParticles > 0);
    m_numberOfParticles = numberOfParticles;
    initializeMPI();
}

void System::setNumberOfDimensions(const int numberOfDimensions)
{
    /* Specify the number of dimensions, d, of the system. Supported 2 and 3 dimensions. */
    assert(numberOfDimensions > 1);
    assert(numberOfDimensions < 4);
    m_numberOfDimensions = numberOfDimensions;
    setNumberOfFreeDimensions();
}

void System::setNumberOfFreeDimensions()
{
    /* Set the number of free dimensions, F=Nd. This is called automatically after the
     * number of particles and number of dimensions are specified. */
    m_degreesOfFreedom = m_numberOfParticles * m_numberOfDimensions;
}

void System::setNumberOfHiddenUnits(const int numberOfHiddenUnits)
{
    /* Set the number of hidden units used in RBM trial wave function. */
    assert(numberOfHiddenUnits > 0);
    m_numberOfHiddenUnits = numberOfHiddenUnits;
}

void System::setNumberOfMetropolisCycles(const int steps)
{
    /* Determine the number of steps used by each process
     * when the equilibriation (burn-in period) is excluded
     * (a power of 2) and when it is included. */
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

void System::setNumberOfElements(const unsigned long numberOfElements)
{
    /* Set the number of wave function elements. This is called
     * automatically when needed. */
    m_numberOfElements = static_cast<int>(numberOfElements);
    collectAllLabels();
}

void System::setNumberOfIterations(const int numberOfIterations)
{
    /* Set maximum number of iterations used in the Monte Carlo
     * integration. The actual number of iterations used can
     * be overruled by the stop criterion. */
    assert(numberOfIterations > 0);
    m_numberOfIterations = numberOfIterations;
}

void System::setStepLength(const double stepLength)
{
    /* Set step length used in sampling. */
    assert(stepLength > 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(const double equilibrationFraction)
{
    /* Set equilibriation fraction (burn-in period). */
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
    setNumberOfMetropolisCycles(m_totalStepsWOEqui);
}

void System::setFrequency(const double omega)
{
    /* Set frequency of harmonic oscillator potential. */
    assert(omega > 0);
    m_omega = omega;
}

void System::setScreeningTools(const double screeningStrength,
                               const double dsl)
{
    /* By calling this function, the screening is enabled.
     * The screening strength and Debye screening length (DSL)
     * need to be specified. */
    assert(screeningStrength >= 1);
    assert(dsl > 0);
    m_screening = true;
    m_screeningStrength = screeningStrength;
    m_dsl = dsl;
}

void System::setAtomicNumber(const int Z)
{
    /* Set atomic number, Z, when simulating atoms. */
    assert(Z > 0);
    m_Z = Z;
}

void System::setLearningRate(const double eta)
{
    /* Set learning rate of simulation. */
    assert(eta > 0);
    m_eta = eta;
}

void System::setWidth(const double sigma)
{
    /* Set distribution width of Gaussian distribution
     * when using a Gaussian-binary restricted Boltzmann
     * machine as the trial wave fucntion guess. */
    assert(sigma > 0);
    m_sigma = sigma;
}

void System::setInteraction(const bool interaction)
{
    /* Decided whether or not the electrons should be
     * interacting. */
    m_interaction = interaction;
}

void System::setConvergenceTools(int numberOfEnergies, double tolerance)
{
    /* Specify convergence. Convergence is enabled when
     * this function is called. */
    m_checkConvergence = true;
    m_tolerance = tolerance;
    m_numberOfEnergies = numberOfEnergies;
    m_energies = Eigen::VectorXd::Zero(numberOfEnergies);
}

void System::setAdaptiveStepTools(int rangeOfAdaptiveSteps,
                                  int additionalSteps,
                                  int additionalStepsLastIteration)
{
    /* Set adaptive steps. Adaptive steps are enabled when
     * this function is called. Default is 4 times as many
     * cycles when 10 iterations are remaining and 8 times
     * as many cycles for the last iteration. */
    m_applyAdaptiveSteps = true;
    m_rangeOfAdaptiveSteps = rangeOfAdaptiveSteps;
    m_additionalSteps = additionalSteps;
    m_additionalStepsLastIter = additionalStepsLastIteration;
}

void System::computeRadialOneBodyDensity(int numberOfBins, double maxRadius)
{
    /* Compute the radial one-body density. Default number of bins is
     * 1000, default max radius is 50 given in natural units. */
    m_computeOneBodyDensity = true;
    m_numberOfBins = numberOfBins;
    m_maxRadius = maxRadius;
}

void System::computeSpatialOneBodyDensity(int numberOfBins, double maxRadius)
{
    /* Compute the spatial one-body density. Default number of bins is
     * 1000 x 1000, default max radius is 50 given in natural units. */
    m_computeOneBodyDensity2 = true;
    m_numberOfBins = numberOfBins;
    m_maxRadius = maxRadius;
}

void System::computeTwoBodyDensity(int numberOfBins, double maxRadius)
{
    /* Compute the radial two-body density. Default number of bins is
     * 1000, default max radius is 50 given in natural units. */
    m_computeTwoBodyDensity = true;
    m_numberOfBins = numberOfBins;
    m_maxRadius = maxRadius;
}

void System::setTotalSpin(const double totalSpin)
{
    /* Set total spin of the system. To be implemented. */
    double intpart;
    assert(std::modf(m_numberOfParticles / 2 - abs(totalSpin), &intpart) == 0.0);
    m_totalSpin = totalSpin;
}

void::System::dumpEnergyToFile(bool printEnergyFile) {
    /* Dump the energy expectation value after each iteration
     * to file. Default setting is to print energy file. */
    m_printEnergyToFile = printEnergyFile;
}

void::System::doBlocking(bool printInstantEnergyFile) {
    /* Dump the instant energy to file (energies from all
     * the cycles). This will be used in the resampling.
     * Default setting is to print instant file for the
     * last iteration. */
    m_doResampling = printInstantEnergyFile;
}

void System::dumpParametersToFile(bool printParametersToFile)
{
    /* Print parameters to file after each iterations. Makes it
     * possible to start where the previous simulation ended.
     * Default setting is to print parameters to file. */
    m_printParametersToFile = printParametersToFile;
}

void System::setPath(const std::string path)
{
    /* Set path to where to save the files. */
    m_path = path;
}


// CLASSES SPECIFIED IN MAIN

void System::setHamiltonian(Hamiltonian *hamiltonian)
{
    /* Specify which Hamiltonian to use. */
    m_hamiltonian = hamiltonian;
}

void System::setBasis(Basis *basis)
{
    /* Specify the basis to be used in the Slater determinant. */
    m_basis = basis;
}

void System::setWaveFunctionElements(std::vector<class WaveFunction *> waveFunctionElements)
{
    /* Set all wave function elements, based on a std::vector. */
    m_waveFunctionElements = waveFunctionElements;
    setNumberOfElements(waveFunctionElements.size());
    setAllConstants();
}

void System::setWaveFunctionElement(WaveFunction *waveFunction)
{
    /* Add a wave function element. */
    m_waveFunctionElements.push_back(waveFunction);
    setNumberOfElements(m_waveFunctionElements.size());
}

void System::setInputLayer(int numberOfUnits)
{
    /* Set the input layer when a feed-forward neural network
     * (FNN) is used as a trial wave function element. */
    m_layers.insert(m_layers.begin(), new Input(this, numberOfUnits));
}

void System::addDenseLayer(int numberOfUnits, Activation *activation)
{
    /* Add dense layer to a feed-forward neural network. */
    m_layers.push_back(new Dense(this, numberOfUnits, activation));
    setNumberOfElements(m_waveFunctionElements.size());
}

void System::setOutputLayer(Activation *activation)
{
    /* Add dense layer to a feed-forward neural network. */
    m_layers.push_back(new Output(this, activation));
    setNumberOfElements(m_waveFunctionElements.size());
}

void System::setInitialState(InitialState *initialState)
{
    /* Initialize particle positions.
     * Possible choices are:
     *  - Random normal (default)
     *  - Random uniform */
    m_initialState = initialState;
}

void System::setInitialWeights(InitialWeights *initialWeights)
{
    /* Initialize parameters.
     * Possible choices are:
     *  - Automatize (default)
     *  - Constant
     *  - Random uniform */
    m_initialWeights = initialWeights;
}

void System::setMetropolis(Metropolis *metropolis)
{
    /* Specify sampling algorithm.
     * Possible choices:
     *  - Brute force
     *  - Importance sampling */
    m_metropolis = metropolis;
}

void System::setOptimization(Optimization *optimization)
{
    /* Specify optimization tool when updating parameters.
     * Possible choices:
     *  - Gradient descent
     *  - Stochastic gradient descent
     *  - Adam */
    m_optimization = optimization;
}

void System::setRandomNumberGenerator(RandomNumberGenerator *randomNumberGenerator)
{
    /* Specify which random number generator to be used.
     * Currently, only Mersenne Twister 19337 is available. */
    m_randomNumberGenerator = randomNumberGenerator;
}

void System::setGradients()
{
    /* Initialize gradient matrix. */
    m_gradients = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxParameters);
}

void System::initializeFromConfig(int argc, char** argv) {
    /* Initialize variables from configuration file. Will
     * overwrite main and default settings. */
    if (argc >= 2) {
        m_configFile = argv[1];
        m_args = argc;
    }
}

void System::parser(const std::string configFile)
{
    /* Parse from configuration file. */
    std::ifstream infile;
    infile.open(configFile.c_str());
    if (!infile.is_open() && m_args >= 2) {
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
                    value = trim(value);
                    if (key == "numParticles") {
                        m_numberOfParticles = std::stoi(value);
                        m_numberOfHiddenUnits = m_numberOfParticles;
                        m_Z = m_numberOfParticles;
                    } else if (key == "numDimensions") {
                        m_numberOfDimensions = std::stoi(value);
                    } else if (key == "omega") {
                        m_omega = std::stod(value);
                        m_stepLength = 0.1 / sqrt(m_omega);
                        m_sigma = 1.0 / sqrt(m_omega);
                    } else if (key == "learningRate") {
                        m_eta = std::stod(value);
                    } else if (key == "maxRadius") {
                        m_maxRadius = std::stod(value);
                    } else if (key == "numIterations") {
                        m_numberOfIterations = std::stoi(value);
                    } else if (key == "numSteps") {
                        setNumberOfMetropolisCycles(std::stoi(value));
                    } else if (key == "numHiddenNodes") {
                        m_numberOfHiddenUnits = std::stoi(value);
                    } else if (key == "totalSpin") {
                        m_totalSpin = std::stod(value);
                    } else if (key == "stepLength") {
                        m_stepLength = std::stod(value);
                    } else if (key == "equilibration") {
                        m_equilibrationFraction = std::stod(value);
                    } else if (key == "interaction") {
                        m_interaction = std::stoi(value);
                        // istringstream(value) >> std::boolalpha >> m_interaction;
                    } else if (key == "checkConvergence") {
                        m_checkConvergence = std::stoi(value);
                    } else if (key == "applyAdaptiveSteps") {
                        m_applyAdaptiveSteps = std::stoi(value);
                    } else if (key == "computeOneBodyDensity") {
                        m_computeOneBodyDensity = std::stoi(value);
                    } else if (key == "computeTwoBodyDensity") {
                        m_computeTwoBodyDensity = std::stoi(value);
                    } else if (key == "printEnergyToFile") {
                        m_printEnergyToFile = std::stoi(value);
                    } else if (key == "printParametersToFile") {
                        m_printParametersToFile = std::stoi(value);
                    } else if (key == "doResampling") {
                        m_doResampling = std::stoi(value);
                    } else if (key == "path") {
                        m_path = value;
                    } else if (key == "numberOfEnergies") {
                        m_numberOfEnergies = std::stoi(value);
                    } else if (key == "tolerance") {
                        m_tolerance = std::stod(value);
                    } else if (key == "rangeOfAdaptiveSteps") {
                        m_rangeOfAdaptiveSteps = std::stoi(value);
                    } else if (key == "additionalSteps") {
                        m_additionalSteps = std::stoi(value);
                    } else if (key == "additionalStepsLastIter") {
                        m_additionalStepsLastIter = std::stoi(value);
                    } else if (key == "numberOfBins") {
                        m_numberOfBins = std::stoi(value);
                    } else if (key == "basis") {
                        if (value == "hermite") {
                            setBasis(new Hermite(this));
                        } else if (value == "hermiteExpansion") {
                            setBasis(new HermiteExpansion(this));
                        } else {
                            std::cout << value << " is not a known basis" << std::endl;
                            MPI_Finalize();
                            exit(0);
                        }
                    } else if (key == "hamiltonian") {
                        if (value == "harmonicOscillator") {
                            setHamiltonian(new HarmonicOscillator(this));
                        } else if (value == "doubleWell") {
                            setHamiltonian(new DoubleWell(this, 2));
                        } else {
                            std::cout << value << " is not a known Hamiltonian" << std::endl;
                            MPI_Finalize();
                            exit(0);
                        }
                    } else if (key == "optimization") {
                        if (value == "adam") {
                            setOptimization(new ADAM(this));
                        } else if (value == "gd") {
                            setOptimization(new GradientDescent(this, 0.0, 0.0));
                        } else if (value == "sgd") {
                            setOptimization(new SGD(this, 0.0, 0.0));
                        } else {
                            std::cout << value << " is not a known optimization tool" << std::endl;
                            MPI_Finalize();
                            exit(0);
                        }
                    } else if (key == "initialWeights") {
                        if (value == "automatize") {
                            setInitialWeights(new Automatize(this));
                        } else if (value == "randomuniform") {
                            setInitialWeights(new RandomUniformWeights(this, 0.1));
                        } else if (value == "constant") {
                            setInitialWeights(new Constant(this, 1.0));
                        } else {
                            std::cout << value << " is not a known initial weight configuration"
                                      << std::endl;
                            MPI_Finalize();
                            exit(0);
                        }
                    } else if (key == "initialState") {
                        if (value == "randomNormal") {
                            setInitialState(new RandomNormal(this));
                        } else if (value == "randomUniform") {
                            setInitialState(new RandomUniform(this));
                        } else {
                            std::cout << value << " is not a known initial state configuration"
                                      << std::endl;
                            MPI_Finalize();
                            exit(0);
                        }
                    } else if (key == "sampling") {
                        if (value == "importanceSampling") {
                            setMetropolis(new ImportanceSampling(this));
                        } else if (value == "bruteForce") {
                            setMetropolis(new BruteForce(this));
                        } else {
                            std::cout << value << " is not a known sampling tool" << std::endl;
                            MPI_Finalize();
                            exit(0);
                        }
                    } else if (key == "waveFunction") {
                        std::vector<class WaveFunction *> waveFunctionElements;
                        if (value == "VMC") {
                            waveFunctionElements.push_back(new class Gaussian(this));
                            waveFunctionElements.push_back(new class SlaterDeterminant(this));
                            waveFunctionElements.push_back(new class PadeJastrow(this));
                        } else if (value == "RBM") {
                            waveFunctionElements.push_back(new class SlaterDeterminant(this));
                            waveFunctionElements.push_back(new class RBMGaussian(this));
                            waveFunctionElements.push_back(new class RBMProduct(this));
                        } else if (value == "RBMPJ") {
                            waveFunctionElements.push_back(new class SlaterDeterminant(this));
                            waveFunctionElements.push_back(new class RBMGaussian(this));
                            waveFunctionElements.push_back(new class RBMProduct(this));
                            waveFunctionElements.push_back(new class PadeJastrow(this));
                        } else if (value == "RBMSJ") {
                            waveFunctionElements.push_back(new class SlaterDeterminant(this));
                            waveFunctionElements.push_back(new class RBMGaussian(this));
                            waveFunctionElements.push_back(new class RBMProduct(this));
                            waveFunctionElements.push_back(new class SimpleJastrow(this));
                        } else if (value == "PRBM") {
                            waveFunctionElements.push_back(new class SlaterDeterminant(this));
                            waveFunctionElements.push_back(new class RBMGaussian(this));
                            waveFunctionElements.push_back(new class RBMProduct(this));
                            waveFunctionElements.push_back(new class PartlyRestricted(this));
                        } else {
                            std::cout << "Error: " << value << " is not a known wave function configuration"
                                      << std::endl;
                            MPI_Finalize();
                            exit(0);
                        }
                        setWaveFunctionElements(waveFunctionElements);
                    } else if (key == "waveFunctionElement") {
                        if (value == "gaussian") {
                            setWaveFunctionElement(new class Gaussian(this));
                        } else if (value == "slaterDeterminant") {
                            setWaveFunctionElement(new class SlaterDeterminant(this));
                        } else if (value == "padeJastrow") {
                            setWaveFunctionElement(new class PadeJastrow(this));
                        } else if (value == "simpleJastrow") {
                            setWaveFunctionElement(new class SimpleJastrow(this));
                        } else if (value == "RBMGaussian") {
                            setWaveFunctionElement(new class RBMGaussian(this));
                        } else if (value == "RBMProduct") {
                            setWaveFunctionElement(new class RBMProduct(this));
                        } else {
                            std::cerr << "Wave function element does not exist" << std::endl;
                            MPI_Abort(MPI_COMM_WORLD, 143);
                        }
                    } else {
                      std::cerr << "Invalid key is passed to configuration file" << std::endl;
                      MPI_Abort(MPI_COMM_WORLD, 143);
                    }
                }
            } else {
                std::cerr << "Invalid key is passed to configuration file" << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 143);
            }
        }
    }
    m_degreesOfFreedom = m_numberOfParticles * m_numberOfDimensions;
}

void System::searchShortning(const std::vector<std::string> labels,
                             const std::string newLabel,
                             std::string &allLabels)
{
    /* Search if a shortning of the specified wave function element
     * configuration exist, and eventually label it as the shortning. */
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

void System::collectAllLabels()
{
    /* Possible shortnings. */
    m_trialWaveFunction = "";
    for (auto &i : m_waveFunctionElements) {
        m_trialWaveFunction += "_";
        m_trialWaveFunction += i->getLabel();
    }

    std::vector<std::string> testVMC1;
    testVMC1.push_back("gaussian");
    searchShortning(testVMC1, "VMC", m_trialWaveFunction);

    std::vector<std::string> testVMC2;
    testVMC2.push_back("gaussian");
    testVMC2.push_back("padejastrow");
    searchShortning(testVMC2, "VMC", m_trialWaveFunction);

    std::vector<std::string> testVMC3;
    testVMC3.push_back("gaussian");
    testVMC3.push_back("padejastrow");
    testVMC3.push_back("slaterdeterminant");
    searchShortning(testVMC3, "VMC", m_trialWaveFunction);

    std::vector<std::string> testRBM1;
    testRBM1.push_back("rbmgaussian");
    testRBM1.push_back("rbmproduct");
    searchShortning(testRBM1, "RBM", m_trialWaveFunction);

    std::vector<std::string> testRBM2;
    testRBM2.push_back("rbmgaussian");
    testRBM2.push_back("rbmproduct");
    testRBM2.push_back("slaterdeterminant");
    searchShortning(testRBM2, "RBM", m_trialWaveFunction);

    std::vector<std::string> testRBMPJ1;
    testRBMPJ1.push_back("rbmgaussian");
    testRBMPJ1.push_back("rbmproduct");
    testRBMPJ1.push_back("padejastrow");
    searchShortning(testRBMPJ1, "RBMPJ", m_trialWaveFunction);

    std::vector<std::string> testRBMPJ2;
    testRBMPJ2.push_back("rbmgaussian");
    testRBMPJ2.push_back("rbmproduct");
    testRBMPJ2.push_back("padejastrow");
    testRBMPJ2.push_back("slaterdeterminant");
    searchShortning(testRBMPJ2, "RBMPJ", m_trialWaveFunction);

    std::vector<std::string> testPRBM1;
    testPRBM1.push_back("rbmgaussian");
    testPRBM1.push_back("rbmproduct");
    testPRBM1.push_back("partlyrestricted");
    searchShortning(testPRBM1, "PRBM", m_trialWaveFunction);

    std::vector<std::string> testPRBM2;
    testPRBM2.push_back("rbmgaussian");
    testPRBM2.push_back("rbmproduct");
    testPRBM2.push_back("partlyrestricted");
    testPRBM2.push_back("slaterdeterminant");
    searchShortning(testPRBM2, "PRBM", m_trialWaveFunction);

    std::vector<std::string> testDRBM1;
    testDRBM1.push_back("rbmgaussian");
    testDRBM1.push_back("drbmproduct");
    searchShortning(testDRBM1, "DRBM", m_trialWaveFunction);

    std::vector<std::string> testDRBM2;
    testDRBM2.push_back("rbmgaussian");
    testDRBM2.push_back("drbmproduct");
    testDRBM2.push_back("slaterdeterminant");
    searchShortning(testDRBM2, "DRBM", m_trialWaveFunction);

    std::vector<std::string> testRBMSJ1;
    testRBMSJ1.push_back("rbmgaussian");
    testRBMSJ1.push_back("rbmproduct");
    testRBMSJ1.push_back("simplejastrow");
    searchShortning(testRBMSJ1, "RBMSJ", m_trialWaveFunction);

    std::vector<std::string> testRBMSJ2;
    testRBMSJ2.push_back("rbmgaussian");
    testRBMSJ2.push_back("rbmproduct");
    testRBMSJ2.push_back("simplejastrow");
    testRBMSJ2.push_back("slaterdeterminant");
    searchShortning(testRBMSJ2, "RBMSJ", m_trialWaveFunction);

    std::vector<std::string> testVMC4;
    testVMC4.push_back("hydrogenlike");
    searchShortning(testVMC4, "VMC", m_trialWaveFunction);

    std::vector<std::string> testVMC5;
    testVMC5.push_back("slaterdeterminant");
    searchShortning(testVMC5, "VMC", m_trialWaveFunction);

    std::vector<std::string> testVMC6;
    testVMC6.push_back("hydrogenlike");
    testVMC6.push_back("padejastrow");
    searchShortning(testVMC6, "VMC", m_trialWaveFunction);

    std::vector<std::string> testVMC7;
    testVMC7.push_back("slaterdeterminant");
    testVMC7.push_back("padejastrow");
    searchShortning(testVMC7, "VMC", m_trialWaveFunction);

    std::vector<std::string> testVMC8;
    testVMC8.push_back("slaterdeterminant");
    testVMC8.push_back("gaussian");
    searchShortning(testVMC8, "VMC", m_trialWaveFunction);

    std::vector<std::string> testSSJ1;
    testSSJ1.push_back("slaterdeterminant");
    testSSJ1.push_back("gaussian");
    testSSJ1.push_back("simplejastrow");
    searchShortning(testSSJ1, "SSJ", m_trialWaveFunction);

    std::vector<std::string> testSSJ2;
    testSSJ2.push_back("gaussian");
    testSSJ2.push_back("simplejastrow");
    searchShortning(testSSJ2, "SSJ", m_trialWaveFunction);

    std::vector<std::string> testFNN;
    testFNN.push_back("fnn");
    searchShortning(testFNN, "FNN", m_trialWaveFunction);
}
