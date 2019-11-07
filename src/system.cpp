#include "system.h"

void System::runSimulation()
{
    initializeSystem();
    m_lastIteration = m_numberOfIterations - m_rangeOfAdaptiveSteps - 1;
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
        if (m_iter < m_lastIteration) {
            m_globalTime += m_totalTime;
        }
        m_sampler->computeTotals();

        if (m_myRank == 0) {
            m_sampler->computeAverages();
            m_parameters -= m_optimization->updateParameters();
        }
        m_sampler->printParametersToFile();
        m_sampler->printEnergyToFile();
        if (m_iter == m_lastIteration + m_rangeOfAdaptiveSteps) {
            m_sampler->printOneBodyDensityToFile();
            m_sampler->printOneBodyDensity2ToFile();
            m_sampler->printTwoBodyDensityToFile();
        }
        printToTerminal();

        if (m_checkConvergence && m_myRank == 0) {
            checkingConvergence();
        }
        MPI_Bcast(m_parameters.data(),
                  int(m_numberOfElements * m_maxParameters),
                  MPI_DOUBLE,
                  0,
                  MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        updateAllParameters(m_parameters);
    }
}

void System::initializeSystem()
{
    initializeMPI();
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

    parser(m_configFile);

    m_sampler = new Sampler(this);
    m_sampler->openOutputFiles();
}

void System::runMetropolisCycles()
{
    for (int i = 0; i < m_stepsWEqui; i++) {
        bool acceptedStep = m_metropolis->acceptMove();
        m_positions = m_metropolis->updatePositions();
        m_distanceMatrix = m_metropolis->updateDistanceMatrix();
        m_radialVector = m_metropolis->updateRadialVector();
        if (i >= m_equilibriationSteps) {
            m_sampler->sample(acceptedStep, i);
            if (m_iter == m_lastIteration + m_rangeOfAdaptiveSteps) {
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
    if (m_iter == m_lastIteration + m_rangeOfAdaptiveSteps) {
        m_sampler->closeOutputFiles();
        if (m_myRank == 0) {
            m_sampler->doResampling();
            m_sampler->printFinalOutputToTerminal();
            std::cout << std::endl;
            std::cout << "Average CPU time: " << m_globalTime / m_lastIteration << std::endl;
            //std::cout << "Finalized successfully" << std::endl;
        }

        MPI_Finalize();
        exit(0);
    } else {
        if (m_myRank == 0) {
            m_sampler->printOutputToTerminal(m_numberOfIterations, m_totalTime);
        }
    }
}

void System::checkingConvergence()
{
    m_energies.head(m_numberOfEnergies - 1) = m_energies.tail(m_numberOfEnergies - 1);
    m_energies(m_numberOfEnergies - 1) = m_sampler->getAverageEnergy();
    if (fabs(m_energies(0) - m_energies(m_numberOfEnergies - 1)) < m_tolerance) {
        std::cout << "The system has converged! Let's run one more cycle to collect data"
                  << std::endl;
        m_lastIteration = m_iter + 1;
    }
}

int System::adaptiveSteps()
{
    int stepRatio = 1;
    if (m_iter == m_lastIteration + m_rangeOfAdaptiveSteps) {
        stepRatio = int(pow(2, m_additionalStepsLastIter));
    } else if (m_iter >= m_lastIteration) {
        stepRatio = int(pow(2, m_additionalSteps));
    }
    return stepRatio;
}

void System::setAllConstants()
{
    for (int i = 0; i < m_numberOfElements; i++) {
        m_waveFunctionElements[unsigned(i)]->setConstants(i);
    }
}

void System::initializeAllArrays(const Eigen::VectorXd positions,
                                 const Eigen::VectorXd radialVector,
                                 const Eigen::MatrixXd distanceMatrix)
{
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
    for (auto &i : m_waveFunctionElements) {
        i->setArrays();
        i->updateArrays(positions, radialVector, distanceMatrix, changedCoord);
    }
}

void System::resetAllArrays()
{
    for (auto &i : m_waveFunctionElements) {
        i->resetArrays();
    }
}

void System::updateAllParameters(const Eigen::MatrixXd parameters)
{
    for (auto &i : m_waveFunctionElements) {
        i->updateParameters(parameters);
    }
}

double System::evaluateProbabilityRatio()
{
    double ratio = 1;
    for (auto &i : m_waveFunctionElements) {
        ratio *= i->evaluateRatio();
    }
    return ratio;
}

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

Eigen::MatrixXd System::getAllParameterGradients()
{
    for (int i = 0; i < m_numberOfElements; i++) {
        m_gradients.row(i) = m_waveFunctionElements[unsigned(i)]->computeParameterGradient();
    }
    return m_gradients;
}

void System::setGlobalArraysToCalculate()
{
    // Check if the elements need distance matrix or radial distance vector
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
    // Check if the Hamiltonian needs distance matrix or radial distance vector
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

void System::setNumberOfParticles(const int numberOfParticles)
{
    assert(numberOfParticles > 0);
    m_numberOfParticles = numberOfParticles;
    initializeMPI();
}

void System::setNumberOfDimensions(const int numberOfDimensions)
{
    assert(numberOfDimensions > 1);
    assert(numberOfDimensions < 4);
    m_numberOfDimensions = numberOfDimensions;
    setNumberOfFreeDimensions();
}

void System::setNumberOfFreeDimensions()
{
    m_degreesOfFreedom = m_numberOfParticles * m_numberOfDimensions;
}

void System::setNumberOfHiddenNodes(const int numberOfHiddenNodes)
{
    assert(numberOfHiddenNodes > 0);
    m_numberOfHiddenNodes = numberOfHiddenNodes;
}

void System::setNumberOfMetropolisSteps(const int steps)
{
    // Calculate number of steps without equilibriation (power of 2)
    m_totalStepsWOEqui = steps;
    if (m_myRank == 0) {
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
    m_numberOfElements = static_cast<int>(numberOfElements);
    collectAllLabels();
}

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

void System::setMaxParameter(WaveFunction *waveFunction)
{
    int numberOfParameters = waveFunction->getNumberOfParameters();
    if (numberOfParameters > m_maxParameters) {
        m_maxParameters = numberOfParameters;
    }
    m_totalNumberOfParameters += numberOfParameters;
}

void System::setStepLength(const double stepLength)
{
    assert(stepLength > 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(const double equilibrationFraction)
{
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setFrequency(const double omega)
{
    assert(omega > 0);
    m_omega = omega;
}

void System::setScreeningTools(const bool screening,
                               const double screeningStrength,
                               const double dsl)
{
    assert(screeningStrength >= 1);
    assert(dsl > 0);
    m_screening = screening;
    m_screeningStrength = screeningStrength;
    m_dsl = dsl;
}

void System::setTotalSpin(const double totalSpin)
{
    double intpart;
    assert(std::modf(m_numberOfParticles / 2 - abs(totalSpin), &intpart) == 0.0);
    m_totalSpin = totalSpin;
}

void System::setAtomicNumber(const int Z)
{
    assert(Z > 0);
    m_Z = Z;
}

void System::setLearningRate(const double eta)
{
    assert(eta > 0);
    m_eta = eta;
}

void System::setWidth(const double sigma)
{
    assert(sigma > 0);
    m_sigma = sigma;
}

void System::setInteraction(const bool interaction)
{
    m_interaction = interaction;
}

void System::setConvergenceTools(bool checkConvergence, int numberOfEnergies, double tolerance)
{
    m_checkConvergence = checkConvergence;
    m_tolerance = tolerance;
    m_numberOfEnergies = numberOfEnergies;
    m_energies = Eigen::VectorXd::Zero(numberOfEnergies);
}

void System::setAdaptiveStepTools(bool applyAdaptiveSteps,
                                  int rangeOfAdaptiveSteps,
                                  int additionalSteps,
                                  int additionalStepsLastIteration)
{
    m_applyAdaptiveSteps = applyAdaptiveSteps;
    if (m_applyAdaptiveSteps) {
        m_rangeOfAdaptiveSteps = rangeOfAdaptiveSteps;
        m_additionalSteps = additionalSteps;
        m_additionalStepsLastIter = additionalStepsLastIteration;
    } else {
        m_rangeOfAdaptiveSteps = 0;
        m_additionalSteps = 0;
        m_additionalStepsLastIter = 0;
    }
}

void System::setDensityTools(bool computeOneBodyDensity,
                             bool computeOneBodyDensity2,
                             bool computeTwoBodyDensity,
                             int numberOfBins,
                             double maxRadius)
{
    m_computeOneBodyDensity = computeOneBodyDensity;
    m_computeOneBodyDensity2 = computeOneBodyDensity2;
    m_computeTwoBodyDensity = computeTwoBodyDensity;
    m_numberOfBins = numberOfBins;
    m_maxRadius = maxRadius;
}

void System::setEnergyPrintingTools(bool printEnergyFile, bool printInstantEnergyFile)
{
    m_printEnergyToFile = printEnergyFile;
    m_doResampling = printInstantEnergyFile;
}

void System::setParameterPrintingTools(bool printParametersToFile)
{
    m_printParametersToFile = printParametersToFile;
}

void System::initializeMPI()
{
    int initialized;
    MPI_Initialized(&initialized);
    if (!initialized)
        MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &m_numberOfProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &m_myRank);
}

void System::setPath(const std::string path)
{
    m_path = path;
}

void System::setHamiltonian(Hamiltonian *hamiltonian)
{
    m_hamiltonian = hamiltonian;
}

void System::setBasis(Basis *basis)
{
    m_basis = basis;
}

void System::setWaveFunctionElements(std::vector<class WaveFunction *> waveFunctionElements)
{
    m_waveFunctionElements = waveFunctionElements;
    //setMaxParameters();
    setNumberOfElements(waveFunctionElements.size());
    setAllConstants();
}

void System::setWaveFunctionElement(WaveFunction *waveFunction)
{
    m_waveFunctionElements.push_back(waveFunction);
    //setMaxParameter(waveFunction);
    setNumberOfElements(m_waveFunctionElements.size());
}

void System::setInitialState(InitialState *initialState)
{
    m_initialState = initialState;
}

void System::setInitialWeights(InitialWeights *initialWeights)
{
    m_initialWeights = initialWeights;
}

void System::setMetropolis(Metropolis *metropolis)
{
    m_metropolis = metropolis;
}

void System::setOptimization(Optimization *optimization)
{
    m_optimization = optimization;
}

void System::setRandomNumberGenerator(RandomNumberGenerator *randomNumberGenerator)
{
    m_randomNumberGenerator = randomNumberGenerator;
}

void System::setGradients()
{
    m_gradients = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxParameters);
}

void System::initializeFromConfig(int argc, char** argv) {
    if (argc >= 2) {
        m_configFile = argv[1];
        m_args = argc;
    }
}

void System::parser(const std::string configFile)
{
    std::ifstream infile;
    infile.open(configFile.c_str());
    if (!infile.is_open() && m_args >= 2) {
        perror("File not found");
    }
    else {
        std::string line;
        while (std::getline(infile, line)) {
            std::istringstream is_line(line);
            std::string key;
            if (std::getline(is_line, key, ':')) {
                std::string value;
                if (std::getline(is_line, value)) {
                    if (key == "numParticles") {
                        m_numberOfParticles = std::stoi(value);
                        m_numberOfHiddenNodes = m_numberOfParticles;
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
                        m_initialTotalStepsWOEqui = std::stoi(value);
                    } else if (key == "numHiddenNodes") {
                        m_numberOfHiddenNodes = std::stoi(value);
                    } else if (key == "totalSpin") {
                        m_totalSpin = std::stod(value);
                    } else if (key == "stepLength") {
                        m_stepLength = std::stod(value);
                    } else if (key == "equilibriation") {
                        m_equilibrationFraction = std::stod(value);
                    } else if (key == "interaction") {
                        m_interaction = std::stoi(value);
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
                        } else if (value == "randomize") {
                            setInitialWeights(new Randomize(this, 0.1));
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
                            std::cout << value << " is not a known wave function configuration"
                                      << std::endl;
                            MPI_Finalize();
                            exit(0);
                        }
                        setWaveFunctionElements(waveFunctionElements);
                    }
                }
            }
        }
    }
    m_degreesOfFreedom = m_numberOfParticles * m_numberOfDimensions;
}

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

void System::collectAllLabels()
{
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
}
