#include "system.h"
#include "sampler.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "Basis/basis.h"
#include "InitialStates/initialstate.h"
#include "InitialWeights/initialweights.h"
#include "Metropolis/metropolis.h"
#include "Optimization/optimization.h"
#include "RNG/rng.h"
//#include "tqdm/tqdm.h"

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <ctime>
#include <string>

void System::runIterations(const unsigned int numberOfIterations) {
    m_positions                 = m_initialState->getParticles();
    m_distanceMatrix            = m_initialState->getDistanceMatrix();
    m_radialVector              = m_initialState->getRadialVector();
    m_parameters                = m_initialWeights->getWeights();
    m_sampler                   = new Sampler(this);
    m_sampler->openOutputFiles();
    m_lastIteration = numberOfIterations - m_rangeOfDynamicSteps - 1;

    for(m_iter = 0; m_iter < numberOfIterations; m_iter++) {
    //for(m_iter : tqdm::range(numberOfIterations)) {
        if(m_applyAdaptiveSteps) {
            m_numberOfStepsWOEqui      = m_initialNumberOfStepsWOEqui * adaptiveSteps();
            m_numberOfStepsWEqui       = m_numberOfStepsWOEqui + m_numberOfEquilibriationSteps;
            m_totalNumberOfStepsWOEqui = m_initialTotalNumberOfStepsWOEqui * adaptiveSteps();
            m_totalNumberOfStepsWEqui  = m_totalNumberOfStepsWOEqui + m_totalNumberOfEquilibriationSteps;
        }
        m_sampler->setNumberOfSteps(m_numberOfStepsWOEqui, m_totalNumberOfStepsWOEqui, m_totalNumberOfStepsWEqui);
        double startTime = MPI_Wtime();
        runMetropolisCycles();
        double endTime = MPI_Wtime();
        double time = endTime - startTime;

        MPI_Reduce(&time, &m_totalTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        m_sampler->computeTotals();

        if(m_myRank == 0) {
            m_sampler->computeAverages();
            m_parameters -= m_optimization->updateParameters();
            std::cout << m_optimization->updateParameters() << std::endl;
        }
        m_sampler->printEnergyToFile();
        m_sampler->printParametersToFile();
        if(m_iter == m_lastIteration + m_rangeOfDynamicSteps) {
            m_sampler->printOneBodyDensityToFile();
            m_sampler->printTwoBodyDensityToFile();
        }
        printToTerminal(numberOfIterations);

        if(m_checkConvergence && m_myRank == 0) {
            checkingConvergence();
        }

        MPI_Bcast(m_parameters.data(), int(m_numberOfWaveFunctionElements*m_maxNumberOfParametersPerElement), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        updateAllParameters(m_parameters);
    }
}

void System::runMetropolisCycles() {
    for(unsigned long i=0; i < m_numberOfStepsWEqui; i++) {
        bool acceptedStep = m_metropolis->acceptMove();
        m_positions       = m_metropolis->updatePositions();
        m_distanceMatrix  = m_metropolis->updateDistanceMatrix();
        m_radialVector    = m_metropolis->updateRadialVector();
        if(i >= m_numberOfEquilibriationSteps) {
            m_sampler->sample(acceptedStep, i);
            if(m_iter == m_lastIteration + m_rangeOfDynamicSteps) {
                m_sampler->printInstantValuesToFile();
                m_sampler->computeOneBodyDensity(m_positions);
                m_sampler->computeTwoBodyDensity(m_positions);
            }
        }
    }
}

void System::printToTerminal(unsigned int numberOfIterations) {
    if(m_iter == m_lastIteration + m_rangeOfDynamicSteps) {
        m_sampler->closeOutputFiles();
        if(m_myRank == 0) {
            m_sampler->doResampling();
            m_sampler->printFinalOutputToTerminal();
        }
        MPI_Finalize();
        exit(0);
    }
    else {
        if(m_myRank == 0) {
            m_sampler->printOutputToTerminal(numberOfIterations, m_totalTime);
        }
    }
}

void System::checkingConvergence() {
    m_energies.head(m_numberOfEnergies-1) = m_energies.tail(m_numberOfEnergies-1);
    m_energies(m_numberOfEnergies-1) = m_sampler->getAverageEnergy();
    if(fabs(m_energies(0) - m_energies(m_numberOfEnergies-1)) < m_tolerance) {
        std::cout << "The system has converged! Let's run one more cycle to collect data" << std::endl;
        m_lastIteration = m_iter + 1;
    }
}

unsigned long System::adaptiveSteps() {
    unsigned long stepRatio = 1;
    if(m_iter == m_lastIteration+m_rangeOfDynamicSteps) {
        stepRatio = unsigned(pow(2,m_additionalStepsLastIteration));
        stepRatio = unsigned(pow(2,m_additionalStepsLastIteration));
    }
    else if(m_iter >= m_lastIteration) {
        stepRatio = unsigned(pow(2,m_additionalSteps));
        stepRatio = unsigned(pow(2,m_additionalSteps));
    }
    return stepRatio;
}

void System::updateAllArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix, const unsigned int changedCoord) {
    for(auto& i : m_waveFunctionElements) {
        i->setArrays();
        i->updateArrays(positions, radialVector, distanceMatrix, changedCoord);
    }
}

void System::resetAllArrays() {
    for(auto& i : m_waveFunctionElements) {
        i->resetArrays();
    }
}

void System::updateAllParameters(const Eigen::MatrixXd parameters) {
    for(unsigned long i=0; i<m_numberOfWaveFunctionElements; i++) {
        m_waveFunctionElements[i]->updateParameters(parameters, i);
    }
}

double System::evaluateWaveFunctionRatio() {
    double ratio = 1;
    for(auto& i : m_waveFunctionElements) {
        ratio *= i->evaluateRatio();
    }
    return ratio;
}

double System::getKineticEnergy() {
    double kineticEnergy = 0;
    for(auto& i : m_waveFunctionElements) {
        kineticEnergy += i->computeLaplacian();
    }
    for(unsigned int k = 0; k < m_numberOfFreeDimensions; k++) {
        double nablaLnPsi = 0;
        for(auto& i : m_waveFunctionElements) {
            nablaLnPsi += i->computeGradient(k);
        }
        kineticEnergy += nablaLnPsi * nablaLnPsi;
    }
    return - 0.5 * kineticEnergy;
}

Eigen::MatrixXd System::getAllInstantGradients() {
    Eigen::MatrixXd gradients = Eigen::MatrixXd::Zero(Eigen::Index(m_numberOfWaveFunctionElements), m_maxNumberOfParametersPerElement);
    for(unsigned long i = 0; i < m_numberOfWaveFunctionElements; i++) {
        gradients.row(Eigen::Index(i)) = m_waveFunctionElements[i]->computeParameterGradient();
    }
    return gradients;
}

void System::setGlobalArraysToCalculate() {
    // Check if the elements need distance matrix or radial distance vector
    for(auto& i : m_waveFunctionElements) {
        int need = i->getGlobalArrayNeed();
        if(need == 1) {
            m_calculateDistanceMatrix   = true;
        }
        if(need == 2) {
            m_calculateRadialVector     = true;
        }
        if(need == 3) {
            m_calculateDistanceMatrix   = true;
            m_calculateRadialVector     = true;
        }
    }
    // Check if the Hemiltonian needs distance matrix or radial distance vector
    unsigned int need = m_hamiltonian->getGlobalArrayNeed();
    if(need == 1) {
        m_calculateDistanceMatrix   = true;
    }
    if(need == 2) {
        m_calculateRadialVector     = true;
    }
    if(need == 3) {
        m_calculateDistanceMatrix   = true;
        m_calculateRadialVector     = true;
    }
}

void System::setNumberOfParticles(const unsigned int numberOfParticles) {
    assert(numberOfParticles > 0); // Check if the elements need distance matrix or radial distance vector
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(const unsigned short numberOfDimensions) {
    assert(numberOfDimensions > 0);
    assert(m_numberOfParticles > 0);
    m_numberOfDimensions = numberOfDimensions;
    setNumberOfFreeDimensions();
}

void System::setNumberOfFreeDimensions() {
    m_numberOfFreeDimensions = m_numberOfParticles * m_numberOfDimensions;
}

void System::setNumberOfHiddenNodes(const unsigned int numberOfHiddenNodes) {
    assert(numberOfHiddenNodes > 0);
    m_numberOfHiddenNodes = numberOfHiddenNodes;
}

void System::setNumberOfMetropolisSteps(const unsigned long steps) {
    // Calculate number of steps without equilibriation (power of 2)
    m_totalNumberOfStepsWOEqui         = steps;
    if(m_myRank == 0) {
        m_numberOfStepsWOEqui          = steps / m_numberOfProcesses + steps % m_numberOfProcesses;
    }
    else {
        m_numberOfStepsWOEqui          = steps / m_numberOfProcesses;
    }

    // Store the initial steps in case adaptive step is chosen
    m_initialNumberOfStepsWOEqui       = m_numberOfStepsWOEqui;
    m_initialTotalNumberOfStepsWOEqui  = m_totalNumberOfStepsWOEqui;

    // Calculate the number of equilibriation steps (needs to be unaffected by the number of processes)
    m_numberOfEquilibriationSteps      = unsigned(m_totalNumberOfStepsWOEqui * m_equilibrationFraction);
    m_totalNumberOfEquilibriationSteps = unsigned(m_totalNumberOfStepsWOEqui * m_equilibrationFraction * m_numberOfProcesses);

    // Calculate the number of steps included equilibriation
    m_totalNumberOfStepsWEqui          = m_totalNumberOfStepsWOEqui + m_totalNumberOfEquilibriationSteps;
    m_numberOfStepsWEqui               = m_numberOfStepsWOEqui + m_numberOfEquilibriationSteps;
}

void System::setNumberOfWaveFunctionElements(const unsigned long numberOfWaveFunctionElements) {
    m_numberOfWaveFunctionElements = numberOfWaveFunctionElements;
}

void System::setMaxNumberOfParametersPerElement() {
    unsigned int maxNumberOfWaveFunctionElements = 0;
    unsigned int counter = 0;
    for(auto& i : m_waveFunctionElements) {
        unsigned int numberOfParameters = i->getNumberOfParameters();
        if(numberOfParameters > maxNumberOfWaveFunctionElements) {
            maxNumberOfWaveFunctionElements = numberOfParameters;
        }
        counter += numberOfParameters;
    }
    m_maxNumberOfParametersPerElement = maxNumberOfWaveFunctionElements;
    m_totalNumberOfParameters = counter;
}

void System::setStepLength(const double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(const double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setFrequency(const double omega) {
    assert(omega > 0);
    m_omega = omega;
}

void System::setAtomicNumber(const unsigned int Z) {
    assert(Z > 0);
    m_Z = Z;
}

void System::setLearningRate(const double learningRate) {
    assert(learningRate > 0);
    m_eta = learningRate;
}

void System::setWidth(const double sigma) {
    assert(sigma > 0);
    m_sigma = sigma;
}

void System::setInteraction(const bool interaction) {
    m_interaction = interaction;
}

void System::setConvergenceTools(bool checkConvergence, unsigned int numberOfEnergies, double tolerance) {
    m_checkConvergence = checkConvergence;
    m_tolerance        = tolerance;
    m_numberOfEnergies = numberOfEnergies;
    m_energies         = Eigen::VectorXd::Zero(numberOfEnergies);
}

void System::setDynamicStepTools(bool applyAdaptiveSteps, unsigned int rangeOfDynamicSteps, unsigned int additionalSteps, unsigned int additionalStepsLastIteration) {
    m_applyAdaptiveSteps = applyAdaptiveSteps;
    m_rangeOfDynamicSteps = rangeOfDynamicSteps;
    m_additionalSteps = additionalSteps;
    m_additionalStepsLastIteration = additionalStepsLastIteration;
}

void System::setDensityTools(bool computeOneBodyDensity, bool computeTwoBodyDensity, int numberOfBins, double maxRadius) {
    m_computeOneBodyDensity = computeOneBodyDensity;
    m_computeTwoBodyDensity = computeTwoBodyDensity;
    m_numberOfBins          = numberOfBins;
    m_maxRadius             = maxRadius;
}

void System::setEnergyPrintingTools(bool printEnergyFile, bool printInstantEnergyFile) {
    m_printEnergyFile        = printEnergyFile;
    m_printInstantEnergyFile = printInstantEnergyFile;
}

void System::setParameterPrintingTools(bool printParametersToFile) {
    m_printParametersToFile        = printParametersToFile;
}

void System::setMPITools(int myRank, int numberOfProcesses) {
    m_myRank = myRank;
    m_numberOfProcesses = numberOfProcesses;
}

void System::setPath(const std::string path) {
    m_path = path;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setBasis(Basis* basis) {
    m_basis = basis;
}

void System::setWaveFunctionElements(std::vector<class WaveFunction *> waveFunctionElements) {
    m_waveFunctionElements = waveFunctionElements;
    setMaxNumberOfParametersPerElement();
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}

void System::setInitialWeights(InitialWeights* initialWeights) {
    m_initialWeights = initialWeights;
}

void System::setMetropolis(Metropolis* metropolis) {
    m_metropolis = metropolis;
}

void System::setOptimization(Optimization* optimization) {
    m_optimization = optimization;
}

void System::setRandomNumberGenerator(RandomNumberGenerator* randomNumberGenerator) {
    m_randomNumberGenerator = randomNumberGenerator;
}

void System::setGradients() {
    m_gradients = Eigen::MatrixXd::Zero(Eigen::Index(m_numberOfWaveFunctionElements), m_maxNumberOfParametersPerElement);
}

std::string System::getAllLabels() {
    std::string total_string = "";
    for(auto& i : m_waveFunctionElements) {
        total_string += "_";
        total_string += i->getLabel();
    }

    if((total_string == "_gaussian_padejastrow") || \
       (total_string == "_padejastrow_gaussian") || \
       (total_string == "_gaussian_padejastrow_slaterdeterminant") || \
       (total_string == "_gaussian_slaterdeterminant_padejastrow") || \
       (total_string == "_padejastrow_gaussian_slaterdeterminant") || \
       (total_string == "_padejastrow_slaterdeterminant_gaussian") || \
       (total_string == "_slaterdeterminant_gaussian_padejastrow") || \
       (total_string == "_slaterdeterminant_padejastrow_gaussian")) {
        total_string = "VMC";
    }
    else if((total_string == "_rbmgaussian_rbmjastrow") || \
            (total_string == "_rbmjastrow_rbmgaussian") || \
            (total_string == "_rbmgaussian_rbmjastrow_slaterdeterminant") || \
            (total_string == "_rbmgaussian_slaterdeterminant_rbmjastrow") || \
            (total_string == "_rbmjastrow_rbmgaussian_slaterdeterminant") || \
            (total_string == "_rbmjastrow_slaterdeterminant_rbmgaussian") || \
            (total_string == "_slaterdeterminant_rbmgaussian_rbmjastrow") || \
            (total_string == "_slaterdeterminant_rbmjastrow_rbmgaussian")) {
        total_string = "RBM";
    }
    else if((total_string == "_rbmgaussian_rbmjastrow_padejastrow") || \
            (total_string == "_rbmjastrow_rbmgaussian_padejastrow") || \
            (total_string == "_padejastrow_rbmjastrow_rbmgaussian") || \
            (total_string == "_padejastrow_rbmgaussian_rbmjastrow") || \
            (total_string == "_rbmjastrow_padejastrow_rbmgaussian") || \
            (total_string == "_rbmgaussian_padejastrow_rbmjastrow") || \
            (total_string == "_rbmgaussian_rbmjastrow_slaterdeterminant_padejastrow") || \
            (total_string == "_rbmgaussian_slaterdeterminant_rbmjastrow_padejastrow") || \
            (total_string == "_rbmjastrow_rbmgaussian_slaterdeterminant_padejastrow") || \
            (total_string == "_rbmjastrow_slaterdeterminant_rbmgaussian_padejastrow") || \
            (total_string == "_slaterdeterminant_rbmgaussian_rbmjastrow_padejastrow") || \
            (total_string == "_slaterdeterminant_rbmjastrow_rbmgaussian_padejastrow")) {
        total_string = "RBMPJ";
    }
    /*
    std::vector<std::string> elementLabels;
    elementLabels.push_back("_rbmgaussian");
    elementLabels.push_back("_rbmjastrow");
    elementLabels.push_back("_padejastrow");
    for(auto& i : elementLabels) {

    }
    */
    return total_string;
}
