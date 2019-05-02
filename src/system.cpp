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

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <ctime>
#include <string>

void System::runIterations(const int numberOfIterations) {
    m_positions                 = m_initialState->getParticles();
    m_distanceMatrix            = m_initialState->getDistanceMatrix();
    m_radialVector              = m_initialState->getRadialVector();
    m_parameters                = m_initialWeights->getWeights();
    m_sampler                   = new Sampler(this);
    m_sampler->openOutputFiles();
    m_lastIteration = numberOfIterations - m_rangeOfDynamicSteps - 1;

    for(m_iter = 0; m_iter < numberOfIterations; m_iter++) {
        if(m_applyDynamicSteps) {
            m_numberOfStepsWOEqui      = m_initialNumberOfStepsWOEqui * dynamicSteps();
            m_numberOfStepsWEqui       = m_numberOfStepsWOEqui + m_numberOfEquilibriationSteps;
            m_totalNumberOfStepsWOEqui = m_initialTotalNumberOfStepsWOEqui * dynamicSteps();
            m_totalNumberOfStepsWEqui  = m_totalNumberOfStepsWOEqui + m_numberOfEquilibriationSteps;
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
        }
        m_sampler->printEnergyToFile();
        if(m_iter == m_lastIteration + m_rangeOfDynamicSteps) {
            m_sampler->printOneBodyDensityToFile();
        }
        printToTerminal(numberOfIterations);

        if(m_checkConvergence && m_myRank == 0) {
            checkingConvergence();
        }

        for(int i=0; i<m_numberOfWaveFunctionElements; i++) {
            for(int j=0; j<m_maxNumberOfParametersPerElement; j++) {
                MPI_Bcast(&m_parameters(i,j), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        updateAllParameters(m_parameters);
    }
}

void System::runMetropolisCycles() {
    for(int i=0; i < m_numberOfStepsWEqui; i++) {
        bool acceptedStep = m_metropolis->acceptMove();
        m_positions       = m_metropolis->updatePositions();
        m_distanceMatrix  = m_metropolis->updateDistanceMatrix();
        m_radialVector    = m_metropolis->updateRadialVector();
        if(i >= m_numberOfEquilibriationSteps) {
            m_sampler->sample(acceptedStep, i);
            if(m_iter == m_lastIteration + m_rangeOfDynamicSteps) {
                m_sampler->printInstantValuesToFile();
                m_sampler->calculateOneBodyDensities(m_positions);
            }
        }
    }
}

void System::printToTerminal(int numberOfIterations) {
    if(m_iter == m_lastIteration + m_rangeOfDynamicSteps) {
        m_sampler->closeOutputFiles();
        if(m_myRank == 0) {
            m_sampler->printFinalOutputToTerminal();
            exit(0);
        }
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

int System::dynamicSteps() {
    int stepRatio = 1;
    if(m_iter == m_lastIteration+m_rangeOfDynamicSteps) {
        stepRatio = int(pow(2,m_additionalStepsLastIteration));
        stepRatio = int(pow(2,m_additionalStepsLastIteration));
    }
    else if(m_iter >= m_lastIteration) {
        stepRatio = int(pow(2,m_additionalSteps));
        stepRatio = int(pow(2,m_additionalSteps));
    }
    return stepRatio;
}

void System::updateAllArrays(const Eigen::VectorXd positions, const Eigen::VectorXd radialVector, const Eigen::MatrixXd distanceMatrix, const int changedCoord) {
    for(auto& i : m_waveFunctionElements) {
        i->updateArrays(positions, radialVector, distanceMatrix, changedCoord);
    }
}

void System::resetAllArrays() {
    for(auto& i : m_waveFunctionElements) {
        i->resetArrays();
    }
}

void System::updateAllParameters(const Eigen::MatrixXd parameters) {
    for(int i=0; i<m_numberOfWaveFunctionElements; i++) {
        m_waveFunctionElements[unsigned(i)]->updateParameters(parameters, i);
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
    for(int k = 0; k < m_numberOfFreeDimensions; k++) {
        double nablaLnPsi = 0;
        for(auto& i : m_waveFunctionElements) {
            nablaLnPsi += i->computeGradient(k);
        }
        kineticEnergy += nablaLnPsi * nablaLnPsi;
    }
    return - 0.5 * kineticEnergy;
}

Eigen::MatrixXd System::getAllInstantGradients() {
    Eigen::MatrixXd gradients = Eigen::MatrixXd::Zero(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
    for(int i = 0; i < m_numberOfWaveFunctionElements; i++) {
        gradients.row(i) = m_waveFunctionElements[unsigned(i)]->computeParameterGradient();
    }
    return gradients;
}

void System::setGlobalArraysToCalculate() {
    // Check if the elements need distance matrix or radial distance vector
    for(auto& i : m_waveFunctionElements) {
        int need = i->getGlobalArrayNeed();
        if(need == 1) {
            m_calculateDistanceMatrix = true;
        }
        if(need == 2) {
            m_calculateRadialVector = true;
        }
        if(need == 3) {
            m_calculateDistanceMatrix = true;
            m_calculateRadialVector = true;
        }
    }
    // Check if the Hemiltonian needs distance matrix or radial distance vector
    int need = m_hamiltonian->getGlobalArrayNeed();
    if(need == 1) {
        m_calculateDistanceMatrix = true;
    }
    if(need == 2) {
        m_calculateRadialVector = true;
    }
    if(need == 3) {
        m_calculateDistanceMatrix = true;
        m_calculateRadialVector = true;
    }
}

void System::setNumberOfParticles(const int numberOfParticles) {
    assert(numberOfParticles > 0);// Check if the elements need distance matrix or radial distance vector
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(const int numberOfDimensions) {
    assert(numberOfDimensions > 0);
    assert(m_numberOfParticles > 0);
    m_numberOfDimensions = numberOfDimensions;
    setNumberOfFreeDimensions();
}

void System::setNumberOfFreeDimensions() {
    m_numberOfFreeDimensions = m_numberOfParticles * m_numberOfDimensions;
}

void System::setNumberOfHiddenNodes(const int numberOfHiddenNodes) {
    assert(numberOfHiddenNodes > 0);
    m_numberOfHiddenNodes = numberOfHiddenNodes;
}

void System::setNumberOfMetropolisSteps(const int steps) {
    m_totalNumberOfStepsWOEqui      = steps;
    m_totalNumberOfStepsWEqui       = int(steps * (1 + m_equilibrationFraction));
    if(m_myRank == 0) {
        m_numberOfStepsWOEqui       = steps / m_numberOfProcesses + steps % m_numberOfProcesses;
    }
    else {
        m_numberOfStepsWOEqui       = steps / m_numberOfProcesses;
    }
    m_numberOfStepsWEqui            = int(m_numberOfStepsWOEqui * (1 + m_equilibrationFraction));
    m_numberOfEquilibriationSteps   = m_numberOfStepsWEqui - m_numberOfStepsWOEqui;
    m_initialNumberOfStepsWOEqui    = m_numberOfStepsWOEqui;
    m_initialTotalNumberOfStepsWOEqui = m_totalNumberOfStepsWOEqui;
}

void System::setNumberOfWaveFunctionElements(const int numberOfWaveFunctionElements) {
    m_numberOfWaveFunctionElements = numberOfWaveFunctionElements;
}

void System::setMaxNumberOfParametersPerElement() {
    int maxNumberOfWaveFunctionElements = 0;
    int counter = 0;
    for(auto& i : m_waveFunctionElements) {
        int numberOfParameters = i->getNumberOfParameters();
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

void System::setAtomicNumber(const int Z) {
    assert(Z > 0);
    m_Z = Z;
}

void System::setLearningRate(const double eta) {
    assert(eta > 0);
    m_eta = eta;
}

void System::setWidth(const double sigma) {
    assert(sigma > 0);
    m_sigma = sigma;
}

void System::setInteraction(const bool interaction) {
    m_interaction = interaction;
}

void System::setConvergenceTools(bool checkConvergence, int numberOfEnergies, double tolerance) {
    m_checkConvergence = checkConvergence;
    m_tolerance        = tolerance;
    m_numberOfEnergies = numberOfEnergies;
    m_energies         = Eigen::VectorXd::Zero(numberOfEnergies);
}

void System::setDynamicStepTools(bool applyDynamicSteps, int rangeOfDynamicSteps, int additionalSteps, int additionalStepsLastIteration) {
    m_applyDynamicSteps = applyDynamicSteps;
    m_rangeOfDynamicSteps = rangeOfDynamicSteps;
    m_additionalSteps = additionalSteps;
    m_additionalStepsLastIteration = additionalStepsLastIteration;
}

void System::setDensityTools(bool computeDensity, int numberOfBins, double maxRadius) {
    m_computeDensity    = computeDensity;
    m_numberOfBins      = numberOfBins;
    m_maxRadius         = maxRadius;
}

void System::setEnergyPrintingTools(bool printEnergyFile, bool printInstantEnergyFile) {
    m_printEnergyFile        = printEnergyFile;
    m_printInstantEnergyFile = printInstantEnergyFile;
}

void System::setMPITools(int myRank, int numberOfProcesses) {
    m_myRank = myRank;
    m_numberOfProcesses = numberOfProcesses;
}

void System::setPath(std::string path) {
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

void System::setRandomNumberGenerator(RandomNumberGenerator* randomnumbergenerator) {
    m_randomnumbergenerator = randomnumbergenerator;
}

void System::setGradients() {
    m_gradients = Eigen::MatrixXd::Zero(m_numberOfWaveFunctionElements, m_maxNumberOfParametersPerElement);
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
