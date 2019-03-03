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

#include <iostream>
#include <fstream>
#include <cassert>
#include <ctime>
#include <string>

void System::runMetropolisSteps(const int numberOfIterations) {
    m_positions                 = m_initialState->getParticles();
    m_parameters                = m_initialWeights->getWeights();
    m_sampler                   = new Sampler(this);
    m_sampler->openOutputFiles("../data/");

    for (int iter = 0; iter < numberOfIterations; iter++) {
        clock_t start_time = clock();
        for (int i=0; i < m_totalNumberOfSteps; i++) {
            bool acceptedStep = m_metropolis->acceptMove();
            m_positions       = m_metropolis->updatePositions();
            if(i >= (m_totalNumberOfSteps - m_numberOfMetropolisSteps)) {
                m_sampler->sample(acceptedStep, i);
                if(iter == numberOfIterations-1) {
                    m_sampler->printImmediatelyToFile(m_positions);
                }
            }
        }
        clock_t end_time = clock();
        m_sampler->computeAverages();
        m_sampler->printOutputToTerminal(iter, numberOfIterations, double(end_time - start_time)/CLOCKS_PER_SEC);
        m_sampler->printOutputToFile();
        m_parameters -= m_optimization->updateParameters();
        updateAllParameters(m_parameters);
    }
    m_sampler->closeOutputFiles();
}

void System::setNumberOfParticles(const int numberOfParticles) {
    assert(numberOfParticles > 0);
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(const int numberOfDimensions) {
    assert(numberOfDimensions > 0);
    m_numberOfDimensions = numberOfDimensions;
}

void System::setNumberOfFreeDimensions() {
    m_numberOfFreeDimensions = m_numberOfParticles * m_numberOfDimensions;
}

void System::setNumberOfHiddenNodes(const int numberOfHiddenNodes) {
    assert(numberOfHiddenNodes > 0);
    m_numberOfHiddenNodes = numberOfHiddenNodes;
}

void System::setNumberOfMetropolisSteps(const int steps) {
    m_numberOfMetropolisSteps = steps;
}

void System::setNumberOfWaveFunctionElements(const int numberOfWaveFunctionElements) {
    m_numberOfWaveFunctionElements = numberOfWaveFunctionElements;
}

void System::setMaxNumberOfParametersPerElement(const int maxNumberOfParametersPerElement) {
    m_maxNumberOfParametersPerElement = maxNumberOfParametersPerElement;
}

void System::setStepLength(const double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setTotalNumberOfSteps() {
    m_totalNumberOfSteps = int(m_numberOfMetropolisSteps*(1 + m_equilibrationFraction));
}

void System::setEquilibrationFraction(const double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setInteraction(const bool interaction) {
    m_interaction = interaction;
}

void System::setFrequency(const double omega) {
    assert(omega > 0);
    m_omega = omega;
}

void System::setLearningRate(const double eta) {
    assert(eta > 0);
    m_eta = eta;
}

void System::setWidth(const double sigma) {
    assert(sigma > 0);
    m_sigma = sigma;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setBasis(Basis* basis) {
    m_basis = basis;
}

void System::setWaveFunction(std::vector<class WaveFunction *> waveFunctionVector) {
    m_waveFunctionVector = waveFunctionVector;
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

void System::updateAllArrays(const Eigen::VectorXd positions, const int pRand) {
    for(auto& i : m_waveFunctionVector) {
        i->updateArrays(positions, pRand);
    }
}

void System::resetAllArrays() {
    for(auto& i : m_waveFunctionVector) {
        i->resetArrays();
    }
}

void System::updateAllParameters(const Eigen::MatrixXd parameters) {
    for(int i=0; i<m_numberOfWaveFunctionElements; i++) {
        m_waveFunctionVector[unsigned(i)]->updateParameters(parameters, i);
    }
}

double System::evaluateWaveFunction() {
    double WF = 1;
    for(auto& i : m_waveFunctionVector) {
        WF *= i->evaluate();
    }
    return WF;
}

double System::evaluateWaveFunctionSqrd() {
    double WF = 1;
    for(auto& i : m_waveFunctionVector) {
        WF *= i->evaluateSqrd();
    }
    return WF;
}

double System::getKineticEnergy() {
    double KineticEnergy = 0;
    for(auto& i : m_waveFunctionVector) {
        KineticEnergy += i->computeSecondDerivative();
    }
    for(int k = 0; k < m_numberOfFreeDimensions; k++) {
        double NablaLnPsi = 0;
        for(auto& i : m_waveFunctionVector) {
            NablaLnPsi += i->computeFirstDerivative(k);
        }
        KineticEnergy += NablaLnPsi * NablaLnPsi;
    }

    return -0.5 * KineticEnergy;
}
