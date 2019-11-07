#include "importancesampling.h"
#include "../InitialStates/initialstate.h"
#include "../RNG/mersennetwister.h"
#include "../WaveFunctions/wavefunction.h"
#include "../system.h"
#include <cassert>
#include <iostream>

ImportanceSampling::ImportanceSampling(System *system)
    : Metropolis(system)
{}

void ImportanceSampling::initialize() {
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_degreesOfFreedom = m_system->getNumberOfFreeDimensions();
    m_stepLength = m_system->getStepLength();
    m_positions = m_system->getInitialState()->getParticles();
    m_radialVector = m_system->getInitialState()->getRadialVector();
    m_distanceMatrix = m_system->getInitialState()->getDistanceMatrix();
    m_system->setGlobalArraysToCalculate();
    m_calculateDistanceMatrix = m_system->getCalculateDistanceMatrix();
    m_calculateRadialVector = m_system->getCalculateRadialVector();
    m_RNG = m_system->getRandomNumberGenerator();
    m_waveFunctionVector = m_system->getWaveFunctionElements();
    m_dtD = m_stepLength * m_diff;
    m_sqrtStep = sqrt(m_stepLength);
    initializeQuantumForce();
}

bool ImportanceSampling::acceptMove()
{
    arma::uword i = m_RNG->nextInt(m_degreesOfFreedom);

    m_quantumForceOld = m_quantumForceNew;
    m_positionsOld = m_positions;
    m_radialVectorOld = m_radialVector;
    m_distanceMatrixOld = m_distanceMatrix;

    m_quantumForceNew(i) = QuantumForce(i);
    m_dx = m_dtD * m_quantumForceNew(i) + m_RNG->nextGaussian(0, 1) * m_sqrtStep;
    m_positions(i) += m_dx;
    if (m_calculateDistanceMatrix) {
        Metropolis::calculateDistanceMatrixCross(arma::uword(i / m_numberOfDimensions));
    }
    if (m_calculateRadialVector) {
        Metropolis::calculateRadialVectorElement(arma::uword(i / m_numberOfDimensions));
    }

    m_system->updateAllArrays(m_positions, m_radialVector, m_distanceMatrix, i);

    double p = m_system->evaluateProbabilityRatio();
    double w = GreenRatio(i) * p;
    if (w < m_RNG->nextDouble()) {
        m_positions = m_positionsOld;
        m_quantumForceNew = m_quantumForceOld;
        m_distanceMatrix = m_distanceMatrixOld;
        m_radialVector = m_radialVectorOld;
        m_system->resetAllArrays();
        return false;
    }
    return true;
}

void ImportanceSampling::initializeQuantumForce()
{
    m_quantumForceNew.zeros(m_degreesOfFreedom);
    for (arma::uword i = 0; i < m_degreesOfFreedom; i++) {
        m_quantumForceNew(i) = QuantumForce(i);
    }
}

double ImportanceSampling::QuantumForce(const arma::uword i)
{
    double QF = 0;
    for (auto &j : m_waveFunctionVector) {
        QF += j->computeGradient(i);
    }
    return 2 * QF;
}

double ImportanceSampling::GreenRatio(const arma::uword i)
{
    double dQF = m_quantumForceOld(i) - m_quantumForceNew(i);
    return exp(0.5 * dQF * m_dx) + 1;
}
