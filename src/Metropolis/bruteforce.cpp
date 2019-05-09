#include "bruteforce.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"
#include "../RNG/rng.h"
#include "InitialStates/initialstate.h"

BruteForce::BruteForce(System* system) :
        Metropolis(system) {
    m_numberOfDimensions     = m_system->getNumberOfDimensions();
    m_numberOfParticles      = m_system->getNumberOfParticles();
    m_numberOfFreeDimensions = m_system->getNumberOfFreeDimensions();
    m_stepLength             = m_system->getStepLength();
    m_positions              = m_system->getInitialState()->getParticles();
    m_radialVector           = m_system->getInitialState()->getRadialVector();
    m_distanceMatrix         = m_system->getInitialState()->getDistanceMatrix();
    m_calculateDistanceMatrix = m_system->getCalculateDistanceMatrix();
    m_calculateRadialVector  = m_system->getCalculateRadialVector();
}

bool BruteForce::acceptMove() {
    int pRand = m_system->getRandomNumberGenerator()->nextInt(m_numberOfFreeDimensions);

    m_positionsOld      = m_positions;
    m_radialVectorOld   = m_radialVector;
    m_distanceMatrixOld = m_distanceMatrix;

    m_positions(pRand) = m_positionsOld(pRand) + 10 * (m_system->getRandomNumberGenerator()->nextDouble() - 0.5) * m_stepLength;
    if(m_calculateDistanceMatrix) {
        Metropolis::calculateDistanceMatrixCross(int(pRand/m_numberOfDimensions));
    }
    if(m_calculateRadialVector) {
        Metropolis::calculateRadialVectorElement(int(pRand/m_numberOfDimensions));
    }
    m_system->updateAllArrays(m_positions, m_radialVector, m_distanceMatrix, pRand);

    double ratio = m_system->evaluateWaveFunctionRatio();
    double r = m_system->getRandomNumberGenerator()->nextDouble();
    if(ratio < r) {
        m_positions(pRand) = m_positionsOld(pRand);
        m_distanceMatrix   = m_distanceMatrixOld;
        m_radialVector     = m_radialVectorOld;
        m_system->resetAllArrays();
        return false;
    }
    return true;
}

double BruteForce::calculateDistanceMatrixElement(int i, int j) {
    double dist = 0;
    int parti   = m_numberOfDimensions*i;
    int partj   = m_numberOfDimensions*j;
    for(int d=0; d<m_numberOfDimensions; d++) {
        double diff = m_positions(parti+d)-m_positions(partj+d);
        dist += diff*diff;
    }
    return sqrt(dist);
}

void BruteForce::calculateDistanceMatrixCross(int particle) {
    for(int i=0; i<m_numberOfParticles; i++) {
        m_distanceMatrix(particle, i) = calculateDistanceMatrixElement(particle, i);
        m_distanceMatrix(i, particle) = m_distanceMatrix(particle, i);
    }
}

double BruteForce::calculateRadialVectorElement(int particle) {
    double sqrtElementWise = 0;
    int part = particle*m_numberOfDimensions;
    for(int d=0; d<m_numberOfDimensions; d++) {
        sqrtElementWise += m_positions(part + d) * m_positions(part + d);
    }
    return sqrt(sqrtElementWise);
}
