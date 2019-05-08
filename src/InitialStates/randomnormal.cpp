#include "randomnormal.h"
#include <iostream>
#include <cassert>
#include "RNG/rng.h"
#include "../system.h"
#include "WaveFunctions/wavefunction.h"

RandomNormal::RandomNormal(System*    system)  :
        InitialState(system) {
    m_numberOfParticles      = m_system->getNumberOfParticles();
    m_numberOfDimensions     = m_system->getNumberOfDimensions();
    m_numberOfFreeDimensions = m_system->getNumberOfFreeDimensions();
    setupInitialState();
}

double RandomNormal::calculateDistanceMatrixElement(const unsigned int i, const unsigned int j) {
    double dist = 0;
    unsigned int parti   = m_numberOfDimensions*i;
    unsigned int partj   = m_numberOfDimensions*j;
    for(unsigned short d=0; d<m_numberOfDimensions; d++) {
        double diff = m_positions(parti+d)-m_positions(partj+d);
        dist += diff*diff;
    }
    return sqrt(dist);
}

void RandomNormal::calculateDistanceMatrix() {
    m_distanceMatrix = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);
    for(unsigned int i=0; i<m_numberOfParticles; i++) {
        for(unsigned int j=i+1; j<m_numberOfParticles; j++) {
            m_distanceMatrix(i,j) = calculateDistanceMatrixElement(i,j);
            m_distanceMatrix(j,i) = m_distanceMatrix(i,j);
        }
    }
}

double RandomNormal::calculateRadialVectorElement(const unsigned int particle) {
    double sqrtElementWise = 0;
    unsigned int part = particle*m_numberOfDimensions;
    for(unsigned short d=0; d<m_numberOfDimensions; d++) {
        sqrtElementWise += m_positions(part + d) * m_positions(part + d);
    }
    return sqrt(sqrtElementWise);
}

void RandomNormal::calculateRadialVector() {
    m_radialVector = Eigen::VectorXd::Zero(m_numberOfParticles);
    for(unsigned int i=0; i<m_numberOfParticles; i++) {
        m_radialVector(i) = calculateRadialVectorElement(i);
    }
}

void RandomNormal::setupInitialState() {
    m_positions = Eigen::VectorXd::Zero(m_numberOfFreeDimensions);
    for (unsigned int i=0; i < m_numberOfFreeDimensions; i++) {
        m_positions(i) = m_system->getRandomNumberGenerator()->nextGaussian(0,1);
    }
    calculateDistanceMatrix();
    calculateRadialVector();
    for(auto& i : m_system->getWaveFunctionElements()) {
        i->initializeArrays(m_positions, m_radialVector, m_distanceMatrix);
        i->setArrays();
    }
}
