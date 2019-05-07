#include "importancesampling.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"
#include "../RNG/mersennetwister.h"
#include "InitialStates/initialstate.h"

using std::cout;
using std::endl;

ImportanceSampling::ImportanceSampling(System* system) :
        Metropolis(system) {
    m_numberOfParticles      = m_system->getNumberOfParticles();
    m_numberOfDimensions     = m_system->getNumberOfDimensions();
    m_numberOfFreeDimensions = m_system->getNumberOfFreeDimensions();
    m_stepLength             = m_system->getStepLength();
    m_waveFunctionVector     = m_system->getWaveFunctionElements();
    m_positions              = m_system->getInitialState()->getParticles();
    m_radialVector           = m_system->getInitialState()->getRadialVector();
    m_distanceMatrix         = m_system->getInitialState()->getDistanceMatrix();
    m_quantumForceOld        = Eigen::VectorXd::Zero(m_numberOfFreeDimensions);
    for(unsigned int i=0; i<m_numberOfFreeDimensions; i++) {
        m_quantumForceOld(i) = QuantumForce(i);
    }
    m_quantumForceNew        = m_quantumForceOld;
    m_calculateDistanceMatrix = m_system->getCalculateDistanceMatrix();
    m_calculateRadialVector  = m_system->getCalculateRadialVector();
}

double ImportanceSampling::QuantumForce(const unsigned int i) {
    double QF = 0;
    for(auto& j : m_waveFunctionVector) {
        QF += j->computeGradient(i);
    }
    return 2*QF;
}

double ImportanceSampling::GreenFuncSum() {
    double GreenSum  = 0;
    for(unsigned int i=0; i<m_numberOfParticles; i++) {
        double GreenFunc = 0;
        for(unsigned short d=0; d<m_numberOfDimensions; d++) {
            unsigned int l = m_numberOfDimensions*i+d;
            double QForceOld = m_quantumForceOld(l);
            double QForceNew = m_quantumForceNew(l);
            GreenFunc += 0.5 * (QForceOld + QForceNew) * (0.5 * m_diff*m_stepLength*(QForceOld - QForceNew) - m_positions(l)+m_positionsOld(l));
        }
        GreenSum += exp(GreenFunc);
    }
    return GreenSum;
}

bool ImportanceSampling::acceptMove() {
    unsigned int changedCoord  = m_system->getRandomNumberGenerator()->nextInt(m_numberOfFreeDimensions);

    m_quantumForceOld = m_quantumForceNew;
    m_positionsOld    = m_positions;
    m_radialVectorOld = m_radialVector;
    m_distanceMatrixOld = m_distanceMatrix;

    m_positions(changedCoord) += m_diff * QuantumForce(changedCoord) * m_stepLength + m_system->getRandomNumberGenerator()->nextGaussian(0,1) * sqrt(m_stepLength);
    if(m_calculateDistanceMatrix) {
        Metropolis::calculateDistanceMatrixCross((unsigned int)(changedCoord/m_numberOfDimensions));
    }
    if(m_calculateRadialVector) {
        Metropolis::calculateRadialVectorElement((unsigned int)(changedCoord/m_numberOfDimensions));
    }
    m_system->updateAllArrays(m_positions, m_radialVector, m_distanceMatrix, changedCoord);
    m_quantumForceNew(changedCoord) = QuantumForce(changedCoord);

    double ratio = m_system->evaluateWaveFunctionRatio();
    double w = GreenFuncSum() * ratio;
    double r = m_system->getRandomNumberGenerator()->nextDouble();

    if(w < r) {
        m_system->resetAllArrays();
        m_positions       = m_positionsOld;
        m_quantumForceNew = m_quantumForceOld;
        m_distanceMatrix  = m_distanceMatrixOld;
        m_radialVector    = m_radialVectorOld;
        return false;
    }
    return true;
}

double ImportanceSampling::calculateDistanceMatrixElement(const unsigned int i, const unsigned int j) {
    double dist = 0;
    unsigned int parti   = m_numberOfDimensions*i;
    unsigned int partj   = m_numberOfDimensions*j;
    for(unsigned short d=0; d<m_numberOfDimensions; d++) {
        double diff = m_positions(parti+d)-m_positions(partj+d);
        dist += diff*diff;
    }
    return sqrt(dist);
}

void ImportanceSampling::calculateDistanceMatrixCross(const unsigned int particle) {
    for(unsigned int i=0; i<m_numberOfParticles; i++) {
        m_distanceMatrix(particle, i) = calculateDistanceMatrixElement(particle, i);
        m_distanceMatrix(i, particle) = m_distanceMatrix(particle, i);
    }
}

double ImportanceSampling::calculateRadialVectorElement(const unsigned int particle) {
    double sqrtElementWise = 0;
    unsigned int part = particle*m_numberOfDimensions;
    for(unsigned short d=0; d<m_numberOfDimensions; d++) {
        sqrtElementWise += m_positions(part + d) * m_positions(part + d);
    }
    return sqrt(sqrtElementWise);
}
