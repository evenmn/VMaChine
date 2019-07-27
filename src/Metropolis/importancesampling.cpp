#include "importancesampling.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../WaveFunctions/wavefunction.h"
#include "../RNG/mersennetwister.h"
#include "../InitialStates/initialstate.h"

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
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        m_quantumForceOld(i) = QuantumForce(i);
    }
    m_quantumForceNew        = m_quantumForceOld;

    system->setGlobalArraysToCalculate();
    m_calculateDistanceMatrix = m_system->getCalculateDistanceMatrix();
    m_calculateRadialVector  = m_system->getCalculateRadialVector();
    m_dtD = m_stepLength * m_diff;
    m_sqrtStep = sqrt(m_stepLength);
}

double ImportanceSampling::QuantumForce(const int i) {
    double QF = 0;
    for(auto& j : m_waveFunctionVector) {
        QF += j->computeGradient(i);
    }
    return 2*QF;
}

double ImportanceSampling::GreenFuncSum() {
    //double GreenFunc = (m_quantumForceOld(m_changedCoord) - m_quantumForceNew(m_changedCoord)) * (m_positions(m_changedCoord)-m_positionsOld(m_changedCoord));
    //return exp(0.5 * GreenFunc);

    double GreenSum = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        double GreenFunc = 0;
        for(int j=0; j<m_numberOfDimensions; j++) {
            int l = m_numberOfDimensions*i+j;
            double QForceOld = m_quantumForceOld(l);
            double QForceNew = m_quantumForceNew(l);
            //GreenFunc += (QForceOld + QForceNew) * (0.5 * m_dtD*(QForceOld - QForceNew) - m_positions(l)+m_positionsOld(l));
            GreenFunc += (m_quantumForceOld(l) - m_quantumForceNew(l)) * (m_positions(l)-m_positionsOld(l));
        }
        GreenSum += exp(0.5 * GreenFunc);
    }
    return GreenSum;
}

bool ImportanceSampling::acceptMove() {
    m_changedCoord  = m_system->getRandomNumberGenerator()->nextInt(m_numberOfFreeDimensions);

    m_quantumForceOld = m_quantumForceNew;
    m_positionsOld    = m_positions;
    m_radialVectorOld = m_radialVector;
    m_distanceMatrixOld = m_distanceMatrix;

    m_positions(m_changedCoord) += m_dtD * QuantumForce(m_changedCoord) + m_system->getRandomNumberGenerator()->nextGaussian(0,1) * m_sqrtStep;
    if(m_calculateDistanceMatrix) {
        Metropolis::calculateDistanceMatrixCross(int(m_changedCoord/m_numberOfDimensions));
    }
    if(m_calculateRadialVector) {
        Metropolis::calculateRadialVectorElement(int(m_changedCoord/m_numberOfDimensions));
    }

    m_system->updateAllArrays(m_positions, m_radialVector, m_distanceMatrix, m_changedCoord);
    m_quantumForceNew(m_changedCoord) = QuantumForce(m_changedCoord);

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
