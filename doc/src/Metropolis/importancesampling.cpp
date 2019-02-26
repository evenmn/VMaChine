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
}

double ImportanceSampling::QuantumForce(const int i) {
    double QF = 0;
    for(auto& j : m_waveFunctionVector) {
        QF += j->computeFirstDerivative(i);
    }
    return 2*QF;
}

double ImportanceSampling::GreenFuncSum(const Eigen::VectorXd oldPositions, const Eigen::VectorXd QuantumForceOld) {
    double GreenSum  = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        double GreenFunc = 0;
        for(int j=0; j<m_numberOfDimensions; j++) {
            double QForceOld = QuantumForceOld(m_numberOfDimensions*i+j);
            double QForceNew = QuantumForce(m_numberOfDimensions*i+j);
            GreenFunc += 0.5 * (QForceOld + QForceNew) * (0.5 * m_diff*m_stepLength*(QForceOld - QForceNew) - m_positions(m_numberOfDimensions*i+j)+oldPositions(m_numberOfDimensions*i+j));
        }
        GreenSum += exp(GreenFunc);
    }
    return GreenSum;
}

bool ImportanceSampling::acceptMove() {
    double psiOld   = m_system->evaluateWaveFunctionSqrd();
    int pRand       = m_system->getRandomNumberGenerator()->nextInt(m_numberOfFreeDimensions);

    Eigen::VectorXd QuantumForceOld = Eigen::VectorXd::Zero(m_numberOfFreeDimensions);
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        QuantumForceOld(i) = QuantumForce(i);
    }

    Eigen::VectorXd oldPositions = m_positions;
    m_positions(pRand) += m_diff * QuantumForceOld(pRand) * m_stepLength + m_system->getRandomNumberGenerator()->nextGaussian(0,1) * sqrt(m_stepLength);

    m_system->updateAllArrays(m_positions, pRand);
    double psiNew = m_system->evaluateWaveFunctionSqrd();

    double w = GreenFuncSum(oldPositions, QuantumForceOld) * (psiNew/psiOld);
    double r = m_system->getRandomNumberGenerator()->nextDouble();
    if(w < r) {
        m_system->resetAllArrays();
        m_positions = oldPositions;
        return false;
    }
    return true;
}
