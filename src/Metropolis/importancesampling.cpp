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
    m_quantumForceOld        = Eigen::VectorXd::Zero(m_numberOfFreeDimensions);
    for(int i=0; i<m_numberOfFreeDimensions; i++) {
        m_quantumForceOld(i) = QuantumForce(i);
    }
    m_quantumForceNew        = m_quantumForceOld;
}

double ImportanceSampling::QuantumForce(const int i) {
    double QF = 0;
    for(auto& j : m_waveFunctionVector) {
        QF += j->computeGradient(i);
    }
    return 2*QF;
}

double ImportanceSampling::GreenFuncSum() {
    double GreenSum  = 0;
    for(int i=0; i<m_numberOfParticles; i++) {
        double GreenFunc = 0;
        for(int j=0; j<m_numberOfDimensions; j++) {
            int l = m_numberOfDimensions*i+j;
            double QForceOld = m_quantumForceOld(l);
            double QForceNew = m_quantumForceNew(l);
            GreenFunc += 0.5 * (QForceOld + QForceNew) * (0.5 * m_diff*m_stepLength*(QForceOld - QForceNew) - m_positions(l)+m_positionsOld(l));
        }
        GreenSum += exp(GreenFunc);
    }
    return GreenSum;
}

bool ImportanceSampling::acceptMove() {
    int changedCoord  = m_system->getRandomNumberGenerator()->nextInt(m_numberOfFreeDimensions);

    m_quantumForceOld = m_quantumForceNew;
    m_positionsOld    = m_positions;

    m_positions(changedCoord) += m_diff * QuantumForce(changedCoord) * m_stepLength + m_system->getRandomNumberGenerator()->nextGaussian(0,1) * sqrt(m_stepLength);
    m_system->updateAllArrays(m_positions, changedCoord);
    m_quantumForceNew(changedCoord) = QuantumForce(changedCoord);

    double ratio = m_system->evaluateWaveFunctionRatio();
    double w = GreenFuncSum() * ratio;
    double r = m_system->getRandomNumberGenerator()->nextDouble();

    if(w < r) {
        m_system->resetAllArrays();
        m_positions       = m_positionsOld;
        m_quantumForceNew = m_quantumForceOld;
        return false;
    }
    return true;
}