#include "metropolis.h"
#include "../system.h"

Metropolis::Metropolis(System* system) {
    m_system = system;
}

double Metropolis::calculateDistanceMatrixElement(const unsigned int i, const unsigned int j) {
    double dist = 0;
    unsigned int parti   = m_numberOfDimensions*i;
    unsigned int partj   = m_numberOfDimensions*j;
    for(unsigned short d=0; d<m_numberOfDimensions; d++) {
        double diff = m_positions(parti+d)-m_positions(partj+d);
        dist += diff*diff;
    }
    return sqrt(dist);
}

void Metropolis::calculateDistanceMatrixCross(const unsigned int particle) {
    for(unsigned int i=0; i<m_numberOfParticles; i++) {
        m_distanceMatrix(particle, i) = calculateDistanceMatrixElement(particle, i);
        m_distanceMatrix(i, particle) = m_distanceMatrix(particle, i);
    }
}

double Metropolis::calculateRadialVectorElement(const unsigned int particle) {
    double sqrtElementWise = 0;
    unsigned int part = particle*m_numberOfDimensions;
    for(unsigned short d=0; d<m_numberOfDimensions; d++) {
        sqrtElementWise += m_positions(part + d) * m_positions(part + d);
    }
    return sqrt(sqrtElementWise);
}

Metropolis::~Metropolis() {};
