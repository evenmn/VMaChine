#include "metropolis.h"
#include "../system.h"

Metropolis::Metropolis(System* system) {
    m_system = system;
}

double Metropolis::calculateDistanceMatrixElement(const int i, const int j) {
    double dist = 0;
    int parti   = m_numberOfDimensions*i;
    int partj   = m_numberOfDimensions*j;
    for(int d=0; d<m_numberOfDimensions; d++) {
        double diff = m_positions(parti+d)-m_positions(partj+d);
        dist += diff*diff;
    }
    return sqrt(dist);
}

void Metropolis::calculateDistanceMatrixCross(const int particle) {
    for(int i=0; i<m_numberOfParticles; i++) {
        m_distanceMatrix(particle, i) = calculateDistanceMatrixElement(particle, i);
        m_distanceMatrix(i, particle) = m_distanceMatrix(particle, i);
    }
}

double Metropolis::calculateRadialVectorElement(const int particle) {
    double sqrtElementWise = 0;
    int part = particle*m_numberOfDimensions;
    for(int d=0; d<m_numberOfDimensions; d++) {
        sqrtElementWise += m_positions(part + d) * m_positions(part + d);
    }
    return sqrt(sqrtElementWise);
}

Metropolis::~Metropolis() {}
