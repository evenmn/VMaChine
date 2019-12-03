#include "initialstate.h"

InitialState::InitialState(System *system)
{
    m_system = system;
}

double InitialState::calculateDistanceMatrixElement(int i, int j)
{
    double dist = 0;
    int parti = m_numberOfDimensions * i;
    int partj = m_numberOfDimensions * j;
    for (int d = 0; d < m_numberOfDimensions; d++) {
        double diff = m_positions(parti + d) - m_positions(partj + d);
        dist += diff * diff;
    }
    return sqrt(dist);
}

void InitialState::calculateDistanceMatrix()
{
    m_distanceMatrix = Eigen::MatrixXd::Zero(m_numberOfParticles, m_numberOfParticles);
    for (int i = 0; i < m_numberOfParticles; i++) {
        for (int j = i + 1; j < m_numberOfParticles; j++) {
            m_distanceMatrix(i, j) = calculateDistanceMatrixElement(i, j);
            //m_distanceMatrix(j, i) = m_distanceMatrix(i, j);
        }
    }
}

double InitialState::calculateRadialVectorElement(int particle)
{
    double sqrtElementWise = 0;
    int part = particle * m_numberOfDimensions;
    for (int d = 0; d < m_numberOfDimensions; d++) {
        sqrtElementWise += m_positions(part + d) * m_positions(part + d);
    }
    return sqrt(sqrtElementWise);
}

void InitialState::calculateRadialVector()
{
    m_radialVector = Eigen::VectorXd::Zero(m_numberOfParticles);
    for (int i = 0; i < m_numberOfParticles; i++) {
        m_radialVector(i) = calculateRadialVectorElement(i);
    }
}

InitialState::~InitialState() {}