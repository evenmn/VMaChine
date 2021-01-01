#include "coulomb.h"

Coulomb::Coulomb(System *system)
    : Interaction(system)
{}

Eigen::Map<Eigen::VectorXd> flat(Eigen::MatrixXd A)
{
    return Eigen::Map<Eigen::VectorXd>(A.data(), A.size());
}

void Coulomb::initialize()
{
    m_numberOfParticles = m_system->getNumberOfParticles();
}

double Coulomb::getInteractionEnergy()
{
    m_distanceMatrix = m_system->getDistanceMatrix();
    double interactionEnergy = 0;
    if (m_screening) {
        for (int i = 0; i < m_numberOfParticles; i++) {
            for (int j = i + 1; j < m_numberOfParticles; j++) {
                if (m_distanceMatrix(i, j) > m_dsl) {
                    interactionEnergy += exp(-m_distanceMatrix(i, j) / m_screeningStrength)
                                         / m_distanceMatrix(i, j);
                } else {
                    interactionEnergy += 1 / m_distanceMatrix(i, j);
                }
            }
        }
    }
    else {
        for (int i = 0; i < m_numberOfParticles; i++) {
            for (int j = i + 1; j < m_numberOfParticles; j++) {
                interactionEnergy += 1 / m_distanceMatrix(i, j);
            }
        }
        /*
        Eigen::MatrixXd upperTriangular = Eigen::MatrixXd(m_distanceMatrix.triangularView<Eigen::StrictlyUpper>());
        Eigen::MatrixXd distanceMatrixInverse = upperTriangular.cwiseInverse();
        Eigen::VectorXd elements = flat(distanceMatrixInverse);
        elements = (elements.array().isFinite()).select(elements,0);
        interactionEnergy = elements.sum();
        */
    }
    return interactionEnergy;
}
