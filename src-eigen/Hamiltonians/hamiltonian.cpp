#include "hamiltonian.h"
#include "../system.h"

Hamiltonian::Hamiltonian(System *system)
{
    m_system = system;
}

Eigen::Map<Eigen::VectorXd> flatten(Eigen::MatrixXd A)
{
    return Eigen::Map<Eigen::VectorXd>(A.data(), A.size());
}

double Hamiltonian::getInteractionEnergy()
{
    m_distanceMatrix = m_system->getDistanceMatrix();
    double interactionEnergy = 0;
    if (m_interaction) {
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
            Eigen::MatrixXd upperTriangular = Eigen::MatrixXd(m_distanceMatrix.triangularView<Eigen::StrictlyUpper>());
            Eigen::MatrixXd distanceMatrixInverse = upperTriangular.cwiseInverse();
            Eigen::VectorXd elements = flatten(distanceMatrixInverse);
            elements = (elements.array().isFinite()).select(elements,0);
            interactionEnergy = elements.sum();
        }
    }
    return interactionEnergy;
}

Hamiltonian::~Hamiltonian() {}
