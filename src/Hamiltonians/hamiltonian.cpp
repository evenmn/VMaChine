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

Hamiltonian::~Hamiltonian() {}
