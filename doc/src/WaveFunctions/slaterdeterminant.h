#pragma once
#include "wavefunction.h"

class SlaterDeterminant : public WaveFunction {
public:
    SlaterDeterminant(class System* system, const int elementNumber);
    Eigen::MatrixXd list();
    double updateElement(const Eigen::VectorXd positions, double basis(double, int), const int i, const int j);
    Eigen::VectorXd updateRow(const Eigen::VectorXd positions, double basis(double, int), const int i);
    Eigen::MatrixXd updateMatrix(const Eigen::VectorXd positions, double basis(double, int));
    Eigen::VectorXd dA_row(const Eigen::VectorXd positions, const int k);
    Eigen::MatrixXd dA_matrix(const Eigen::VectorXd positions);
    Eigen::VectorXd d2A_row(const Eigen::VectorXd positions, const int k);
    Eigen::MatrixXd d2A_matrix(const Eigen::VectorXd positions);

    void updateArrays(const Eigen::VectorXd positions, const int pRand);
    void resetArrays();
    void initializeArrays(const Eigen::VectorXd positions);
    void updateParameters(const Eigen::MatrixXd parameters);
    double evaluate();
    double evaluateSqrd();
    double computeFirstDerivative(const Eigen::VectorXd positions, const int k);
    double computeSecondDerivative();
    Eigen::VectorXd computeFirstEnergyDerivative(const int k);
    Eigen::VectorXd computeSecondEnergyDerivative();
private:
    int     m_elementNumber     = 2;
    double  m_omega             = 0;
    double  m_alpha             = 0;
    int m_numberOfOrbitals      = 0;
    int m_numberOfParticlesHalf = 0;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_oldPositions;
    Eigen::MatrixXd m_D_up;
    Eigen::MatrixXd m_D_dn;
    Eigen::MatrixXd m_D_upOld;
    Eigen::MatrixXd m_D_dnOld;
    Eigen::MatrixXd m_dD_up;
    Eigen::MatrixXd m_dD_dn;
    Eigen::MatrixXd m_dD_upOld;
    Eigen::MatrixXd m_dD_dnOld;
    Eigen::MatrixXd m_d2D_up;
    Eigen::MatrixXd m_d2D_dn;
    Eigen::MatrixXd m_d2D_upOld;
    Eigen::MatrixXd m_d2D_dnOld;
    Eigen::MatrixXd m_D_up_inv;
    Eigen::MatrixXd m_D_dn_inv;
    Eigen::MatrixXd m_D_up_old_inv;
    Eigen::MatrixXd m_D_dn_old_inv;
};
