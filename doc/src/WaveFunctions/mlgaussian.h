#pragma once
#include "wavefunction.h"

class MLGaussian : public WaveFunction {
public:
    MLGaussian(class System* system, int elementNumber);
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
    int     m_elementNumber = 0;
    double  m_omega         = 1;
    double  m_sigmaSqrd     = 1;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_oldPositions;
    Eigen::VectorXd m_Xa;
    Eigen::VectorXd m_oldXa;
    Eigen::VectorXd m_a;
};
