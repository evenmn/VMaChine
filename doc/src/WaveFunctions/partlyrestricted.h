#pragma once
#include "wavefunction.h"

class PartlyRestricted : public WaveFunction {
public:
    PartlyRestricted(class System* system, const int elementNumber);
    void updateArrays(const Eigen::VectorXd positions, const int pRand);
    void resetArrays();
    void initializeArrays(const Eigen::VectorXd positions);
    void updateParameters(Eigen::MatrixXd parameters);
    double evaluate();
    double evaluateSqrd();
    double computeFirstDerivative(const Eigen::VectorXd positions, const int k);
    double computeSecondDerivative();
    Eigen::VectorXd computeFirstEnergyDerivative(const int k);
    Eigen::VectorXd computeSecondEnergyDerivative();
private:
    int     m_elementNumber = 0;
    double  m_sigmaSqrd2    = 1;
    double  m_oldXCx;
    double  m_xCx;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_oldPositions;
    Eigen::MatrixXd m_c;
};
