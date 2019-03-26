#pragma once
#include "wavefunction.h"

class Gaussian : public WaveFunction {
public:
    Gaussian(class System* system);
    void updateArrays(const Eigen::VectorXd positions, const int pRand);
    void resetArrays();
    void initializeArrays(const Eigen::VectorXd positions);
    void updateParameters(const Eigen::MatrixXd parameters, const int elementNumber);
    double evaluateRatio();
    double computeFirstDerivative(const int k);
    double computeSecondDerivative();
    Eigen::VectorXd computeFirstEnergyDerivative(const int k);
    Eigen::VectorXd computeSecondEnergyDerivative();
private:
    int     m_elementNumber = 0;
    double  m_omega         = 0;
    double  m_alpha         = 0;
    double  m_probabilityRatio         = 0;
    double  m_probabilityRatioOld      = 0;
    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_positionsOld;

};
