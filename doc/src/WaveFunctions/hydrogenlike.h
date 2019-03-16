#pragma once
#include "wavefunction.h"

class HydrogenLike : public WaveFunction {
public:
    HydrogenLike(class System* system);
    void updateArrays(const Eigen::VectorXd positions, const int pRand);
    void resetArrays(int pRand);
    void initializeArrays(const Eigen::VectorXd positions);
    void updateParameters(const Eigen::MatrixXd parameters, const int elementNumber);
    double evaluateRatio();
    double computeFirstDerivative(const int k);
    double computeSecondDerivative();
    Eigen::VectorXd computeFirstEnergyDerivative(const int k);
    Eigen::VectorXd computeSecondEnergyDerivative();

    double calculateRadialVectorElement(const Eigen::VectorXd particles, const int par);
    Eigen::VectorXd calculateRadialVector(const Eigen::VectorXd particles);

private:
    double m_alpha = 0;
    int m_elementNumber = 0;
    int m_Z = 1;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_oldPositions;
    Eigen::VectorXd m_radialVector;
    Eigen::VectorXd m_oldRadialVector;

    double m_ratio;
    double m_oldRatio;
};
