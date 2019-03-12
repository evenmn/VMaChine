#pragma once
#include "wavefunction.h"

class PadeJastrow : public WaveFunction {
public:
    PadeJastrow(class System* system);
    void updateArrays(const Eigen::VectorXd positions, const int pRand);
    void resetArrays();
    void initializeArrays(const Eigen::VectorXd positions);
    void updateParameters(const Eigen::MatrixXd parameters, const int elementNumber);
    double evaluateRatio();
    double computeFirstDerivative(const int k);
    double computeSecondDerivative();
    Eigen::VectorXd computeFirstEnergyDerivative(const int k);
    Eigen::VectorXd computeSecondEnergyDerivative();

    Eigen::MatrixXd calculateDistanceMatrix(const Eigen::VectorXd particles);
    double calculateDistanceMatrixElement(const int i, const int j, const Eigen::VectorXd particles);
    void calculateDistanceMatrixCross(const int par, const Eigen::VectorXd particles, Eigen::MatrixXd &distanceMatrix);

    //double f(int i, int j);
    Eigen::MatrixXd f();
    double g(int i, int j, int k, int l);
    double beta(int i, int j);
    double gamma();

private:
    int m_elementNumber = 1;
    Eigen::MatrixXd m_distanceMatrix;
    Eigen::MatrixXd m_oldDistanceMatrix;
    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_oldPositions;
    Eigen::MatrixXd m_beta;
    Eigen::MatrixXd m_f;
    Eigen::MatrixXd m_oldF;
    Eigen::MatrixXd m_g;
    Eigen::MatrixXd m_oldG;
    double m_gamma;
    double m_oldRatio;
    double m_ratio;
};
