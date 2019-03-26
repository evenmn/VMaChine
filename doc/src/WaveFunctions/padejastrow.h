#pragma once
#include "wavefunction.h"

class PadeJastrow : public WaveFunction {
public:
    PadeJastrow(class System* system);
    void            updateArrays(const Eigen::VectorXd positions, const int pRand);
    void            resetArrays();
    void            initializeArrays(const Eigen::VectorXd positions);
    void            updateParameters(const Eigen::MatrixXd parameters, const int elementNumber);
    double          evaluateRatio();
    double          computeFirstDerivative(const int k);
    double          computeSecondDerivative();
    Eigen::VectorXd computeFirstEnergyDerivative(const int k);
    Eigen::VectorXd computeSecondEnergyDerivative();

    void            setArrays();
    void            calculateF(int particle);
    void            calculateG(int pRand);
    void            calculateH(int particle);
    void            calculateProbabilityRatio(int particle);
    void            initializeBeta();

    void            calculateDistanceMatrix();
    double          calculateDistanceMatrixElement(const int i, const int j);
    void            calculateDistanceMatrixCross(const int par);

private:
    int m_elementNumber = 1;
    Eigen::MatrixXd m_distanceMatrix;
    Eigen::MatrixXd m_distanceMatrixOld;
    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_positionsOld;
    Eigen::MatrixXd m_beta;
    Eigen::MatrixXd m_f;
    Eigen::MatrixXd m_fOld;
    Eigen::MatrixXd m_g;
    Eigen::MatrixXd m_gOld;
    Eigen::MatrixXd m_h;
    Eigen::MatrixXd m_hOld;
    Eigen::MatrixXd m_hOldOld;
    double m_gamma;
    double m_probabilityRatioOld;
    double m_probabilityRatio;
};
