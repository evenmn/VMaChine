#pragma once
#include "wavefunction.h"

class SamsethJastrow : public WaveFunction {
public:
    SamsethJastrow(class System* system);
    void            updateArrays                (const Eigen::VectorXd positions, const int changedCoord);
    void            setArrays                   ();
    void            resetArrays                 ();
    void            initializeArrays            (const Eigen::VectorXd positions);
    void            updateParameters            (const Eigen::MatrixXd parameters, const int elementNumber);
    double          evaluateRatio               ();
    double          computeGradient             (const int k);
    double          computeLaplacian            ();
    Eigen::VectorXd computeParameterGradient    ();

    void            calculateG(int pRand);
    void            calculateProbabilityRatio(int particle);
    void            calculateDistanceMatrix();
    double          calculateDistanceMatrixElement(const int i, const int j);
    void            calculateDistanceMatrixCross(const int par);

private:
    int m_elementNumber = 1;
    Eigen::MatrixXd m_distanceMatrix;
    Eigen::MatrixXd m_distanceMatrixOld;
    Eigen::MatrixXd m_distanceMatrixSqrd;
    Eigen::MatrixXd m_distanceMatrixSqrdOld;
    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_positionsOld;
    Eigen::MatrixXd m_g;
    Eigen::MatrixXd m_gOld;
    Eigen::MatrixXd m_gSqrd;
    Eigen::MatrixXd m_gSqrdOld;
    double m_beta;
    double m_betaSqrd;
    double m_gamma;
    double m_betagamma;
    double m_probabilityRatioOld;
    double m_probabilityRatio;
};
