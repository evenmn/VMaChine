#pragma once
#include "wavefunction.h"

class RBMGaussian2 : public WaveFunction {
public:
    RBMGaussian2(class System* system);
    void            updateArrays                (const Eigen::VectorXd positions, const int changedCoord);
    void            setArrays                   ();
    void            resetArrays                 ();
    void            initializeArrays            (const Eigen::VectorXd positions);
    void            updateParameters            (const Eigen::MatrixXd parameters, const int elementNumber);
    double          evaluateRatio               ();
    double          computeGradient             (const int k);
    double          computeLaplacian            ();
    Eigen::VectorXd computeParameterGradient    ();

    void            calculateDistanceMatrix();
    double          calculateDistanceMatrixElement(const int i, const int j);
    void            calculateDistanceMatrixCross(const int par);
    void            calculateG(int changedCoord);

private:
    int     m_elementNumber = 0;
    double  m_omega         = 1;
    double  m_sigmaSqrd     = 1;

    double  m_probabilityRatio = 0;
    double  m_probabilityRatioOld = 0;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_positionsOld;
    Eigen::MatrixXd m_distanceMatrix;
    Eigen::MatrixXd m_distanceMatrixOld;
    Eigen::MatrixXd m_Xa;
    Eigen::MatrixXd m_XaOld;
    Eigen::MatrixXd m_a;
    Eigen::MatrixXd m_g;
    Eigen::MatrixXd m_gOld;
    Eigen::MatrixXd m_gSqrd;
    Eigen::MatrixXd m_gSqrdOld;
};
