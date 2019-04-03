#pragma once
#include "wavefunction.h"

class NQSJastrow2 : public WaveFunction {
public:
    NQSJastrow2(class System* system);
    void            updateArrays                (const Eigen::VectorXd positions, const int changedCoord);
    void            setArrays                   ();
    void            resetArrays                 ();
    void            initializeArrays            (const Eigen::VectorXd positions);
    void            updateParameters            (const Eigen::MatrixXd parameters, const int elementNumber);
    double          evaluateRatio               ();
    double          computeGradient             (const int k);
    double          computeLaplacian            ();
    Eigen::VectorXd computeParameterGradient    ();

    void            updateVectors               ();
    void            updateRatio                 ();

private:
    int    m_elementNumber          = 1;
    int    m_numberOfHiddenNodes    = 1;
    double m_sigmaSqrd              = 1;
    double m_sigmaQuad              = 1;

    double m_probabilityRatio       = 1;
    double m_probabilityRatioOld    = 1;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_positionsOld;
    Eigen::VectorXd m_b1;
    Eigen::MatrixXd m_W1;
    Eigen::VectorXd m_b2;
    Eigen::MatrixXd m_W2;
    Eigen::VectorXd m_gOld;
    Eigen::VectorXd m_g;
    Eigen::VectorXd m_hOld;
    Eigen::VectorXd m_h;
    Eigen::VectorXd m_nOld;
    Eigen::VectorXd m_n;
    Eigen::VectorXd m_pOld;
    Eigen::VectorXd m_p;
};
