#pragma once
#include "wavefunction.h"

class RBMJastrow5 : public WaveFunction {
public:
    RBMJastrow5(class System* system);
    void            updateArrays                (Eigen::VectorXd positions, const int changedCoord);
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
    void            updateGradient              ();

private:
    int    m_elementNumber          = 1;
    int    m_numberOfHiddenNodes    = 1;
    double m_sigmaSqrd              = 1;
    double m_sigmaQuad              = 1;

    double m_probabilityRatio       = 1;
    double m_probabilityRatioOld    = 1;

    Eigen::MatrixXd m_gradientPart;
    Eigen::MatrixXd m_gradientPartOld;
    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_positionsOld;
    Eigen::VectorXd m_positionsSqrd;
    Eigen::VectorXd m_positionsSqrdOld;
    Eigen::VectorXd m_b;
    Eigen::MatrixXd m_W1;
    Eigen::MatrixXd m_W1Sqrd;
    Eigen::MatrixXd m_W2;
    Eigen::MatrixXd m_W2Sqrd;
    Eigen::VectorXd m_vOld;
    Eigen::VectorXd m_v;
    Eigen::VectorXd m_nOld;
    Eigen::VectorXd m_n;
    Eigen::VectorXd m_pOld;
    Eigen::VectorXd m_p;
    Eigen::VectorXd m_pMinusN;
    Eigen::VectorXd m_pMinusNOld;
    Eigen::VectorXd m_pDotN;
    Eigen::VectorXd m_pDotNOld;
};

