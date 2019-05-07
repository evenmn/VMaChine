#pragma once
#include "wavefunction.h"

class RBMJastrow2 : public WaveFunction {
public:
    RBMJastrow2(class System* system);
    unsigned int    getNumberOfParameters       ()  { return m_numberOfParameters; }
    unsigned short  getGlobalArrayNeed          ()  { return m_globalArrayNeed; }
    std::string     getLabel                    ()  { return m_label; }

    void            updateArrays                (const Eigen::VectorXd positions, \
                                                 const Eigen::VectorXd radialVector, \
                                                 const Eigen::MatrixXd distanceMatrix, \
                                                 const unsigned int changedCoord);
    void            setArrays                   ();
    void            resetArrays                 ();
    void            initializeArrays            (const Eigen::VectorXd positions, \
                                                 const Eigen::VectorXd radialVector, \
                                                 const Eigen::MatrixXd distanceMatrix);
    void            updateParameters            (const Eigen::MatrixXd parameters, \
                                                 const unsigned short elementNumber);
    double          evaluateRatio               ();
    double          computeGradient             (const unsigned int k);
    double          computeLaplacian            ();
    Eigen::VectorXd computeParameterGradient    ();

    void            updateVectors               ();
    void            updateRatio                 ();

private:
    unsigned int    m_numberOfParameters    = 0;
    unsigned int    m_numberOfHiddenNodes   = 0;
    unsigned short  m_globalArrayNeed       = 0;
    unsigned short  m_elementNumber         = 0;

    double          m_sigmaSqrd             = 1;
    double          m_sigmaQuad             = 1;
    double          m_probabilityRatio      = 1;
    double          m_probabilityRatioOld   = 1;

    Eigen::MatrixXd m_distanceMatrix;
    Eigen::MatrixXd m_distanceMatrixOld;

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

    std::string m_label = "rbmjastrow2";
};
