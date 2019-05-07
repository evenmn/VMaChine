#pragma once
#include "wavefunction.h"

class RBMGaussian2 : public WaveFunction {
public:
    RBMGaussian2(class System* system);
    unsigned int    getNumberOfParameters       ()  { return m_numberOfParameters; }
    unsigned short  getGlobalArrayNeed          ()  { return m_globalArrayNeed; }
    std::string     getLabel                    ()  { return m_label; }

    void            updateArrays                (const Eigen::VectorXd positions,      \
                                                 const Eigen::VectorXd radialVector,   \
                                                 const Eigen::MatrixXd distanceMatrix, \
                                                 const unsigned int changedCoord);
    void            setArrays                   ();
    void            resetArrays                 ();
    void            initializeArrays            (const Eigen::VectorXd positions,      \
                                                 const Eigen::VectorXd radialVector,   \
                                                 const Eigen::MatrixXd distanceMatrix);
    void            updateParameters            (const Eigen::MatrixXd parameters,     \
                                                 const unsigned short elementNumber);
    double          evaluateRatio               ();
    double          computeGradient             (const unsigned int k);
    double          computeLaplacian            ();
    Eigen::VectorXd computeParameterGradient    ();

    void            calculateG                  (unsigned int changedCoord);

private:
    unsigned int    m_numberOfHiddenNodes   = 0;
    unsigned int    m_numberOfParameters    = 0;
    unsigned short  m_globalArrayNeed       = 1;
    unsigned short  m_elementNumber         = 0;

    double          m_omega                 = 1;
    double          m_sigmaSqrd             = 1;
    double          m_sigmaQuad             = 1;
    double          m_probabilityRatio      = 1;
    double          m_probabilityRatioOld   = 1;

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

    std::string     m_label = "rbmgaussian2";
};
