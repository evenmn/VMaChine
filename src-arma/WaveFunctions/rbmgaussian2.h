#pragma once
#include "wavefunction.h"

class RBMGaussian2 : public WaveFunction {
public:
    RBMGaussian2(class System* system);
    int             getNumberOfParameters       ()      { return m_numberOfParameters; }
    int             getGlobalArrayNeed          ()      { return m_globalArrayNeed; }
    std::string     getLabel                    ()      { return m_label; }

    void            updateArrays                    (const Eigen::VectorXd positions, \
                                                     const Eigen::VectorXd radialVector, \
                                                     const Eigen::MatrixXd distanceMatrix, \
                                                     const int changedCoord);
    void            setArrays                   ();
    void            resetArrays                 ();
    void            initializeArrays            (const Eigen::VectorXd positions, \
                                                 const Eigen::VectorXd radialVector, \
                                                 const Eigen::MatrixXd distanceMatrix);
    void            updateParameters            (const Eigen::MatrixXd parameters, const int elementNumber);
    double          evaluateRatio               ();
    double          computeGradient             (const int k);
    double          computeLaplacian            ();
    Eigen::VectorXd computeParameterGradient    ();

    void            calculateG(int changedCoord);

private:
    int     m_numberOfParameters       = 1;
    int     m_globalArrayNeed          = 1;
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

    std::string m_label = "rbmgaussian2";
};
