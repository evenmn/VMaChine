#pragma once
#include "wavefunction.h"

class PadeJastrow : public WaveFunction {
public:
    PadeJastrow(class System* system);
    int             getNumberOfParameters       ()      { return m_numberOfParameters; }
    int             getGlobalArrayNeed          ()      { return m_globalArrayNeed; }
    std::string     getLabel                    ()      { return m_label; }

    void            updateParameters            (const Eigen::MatrixXd parameters);
    void            initializeArrays            (const Eigen::VectorXd positions, \
                                                 const Eigen::VectorXd radialVector, \
                                                 const Eigen::MatrixXd distanceMatrix);
    void            updateArrays                (const Eigen::VectorXd positions, \
                                                 const Eigen::VectorXd radialVector, \
                                                 const Eigen::MatrixXd distanceMatrix, \
                                                 const int changedCoord);
    void            setConstants                (const int elementNumber);
    void            setArrays                   ();
    void            resetArrays                 ();
    double          evaluateRatio               ();
    double          computeGradient             (const int k);
    double          computeLaplacian            ();
    Eigen::VectorXd computeParameterGradient    ();


    void            calculateF(int particle);
    void            calculateG(int pRand);
    void            calculateH(int particle);
    void            calculateProbabilityRatio(int particle);
    void            initializeBeta();

private:
    int m_numberOfParameters       = 1;
    int m_globalArrayNeed          = 1;
    int m_elementNumber = 1;
    Eigen::MatrixXd m_distanceMatrix;
    Eigen::MatrixXd m_distanceMatrixOld;
    Eigen::MatrixXd m_distanceMatrixSqrd;
    Eigen::MatrixXd m_distanceMatrixSqrdOld;
    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_positionsOld;
    Eigen::MatrixXd m_beta;
    Eigen::MatrixXd m_f;
    Eigen::MatrixXd m_fOld;
    Eigen::MatrixXd m_fSqrd;
    Eigen::MatrixXd m_fSqrdOld;
    Eigen::MatrixXd m_fCube;
    Eigen::MatrixXd m_fCubeOld;
    Eigen::MatrixXd m_g;
    Eigen::MatrixXd m_gOld;
    Eigen::MatrixXd m_gSqrd;
    Eigen::MatrixXd m_gSqrdOld;
    Eigen::MatrixXd m_h;
    Eigen::MatrixXd m_hOld;
    Eigen::MatrixXd m_hOldOld;
    double m_gamma;
    double m_probabilityRatioOld;
    double m_probabilityRatio;

    std::string m_label = "padejastrow";
};
