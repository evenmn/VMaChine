#pragma once
#include "wavefunction.h"

class PadeJastrow : public WaveFunction {
public:
    PadeJastrow(class System* system);
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


    void            calculateF                  (const unsigned int particle);
    void            calculateG                  (const unsigned int changedCoord);
    void            calculateH                  (const unsigned int particle);
    void            calculateProbabilityRatio   (const unsigned int particle);
    void            initializeBeta              ();

private:
    unsigned int    m_numberOfParameters    = 1;
    unsigned short  m_globalArrayNeed       = 1;
    unsigned short  m_elementNumber         = 0;

    double          m_gamma                 = 1;
    double          m_probabilityRatioOld   = 1;
    double          m_probabilityRatio      = 1;

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

    std::string     m_label = "padejastrow";
};
