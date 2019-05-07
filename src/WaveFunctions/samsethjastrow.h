#pragma once
#include "wavefunction.h"

class SamsethJastrow : public WaveFunction {
public:
    SamsethJastrow(class System* system);
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

    void            calculateG                  (const unsigned int changedCoord);
    void            calculateProbabilityRatio   (const unsigned int particle);

private:
    unsigned int    m_numberOfParameters    = 2;
    unsigned short  m_globalArrayNeed       = 0;
    unsigned short  m_elementNumber         = 0;

    double          m_beta                  = 0;
    double          m_betaSqrd              = 0;
    double          m_gamma                 = 0;
    double          m_betagamma             = 0;
    double          m_probabilityRatioOld   = 0;
    double          m_probabilityRatio      = 0;

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

    std::string m_label = "samsethjastrow";
};
