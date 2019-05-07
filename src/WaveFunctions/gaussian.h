#pragma once
#include "wavefunction.h"

class Gaussian : public WaveFunction {
public:
    Gaussian(class System* system);
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

    void            updateProbabilityRatio      (const unsigned int changedCoord);

private:
    unsigned int    m_numberOfParameters    = 1;
    unsigned short  m_globalArrayNeed       = 0;
    unsigned short  m_elementNumber         = 0;

    double          m_omega                 = 1;
    double          m_alpha                 = 1;
    double          m_probabilityRatio      = 1;
    double          m_probabilityRatioOld   = 1;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_positionsOld;

    std::string     m_label = "gaussian";
};
