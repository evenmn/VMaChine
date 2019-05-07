#pragma once
#include "wavefunction.h"

class PartlyRestricted : public WaveFunction {
public:
    PartlyRestricted(class System* system);
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

    void            calculateProbabilityRatio   (const unsigned int changedCoord);

private:
    unsigned int    m_numberOfParameters    = 0;
    unsigned short  m_globalArrayNeed       = 0;
    unsigned short  m_elementNumber         = 0;

    double          m_xCxOld                = 1;
    double          m_xCx                   = 1;
    double          m_probabilityRatio      = 1;
    double          m_probabilityRatioOld   = 1;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_positionsOld;
    Eigen::MatrixXd m_c;

    std::string     m_label = "partlyrestricted";
};
