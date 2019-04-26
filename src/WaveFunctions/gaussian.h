#pragma once
#include "wavefunction.h"

class Gaussian : public WaveFunction {
public:
    Gaussian(class System* system);
    int             getNumberOfParameters       ()      { return m_numberOfParameters; }
    std::string     getLabel                    ()      { return m_label; }
    void            updateArrays                (const Eigen::VectorXd positions, const int changedCoord);
    void            setArrays                   ();
    void            resetArrays                 ();
    void            initializeArrays            (const Eigen::VectorXd positions);
    void            updateParameters            (const Eigen::MatrixXd parameters, const int elementNumber);
    double          evaluateRatio               ();
    double          computeGradient             (const int k);
    double          computeLaplacian            ();
    Eigen::VectorXd computeParameterGradient    ();

    void            updateProbabilityRatio      (int changedCoord);

private:
    int     m_numberOfParameters       = 1;

    int     m_elementNumber            = 0;
    double  m_omega                    = 0;
    double  m_alpha                    = 0;
    double  m_probabilityRatio         = 0;
    double  m_probabilityRatioOld      = 0;
    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_positionsOld;

    std::string m_label = "gaussian";
};
