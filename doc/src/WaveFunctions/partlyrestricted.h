#pragma once
#include "wavefunction.h"

class PartlyRestricted : public WaveFunction {
public:
    PartlyRestricted(class System* system);
    void            updateArrays                (const Eigen::VectorXd positions, const int changedCoord);
    void            setArrays                   ();
    void            resetArrays                 ();
    void            initializeArrays            (const Eigen::VectorXd positions);
    void            updateParameters            (Eigen::MatrixXd parameters, const int elementNumber);
    double          evaluateRatio();
    double          computeGradient(const int k);
    double          computeLaplacian();
    Eigen::VectorXd computeParameterGradient();

    void            calculateProbabilityRatio(int changedCoord);

private:
    int     m_elementNumber = 0;
    double  m_xCxOld;
    double  m_xCx;
    double  m_probabilityRatio;
    double  m_probabilityRatioOld;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_positionsOld;
    Eigen::MatrixXd m_c;
};
