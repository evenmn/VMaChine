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
    double          evaluate();
    double          evaluateSqrd();
    double          computeFirstDerivative(const int k);
    double          computeSecondDerivative();
    Eigen::VectorXd computeFirstEnergyDerivative(const int k);
    Eigen::VectorXd computeSecondEnergyDerivative();
private:
    int     m_elementNumber = 0;
    double  m_sigmaSqrd2    = 1;
    double  m_oldXCx;
    double  m_xCx;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_oldPositions;
    Eigen::MatrixXd m_c;
};
