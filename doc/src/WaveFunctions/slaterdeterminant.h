#pragma once
#include "wavefunction.h"

class SlaterDeterminant : public WaveFunction {
public:
    SlaterDeterminant(class System* system);
    double          updateSlaterMatrixElement    (const Eigen::VectorXd positions, const int i, const int j);
    Eigen::VectorXd updateSlaterMatrixRow        (const Eigen::VectorXd positions, const int i);
    Eigen::MatrixXd updateSlaterMatrix           (const Eigen::VectorXd positions);
    Eigen::VectorXd updateSlaterMatrixDerRow     (const Eigen::VectorXd positions, const int k);
    Eigen::MatrixXd updateSlaterMatrixDer        (const Eigen::VectorXd positions);

    void            resetArrays                  (int pRand);
    void            updateArrays                 (const Eigen::VectorXd positions, const int pRand);
    void            initializeArrays             (const Eigen::VectorXd positions);
    void            updateParameters             (const Eigen::MatrixXd parameters, const int elementNumber);
    double          evaluateRatio                ();
    double          computeFirstDerivative       (const int k);
    double          computeSecondDerivative      ();
    Eigen::VectorXd computeFirstEnergyDerivative (const int k);
    Eigen::VectorXd computeSecondEnergyDerivative();
private:
    int     m_elementNumber         = 0;
    int     m_numberOfOrbitals      = 0;
    int     m_numberOfParticlesHalf = 0;
    int     m_freeDimensionsHalf    = 0;

    Eigen::MatrixXd m_listOfStates;

    Eigen::VectorXd m_positions;
    Eigen::MatrixXd m_slaterMatrixUp;
    Eigen::MatrixXd m_slaterMatrixDn;
    Eigen::MatrixXd m_slaterMatrixUpDer;
    Eigen::MatrixXd m_slaterMatrixDnDer;
    Eigen::MatrixXd m_slaterMatrixUpInverse;
    Eigen::MatrixXd m_slaterMatrixDnInverse;
    Eigen::VectorXd m_determinantDerivative;
    double m_slaterDeterminantUp = 0;
    double m_slaterDeterminantDn = 0;
    double m_probabilityRatio = 0;

    Eigen::VectorXd m_positionsOld;
    Eigen::MatrixXd m_slaterMatrixUpOld;
    Eigen::MatrixXd m_slaterMatrixDnOld;
    Eigen::MatrixXd m_slaterMatrixUpDerOld;
    Eigen::MatrixXd m_slaterMatrixDnDerOld;
    Eigen::MatrixXd m_slaterMatrixUpInverseOld;
    Eigen::MatrixXd m_slaterMatrixDnInverseOld;
    Eigen::VectorXd m_determinantDerivativeOld;
    double m_slaterDeterminantUpOld = 0;
    double m_slaterDeterminantDnOld = 0;
    double m_probabilityRatioOld = 0;

    double m_slaterDeterminantUpOldOld = 0;
    double m_slaterDeterminantDnOldOld = 0;
};
