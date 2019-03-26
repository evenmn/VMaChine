#pragma once
#include "wavefunction.h"

class SlaterDeterminant : public WaveFunction {
public:
    SlaterDeterminant(class System* system);
    double          updateSlaterMatrixElement    (const Eigen::VectorXd positions, const int i, const int j);
    Eigen::VectorXd updateSlaterMatrixRow        (const Eigen::VectorXd positions, const int i);
    void initializeSlaterMatrix();
    Eigen::VectorXd updateSlaterMatrixDerRow     (const Eigen::VectorXd positions, const int k);
    void initializeSlaterMatrixDer();
    Eigen::VectorXd updateSlaterMatrixSecDerRow  (const Eigen::VectorXd positions, const int k);
    Eigen::MatrixXd initializeSlaterMatrixSecDer (const Eigen::VectorXd positions);

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
    Eigen::MatrixXd m_slaterMatrixUpSecDer;
    Eigen::MatrixXd m_slaterMatrixDnSecDer;
    Eigen::MatrixXd m_slaterMatrixUpInverse;
    Eigen::MatrixXd m_slaterMatrixDnInverse;
    Eigen::VectorXd m_determinantDerivative;
    Eigen::VectorXd m_determinantSecondDerivative;
    double m_probabilityRatio = 0;

    Eigen::VectorXd m_positionsOld;
    Eigen::MatrixXd m_slaterMatrixUpOld;
    Eigen::MatrixXd m_slaterMatrixDnOld;
    Eigen::MatrixXd m_slaterMatrixUpDerOld;
    Eigen::MatrixXd m_slaterMatrixDnDerOld;
    Eigen::MatrixXd m_slaterMatrixUpSecDerOld;
    Eigen::MatrixXd m_slaterMatrixDnSecDerOld;
    Eigen::MatrixXd m_slaterMatrixUpInverseOld;
    Eigen::MatrixXd m_slaterMatrixDnInverseOld;
    Eigen::VectorXd m_determinantDerivativeOld;
    Eigen::VectorXd m_determinantSecondDerivativeOld;
    double m_probabilityRatioOld = 0;

    Eigen::MatrixXd m_slaterMatrix;
    Eigen::MatrixXd m_slaterMatrixOld;
    Eigen::MatrixXd m_slaterMatrixInverse;
    Eigen::MatrixXd m_slaterMatrixInverseOld;
    Eigen::MatrixXd m_slaterMatrixDer;
    Eigen::MatrixXd m_slaterMatrixDerOld;
};
