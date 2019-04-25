#pragma once
#include "wavefunction.h"

class SlaterDeterminant : public WaveFunction {
public:
    SlaterDeterminant(class System* system);
    void            updateArrays                    (const Eigen::VectorXd positions, const int changedCoord);
    void            setArrays                       ();
    void            resetArrays                     ();
    void            initializeArrays                (const Eigen::VectorXd positions);
    void            updateParameters                (const Eigen::MatrixXd parameters, const int elementNumber);
    double          evaluateRatio                   ();
    double          computeGradient                 (const int k);
    double          computeLaplacian                ();
    Eigen::VectorXd computeParameterGradient        ();

    void            initializeSlaterMatrix          ();
    void            initializeSlaterMatrixDer       ();
    void            initializeSlaterMatrixSecDer    ();
    void            initializeSlaterMatrixInverse   ();

    void            updateSlaterMatrixElement       (const int i, const int j);
    void            updateSlaterMatrixRow           (const int i);
    void            updateSlaterMatrixDerRow        (const int k);
    void            updateSlaterMatrixSecDerRow     (const int k);
    void            updateSlaterMatrixInverse       (int start, int end);
    void            updateSlaterDeterminantDerivatives(int start, int end);
    double          updateRatio                     ();

private:
    int             m_elementNumber         = 0;
    int             m_numberOfParticlesHalf = 0;
    int             m_freeDimensionsHalf    = 0;
    int             m_particle              = 0;
    int             m_dimension             = 0;

    Eigen::MatrixXd m_listOfStates;

    Eigen::VectorXd m_positions;
    Eigen::MatrixXd m_slaterMatrix;
    Eigen::MatrixXd m_slaterMatrixDer;
    Eigen::MatrixXd m_slaterMatrixSecDer;
    Eigen::MatrixXd m_slaterMatrixInverse;
    Eigen::VectorXd m_determinantDerivative;
    Eigen::VectorXd m_determinantSecondDerivative;
    double          m_probabilityRatio = 0;

    Eigen::VectorXd m_positionsOld;
    Eigen::MatrixXd m_slaterMatrixOld;
    Eigen::MatrixXd m_slaterMatrixDerOld;
    Eigen::MatrixXd m_slaterMatrixSecDerOld;
    Eigen::MatrixXd m_slaterMatrixInverseOld;
    Eigen::VectorXd m_determinantDerivativeOld;
    Eigen::VectorXd m_determinantSecondDerivativeOld;
    double          m_probabilityRatioOld = 0;
};