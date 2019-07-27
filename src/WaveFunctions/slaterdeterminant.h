#pragma once
#include "wavefunction.h"

class SlaterDeterminant : public WaveFunction {
public:
    SlaterDeterminant(class System* system);
    int             getNumberOfParameters           ()      { return m_numberOfParameters; }
    int             getGlobalArrayNeed              ()      { return m_globalArrayNeed; }
    std::string     getLabel                        ()      { return m_label; }

    void            updateParameters                (const Eigen::MatrixXd parameters);
    void            initializeArrays                (const Eigen::VectorXd positions, \
                                                     const Eigen::VectorXd radialVector, \
                                                     const Eigen::MatrixXd distanceMatrix);
    void            updateArrays                    (const Eigen::VectorXd positions, \
                                                     const Eigen::VectorXd radialVector, \
                                                     const Eigen::MatrixXd distanceMatrix, \
                                                     const int changedCoord);
    void            setConstants                    (const int elementNumber);
    void            setArrays                       ();
    void            resetArrays                     ();
    double          evaluateRatio                   ();
    double          computeGradient                 (const int k);
    double          computeLaplacian                ();
    Eigen::VectorXd computeParameterGradient        ();

    void            initializeSlaterMatrix          ();
    void            initializeSlaterMatrixDer       ();
    void            initializeSlaterMatrixSecDer    ();
    void            initializeSlaterMatrixInverse   ();

    void            updateSlaterMatrixElement       (const int i, const int j);
    void            updateSlaterMatrixRow           (const int row);
    void            updateSlaterMatrixDerRow        (const int row);
    void            updateSlaterMatrixSecDerRow     (const int k);
    void            updateSlaterMatrixInverse       (int start, int end);
    void            updateSlaterDeterminantDerivatives(int start, int end);
    void            updateRatio                     ();

private:
    int             m_numberOfParameters    = 0;
    int             m_globalArrayNeed       = 0;
    int             m_numberOfParticlesHalf = 0;
    int             m_freeDimensionsHalf    = 0;
    int             m_particle              = 0;
    int             m_dimension             = 0;
    int             m_elementNumber         = 0;

    double          m_ratio = 1;
    double          m_ratioOld = 1;

    Eigen::VectorXd m_positions;
    Eigen::MatrixXd m_positionBlock;
    Eigen::MatrixXd m_slaterMatrix;
    Eigen::MatrixXd m_slaterMatrixDer;
    Eigen::MatrixXd m_slaterMatrixSecDer;
    Eigen::MatrixXd m_slaterMatrixInverse;
    Eigen::VectorXd m_determinantDerivative;
    Eigen::VectorXd m_determinantSecondDerivative;
    double          m_probabilityRatio = 0;

    Eigen::VectorXd m_positionsOld;
    Eigen::MatrixXd m_positionBlockOld;
    Eigen::MatrixXd m_slaterMatrixOld;
    Eigen::MatrixXd m_slaterMatrixDerOld;
    Eigen::MatrixXd m_slaterMatrixSecDerOld;
    Eigen::MatrixXd m_slaterMatrixInverseOld;
    Eigen::VectorXd m_determinantDerivativeOld;
    Eigen::VectorXd m_determinantSecondDerivativeOld;
    double          m_probabilityRatioOld = 0;

    std::string     m_label = "slaterdeterminant";
};
