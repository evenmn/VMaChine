#pragma once
#include "wavefunction.h"

class SlaterDeterminant : public WaveFunction {
public:
    SlaterDeterminant(class System* system);
    unsigned int    getNumberOfParameters           ()  { return m_numberOfParameters; }
    unsigned short  getGlobalArrayNeed              ()  { return m_globalArrayNeed; }
    std::string     getLabel                        ()  { return m_label; }

    void            updateArrays                    (const Eigen::VectorXd positions, \
                                                     const Eigen::VectorXd radialVector, \
                                                     const Eigen::MatrixXd distanceMatrix, \
                                                     const unsigned int changedCoord);
    void            setArrays                       ();
    void            resetArrays                     ();
    void            initializeArrays                (const Eigen::VectorXd positions, \
                                                     const Eigen::VectorXd radialVector, \
                                                     const Eigen::MatrixXd distanceMatrix);
    void            updateParameters                (const Eigen::MatrixXd parameters, \
                                                     const unsigned short elementNumber);
    double          evaluateRatio                   ();
    double          computeGradient                 (const unsigned int k);
    double          computeLaplacian                ();
    Eigen::VectorXd computeParameterGradient        ();

    void            initializeSlaterMatrix          ();
    void            initializeSlaterMatrixDer       ();
    void            initializeSlaterMatrixSecDer    ();
    void            initializeSlaterMatrixInverse   ();

    void            updateSlaterMatrixElement       (const unsigned int row, const unsigned int col);
    void            updateSlaterMatrixRow           (const unsigned int row);
    void            updateSlaterMatrixDerRow        (const unsigned int row);
    void            updateSlaterMatrixSecDerRow     (const unsigned int k);
    void            updateSlaterMatrixInverse       (const unsigned int start, const unsigned int end);
    void            updateSlaterDeterminantDer      (const unsigned int start, const unsigned int end);
    double          updateRatio                     ();

private:
    unsigned int    m_numberOfParticlesHalf = 0;
    unsigned int    m_freeDimensionsHalf    = 0;
    unsigned int    m_particle              = 0;
    unsigned int    m_particleHalf          = 0;
    unsigned int    m_numberOfParameters    = 0;
    unsigned short  m_globalArrayNeed       = 0;
    unsigned short  m_elementNumber         = 0;
    unsigned short  m_dimension             = 0;

    double          m_probabilityRatio      = 0;
    double          m_probabilityRatioOld   = 0;

    Eigen::VectorXd m_positions;
    Eigen::MatrixXd m_positionBlock;
    Eigen::MatrixXd m_slaterMatrix;
    Eigen::MatrixXd m_slaterMatrixDer;
    Eigen::MatrixXd m_slaterMatrixSecDer;
    Eigen::MatrixXd m_slaterMatrixInverse;
    Eigen::VectorXd m_determinantDerivative;
    Eigen::VectorXd m_determinantSecondDerivative;

    Eigen::VectorXd m_positionsOld;
    Eigen::MatrixXd m_positionBlockOld;
    Eigen::MatrixXd m_slaterMatrixOld;
    Eigen::MatrixXd m_slaterMatrixDerOld;
    Eigen::MatrixXd m_slaterMatrixSecDerOld;
    Eigen::MatrixXd m_slaterMatrixInverseOld;
    Eigen::VectorXd m_determinantDerivativeOld;
    Eigen::VectorXd m_determinantSecondDerivativeOld;

    std::string     m_label = "slaterdeterminant";
};
