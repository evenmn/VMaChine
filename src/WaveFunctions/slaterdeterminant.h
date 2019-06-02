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

    void            updateSlaterMatrixElementUp         (const int i, const int j);
    void            updateSlaterMatrixRowUp             (const int row);
    void            updateSlaterMatrixDerRowUp          (const int row);
    void            updateSlaterMatrixSecDerRowUp       (const int k);
    void            updateSlaterMatrixInverseUp         ();
    void            updateSlaterDeterminantDerivativesUp();

    void            updateSlaterMatrixElementDn         (const int i, const int j);
    void            updateSlaterMatrixRowDn             (const int row);
    void            updateSlaterMatrixDerRowDn          (const int row);
    void            updateSlaterMatrixSecDerRowDn       (const int k);
    void            updateSlaterMatrixInverseDn         ();
    void            updateSlaterDeterminantDerivativesDn();

    double          updateRatioUp                     ();
    double          updateRatioDn                     ();

    int             spinUp();
    int             spinDn();

private:
    int             m_numberOfParameters    = 0;
    int             m_globalArrayNeed       = 0;
    int             m_freeDimensionsHalf    = 0;
    int             m_particle              = 0;
    int             m_dimension             = 0;
    int             m_elementNumber         = 0;
    int             m_numberOfSpinUp        = 0;
    int             m_numberOfSpinDn        = 0;
    int             m_numberOfFreeDimUp     = 0;
    int             m_numberOfFreeDimDn     = 0;

    double          m_numberOfParticlesHalf = 0;
    double          m_totalSpin             = 0;


    Eigen::VectorXd m_positions;
    Eigen::MatrixXd m_positionBlock;
    Eigen::MatrixXd m_slaterMatrixUp;
    Eigen::MatrixXd m_slaterMatrixDerUp;
    Eigen::MatrixXd m_slaterMatrixSecDerUp;
    Eigen::MatrixXd m_slaterMatrixInverseUp;

    Eigen::MatrixXd m_slaterMatrixDn;
    Eigen::MatrixXd m_slaterMatrixDerDn;
    Eigen::MatrixXd m_slaterMatrixSecDerDn;
    Eigen::MatrixXd m_slaterMatrixInverseDn;

    Eigen::VectorXd m_determinantDerivative;
    Eigen::VectorXd m_determinantSecondDerivative;
    double          m_probabilityRatio = 0;

    Eigen::VectorXd m_positionsOld;
    Eigen::MatrixXd m_positionBlockOld;
    Eigen::MatrixXd m_slaterMatrixUpOld;
    Eigen::MatrixXd m_slaterMatrixDerUpOld;
    Eigen::MatrixXd m_slaterMatrixSecDerUpOld;
    Eigen::MatrixXd m_slaterMatrixInverseUpOld;

    Eigen::MatrixXd m_slaterMatrixDnOld;
    Eigen::MatrixXd m_slaterMatrixDerDnOld;
    Eigen::MatrixXd m_slaterMatrixSecDerDnOld;
    Eigen::MatrixXd m_slaterMatrixInverseDnOld;

    Eigen::VectorXd m_determinantDerivativeOld;
    Eigen::VectorXd m_determinantSecondDerivativeOld;
    double          m_probabilityRatioOld = 0;

    std::string     m_label = "slaterdeterminant";
};
