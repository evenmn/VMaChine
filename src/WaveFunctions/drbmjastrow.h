#pragma once
#include "wavefunction.h"

class DRBMJastrow : public WaveFunction {
public:
    DRBMJastrow(class System* system, int numberOfLayers);
    int             getNumberOfParameters       ()      { return m_numberOfParameters; }
    int             getGlobalArrayNeed          ()      { return m_globalArrayNeed; }
    std::string     getLabel                    ()      { return m_label; }

    void            updateParameters            (const Eigen::MatrixXd parameters);
    void            initializeArrays            (const Eigen::VectorXd positions, \
                                                 const Eigen::VectorXd radialVector, \
                                                 const Eigen::MatrixXd distanceMatrix);
    void            updateArrays                (const Eigen::VectorXd positions, \
                                                 const Eigen::VectorXd radialVector, \
                                                 const Eigen::MatrixXd distanceMatrix, \
                                                 const int changedCoord);
    void            setConstants                (const int elementNumber);
    void            setArrays                   ();
    void            resetArrays                 ();
    double          evaluateRatio               ();
    double          computeGradient             (const int k);
    double          computeLaplacian            ();
    Eigen::VectorXd computeParameterGradient    ();

    void            updateVectors               ();
    void            updateRatio                 ();
    void            updateGradient              ();
    void            updateLaplacian             ();

private:
    int    m_elementNumber          = 0;
    int    m_numberOfHiddenNodes    = 1;
    int    m_globalArrayNeed        = 0;
    int    m_numberOfParameters     = 0;
    int    m_numberOfLayers         = 0;
    double m_sigmaSqrd              = 1;

    double m_probabilityRatio       = 1;
    double m_probabilityRatioOld    = 1;

    Eigen::MatrixXd m_gradientPart;
    Eigen::MatrixXd m_gradientPartOld;
    Eigen::MatrixXd m_laplacianPart;
    Eigen::MatrixXd m_laplacianPartOld;
    Eigen::MatrixXd m_positionsPow;
    Eigen::MatrixXd m_positionsPowOld;
    Eigen::MatrixXd m_W;

    Eigen::VectorXd m_b;
    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_positionsOld;
    Eigen::VectorXd m_vOld;
    Eigen::VectorXd m_v;
    Eigen::VectorXd m_nOld;
    Eigen::VectorXd m_n;
    Eigen::VectorXd m_pOld;
    Eigen::VectorXd m_p;

    std::string m_label = "drbmjastrow";
};

