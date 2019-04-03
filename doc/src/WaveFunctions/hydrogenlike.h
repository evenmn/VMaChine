#pragma once
#include "wavefunction.h"

class HydrogenLike : public WaveFunction {
public:
    HydrogenLike(class System* system);
    void            updateArrays                (const Eigen::VectorXd positions, const int changedCoord);
    void            setArrays                   ();
    void            resetArrays                 ();
    void            initializeArrays            (const Eigen::VectorXd positions);
    void            updateParameters            (const Eigen::MatrixXd parameters, const int elementNumber);
    double          evaluateRatio               ();
    double          computeGradient             (const int k);
    double          computeLaplacian            ();
    Eigen::VectorXd computeParameterGradient    ();

    double          calculateRadialVectorElement(const int particle);
    void            calculateRadialVector       ();

private:
    double  m_alpha = 0;
    int     m_elementNumber = 0;
    int     m_Z = 1;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_positionsOld;
    Eigen::VectorXd m_radialVector;
    Eigen::VectorXd m_radialVectorOld;

    double m_probabilityRatio;
    double m_probabilityRatioOld;
};
