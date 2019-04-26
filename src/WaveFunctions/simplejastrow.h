#pragma once
#include "wavefunction.h"

class SimpleJastrow : public WaveFunction {
public:
    SimpleJastrow(class System* system);
    int             getNumberOfParameters       ()      { return m_numberOfParameters; }
    std::string     getLabel                    ()      { return m_label; }
    void            updateArrays                (const Eigen::VectorXd positions, const int changedCoord);
    void            setArrays                   ();
    void            resetArrays                 ();
    void            initializeArrays            (const Eigen::VectorXd positions);
    void            updateParameters            (const Eigen::MatrixXd parameters, const int elementNumber);
    double          evaluateRatio               ();
    double          computeGradient             (const int k);
    double          computeLaplacian            ();
    Eigen::VectorXd computeParameterGradient    ();

    void            calculateDistanceMatrix();
    double          calculateDistanceMatrixElement(const int i, const int j);
    void            calculateDistanceMatrixCross(const int par);
    void            calculateG(int pRand);
    void            calculateProbabilityRatio(int particle);

private:
    int     m_numberOfParameters       = 1;
    int     m_elementNumber = 1;
    Eigen::MatrixXd m_distanceMatrix;
    Eigen::MatrixXd m_distanceMatrixOld;
    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_positionsOld;
    Eigen::MatrixXd m_beta;
    Eigen::MatrixXd m_g;
    double m_gamma;
    double m_probabilityRatioOld;
    double m_probabilityRatio;

    std::string m_label = "simplejastrow";
};
