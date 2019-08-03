#pragma once
#include "wavefunction.h"

class SimpleJastrow : public WaveFunction
{
public:
    SimpleJastrow(class System *system);
    int getNumberOfParameters() { return m_numberOfParameters; }
    int getGlobalArrayNeed() { return m_globalArrayNeed; }
    std::string getLabel() { return m_label; }

    void updateParameters(const Eigen::MatrixXd parameters);
    void initializeArrays(const Eigen::VectorXd positions,
                          const Eigen::VectorXd radialVector,
                          const Eigen::MatrixXd distanceMatrix);
    void updateArrays(const Eigen::VectorXd positions,
                      const Eigen::VectorXd radialVector,
                      const Eigen::MatrixXd distanceMatrix,
                      const int changedCoord);
    void setConstants(const int elementNumber);
    void setArrays();
    void resetArrays();
    double evaluateRatio();
    double computeGradient(const int k);
    double computeLaplacian();
    Eigen::VectorXd computeParameterGradient();

    void initializePrincipalDistance();
    void updatePrincipalDistance(int i);
    void calculateProbabilityRatio(int particle);

private:
    int m_numberOfParameters = 1;
    int m_globalArrayNeed = 1;
    int m_elementNumber = 0;
    int m_particle = 0;

    double m_gamma;
    double m_probabilityRatioOld;
    double m_probabilityRatio;

    Eigen::MatrixXd m_beta;
    Eigen::VectorXd m_gradients;
    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_positionsOld;
    Eigen::MatrixXd m_distanceMatrix;
    Eigen::MatrixXd m_distanceMatrixOld;
    Eigen::MatrixXd m_principalDistance;
    Eigen::MatrixXd m_principalDistanceOld;

    std::string m_label = "simplejastrow";
};
