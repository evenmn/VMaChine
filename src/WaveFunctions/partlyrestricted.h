#pragma once
#include "wavefunction.h"

class PartlyRestricted : public WaveFunction
{
public:
    PartlyRestricted(class System *system);
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

    void calculateProbabilityRatio(int changedCoord);

private:
    int m_numberOfParameters = 1;
    int m_globalArrayNeed = 0;
    int m_elementNumber = 0;
    double m_xCxOld;
    double m_xCx;
    double m_probabilityRatio;
    double m_probabilityRatioOld;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_positionsOld;
    Eigen::MatrixXd m_c;
    Eigen::VectorXd m_gradients;

    std::string m_label = "partlyrestricted";
};
