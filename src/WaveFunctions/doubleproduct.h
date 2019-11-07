#pragma once
#include "wavefunction.h"

class DoubleProduct : public WaveFunction
{
public:
    DoubleProduct(class System *system);
    int getNumberOfParameters() { return m_numberOfParameters; }
    int getGlobalArrayNeed() { return m_globalArrayNeed; }
    std::string getLabel() { return m_label; }

    void updateParameters(const Eigen::MatrixXd parameters);
    void initializeArrays(const Eigen::VectorXd positions,
                          const Eigen::VectorXd radialVector,
                          const Eigen::MatrixXd distanceMatrix);
    void updateArrays(const Eigen::VectorXd positions,
                      const Eigen::VectorXd,
                      const Eigen::MatrixXd,
                      const int);
    void setConstants(const int elementNumber);
    void setArrays();
    void resetArrays();
    double evaluateRatio();
    double computeGradient(const int k);
    double computeLaplacian();
    Eigen::VectorXd computeParameterGradient();

    void updateVectors();
    void updateRatio();

private:
    int m_numberOfParameters = 1;
    int m_numberOfHiddenNodes = 1;

    double m_sigmaSqrd = 1;
    double m_sigmaQuad = 1;
    double m_probabilityRatio = 1;
    double m_probabilityRatioOld = 1;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_positionsOld;
    Eigen::MatrixXd m_positionBlock;
    Eigen::MatrixXd m_positionBlockOld;
    Eigen::VectorXd m_b;
    Eigen::MatrixXd m_W;
    Eigen::MatrixXd m_WSqrd;
    Eigen::MatrixXd m_vOld;
    Eigen::MatrixXd m_v;
    Eigen::MatrixXd m_nOld;
    Eigen::MatrixXd m_n;
    Eigen::MatrixXd m_pOld;
    Eigen::MatrixXd m_p;
    Eigen::VectorXd m_gradients;

    // Properties of the element (DO NOT TOUCH!)
    std::string m_label = "rbmproduct";
    int m_globalArrayNeed = 0;
    int m_elementNumber = 1;
};
