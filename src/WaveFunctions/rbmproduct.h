#pragma once
#include "wavefunction.h"

class RBMProduct : public WaveFunction
{
public:
    RBMProduct(class System *system);
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
    Eigen::VectorXd m_b;
    Eigen::MatrixXd m_W;
    Eigen::MatrixXd m_WSqrd;
    Eigen::VectorXd m_vOld;
    Eigen::VectorXd m_v;
    Eigen::VectorXd m_nOld;
    Eigen::VectorXd m_n;
    Eigen::VectorXd m_pOld;
    Eigen::VectorXd m_p;
    Eigen::VectorXd m_pDotN;
    Eigen::VectorXd m_pDotNOld;

    // Properties of the element (DO NOT TOUCH!)
    std::string m_label = "rbmproduct";
    int m_globalArrayNeed = 0;
    int m_elementNumber = 1;
};
