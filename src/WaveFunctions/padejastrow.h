#pragma once
#include "wavefunction.h"

class PadeJastrow : public WaveFunction
{
public:
    PadeJastrow(class System *system);
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

    void updateMatrices(int particle);
    void calculateG(int pRand);
    void calculateProbabilityRatio(int particle);
    void initializeBeta();
    void initializeMatrices();
    void initializePrincipalDistance();
    void updatePrincipalDistance(int i, int i_p);

private:
    int m_numberOfParameters = 1;
    int m_globalArrayNeed = 1;
    int m_elementNumber = 0;

    double m_gamma;
    double m_probabilityRatio;
    double m_probabilityRatioOld;

    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_positionsOld;

    Eigen::MatrixXd m_distanceMatrix;
    Eigen::MatrixXd m_distanceMatrixOld;
    Eigen::MatrixXd m_beta;
    Eigen::MatrixXd m_f;
    Eigen::MatrixXd m_fOld;
    Eigen::MatrixXd m_g;
    Eigen::MatrixXd m_gOld;
    Eigen::MatrixXd m_h;
    Eigen::MatrixXd m_hOld;
    Eigen::MatrixXd m_hOldOld;
    Eigen::MatrixXd m_gradients;
    Eigen::MatrixXd m_principalDistance;
    Eigen::MatrixXd m_principalDistanceOld;

    std::string m_label = "padejastrow";
};

