#pragma once
#include "wavefunction.h"

class PadeJastrow : public WaveFunction
{
public:
    PadeJastrow(class System *system);
    arma::uword getNumberOfParameters() { return m_numberOfParameters; }
    arma::uword getGlobalArrayNeed() { return m_globalArrayNeed; }
    std::string getLabel() { return m_label; }

    void updateParameters(const arma::mat parameters);
    void initializeArrays(const arma::vec positions,
                          const arma::vec radialVector,
                          const arma::mat distanceMatrix);
    void updateArrays(const arma::vec positions,
                      const arma::vec radialVector,
                      const arma::mat distanceMatrix,
                      const arma::uword changedCoord);
    void setConstants(const arma::uword elementNumber);
    void setArrays();
    void resetArrays();
    double evaluateRatio();
    double computeGradient(const arma::uword k);
    double computeLaplacian();
    arma::vec computeParameterGradient();

    void updateMatrices(arma::uword particle);
    void calculateProbabilityRatio(arma::uword particle);
    void initializeBeta();
    void initializeMatrices();
    void initializePrincipalDistance();
    void updatePrincipalDistance(arma::uword i, arma::uword i_p);

private:
    arma::uword m_numberOfParameters = 1;
    arma::uword m_globalArrayNeed = 1;
    arma::uword m_elementNumber = 0;

    double m_gamma;
    double m_probabilityRatio;
    double m_probabilityRatioOld;

    arma::vec m_positions;
    arma::vec m_positionsOld;

    arma::mat m_distanceMatrix;
    arma::mat m_distanceMatrixOld;
    arma::mat m_beta;
    arma::mat m_f;
    arma::mat m_fOld;
    arma::mat m_h;
    arma::mat m_hOld;
    arma::mat m_hOldOld;
    arma::mat m_gradients;
    arma::mat m_principalDistance;
    arma::mat m_principalDistanceOld;

    std::string m_label = "padejastrow";
};

