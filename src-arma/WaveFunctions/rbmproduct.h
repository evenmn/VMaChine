#pragma once
#include "wavefunction.h"

class RBMProduct : public WaveFunction
{
public:
    RBMProduct(class System *system);
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

    void updateVectors();
    void updateRatio();

private:
    arma::uword m_numberOfParameters = 1;
    arma::uword m_numberOfHiddenNodes = 1;

    double m_sigmaSqrd = 1;
    double m_sigmaQuad = 1;
    double m_probabilityRatio = 1;
    double m_probabilityRatioOld = 1;

    arma::vec m_positions;
    arma::vec m_positionsOld;
    arma::vec m_b;
    arma::mat m_W;
    arma::mat m_WSqrd;
    arma::vec m_vOld;
    arma::vec m_v;
    arma::vec m_nOld;
    arma::vec m_n;
    arma::vec m_pOld;
    arma::vec m_p;
    arma::vec m_gradients;

    // Properties of the element (DO NOT TOUCH!)
    std::string m_label = "rbmproduct";
    arma::uword m_globalArrayNeed = 0;
    arma::uword m_elementNumber = 1;
};
