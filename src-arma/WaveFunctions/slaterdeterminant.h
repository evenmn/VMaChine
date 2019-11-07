#pragma once
#include "wavefunction.h"

class SlaterDeterminant : public WaveFunction
{
public:
    SlaterDeterminant(class System *system);
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

    void initializeSlaterMatrix();
    void initializeSlaterMatrixDer();
    void initializeSlaterMatrixSecDer();
    void initializeSlaterMatrixInverse();

    void updateSlaterMatrixElement(const arma::uword i, const arma::uword j);
    void updateSlaterMatrixRow(const arma::uword row);
    void updateSlaterMatrixDerRow(const arma::uword row);
    void updateSlaterMatrixSecDerRow(const arma::uword k);
    void updateSlaterMatrixInverse(arma::uword start, arma::uword end);
    void updateSlaterDeterminantDerivatives(arma::uword start, arma::uword end);
    void updateRatio();

private:
    arma::uword m_numberOfParameters = 0;
    arma::uword m_globalArrayNeed = 0;
    arma::uword m_numberOfParticlesHalf = 0;
    arma::uword m_freeDimensionsHalf = 0;
    arma::uword m_particle = 0;
    arma::uword m_dimension = 0;
    arma::uword m_elementNumber = 0;

    double m_ratio = 1;
    double m_ratioOld = 1;

    arma::mat m_positions;
    arma::mat m_slaterMatrix;
    arma::mat m_slaterMatrixDer;
    arma::mat m_slaterMatrixSecDer;
    arma::mat m_slaterMatrixInverse;
    arma::vec m_determinantDerivative;
    arma::vec m_determinantSecondDerivative;
    double m_probabilityRatio = 0;

    arma::mat m_positionsOld;
    arma::mat m_slaterMatrixOld;
    arma::mat m_slaterMatrixDerOld;
    arma::mat m_slaterMatrixSecDerOld;
    arma::mat m_slaterMatrixInverseOld;
    arma::vec m_determinantDerivativeOld;
    arma::vec m_determinantSecondDerivativeOld;
    double m_probabilityRatioOld = 0;

    arma::vec m_gradients;
    std::string m_label = "slaterdeterminant";

    class Basis *m_basis = nullptr;
};
