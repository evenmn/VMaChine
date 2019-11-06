#pragma once
#include "optimization.h"

class BarzilaiBorwein : public Optimization {
public:
    BarzilaiBorwein(System* system);
    arma::uword getNumberOfBatches() { return m_numberOfBatches; }
    std::string getLabel()           { return m_label; }
    arma::mat updateParameters();

private:
    int    m_step = 0;
    int    m_numberOfBatches = 1;
    double m_omega = 0;
    double m_omega_sqrd = 0;
    arma::mat m_oldParameters;
    arma::mat m_parameters;
    arma::mat m_oldGradients;
    arma::mat m_gradients;

    std::string     m_label = "BB";
};
