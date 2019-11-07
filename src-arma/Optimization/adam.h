#pragma once
#include "optimization.h"

class ADAM : public Optimization
{
public:
    ADAM(System *system);
    arma::uword getNumberOfBatches() { return m_numberOfBatches; }
    std::string getLabel() { return m_label; }

    void initialize();
    arma::mat updateParameters();

private:
    arma::uword m_numberOfBatches = 10;
    double m_beta1 = 0.9;
    double m_beta2 = 0.999;
    double m_epsilon = 1e-8;
    arma::uword m_step = 0;
    arma::mat m_g;
    arma::mat m_m;
    arma::mat m_v;
    arma::mat m_mHat;
    arma::mat m_vHat;
    arma::mat m_theta;
    arma::mat m_denom;

    std::string m_label = "ADAM";
};
