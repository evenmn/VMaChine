#pragma once
#include <Eigen/Dense>
#include <vector>

class Optimization {
public:
    Optimization(class System* system);
    virtual unsigned int    getNumberOfBatches()        = 0;
    virtual Eigen::MatrixXd updateParameters()          = 0;
    virtual                 ~Optimization()             = 0;

    Eigen::MatrixXd getEnergyGradient();

protected:
    unsigned int    m_numberOfFreeDimensions            = 0;
    unsigned int    m_numberOfStepsAfterEquilibrium     = 0;
    unsigned short  m_numberOfWaveFunctionElements      = 0;
    unsigned int    m_maxNumberOfParametersPerElement   = 0;

    double          m_eta                               = 1;

    Eigen::VectorXd m_positions;

    std::vector<class WaveFunction*> m_waveFunctionVector;

    class System* m_system = nullptr;
};
