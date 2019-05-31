#pragma once
#include <Eigen/Dense>
#include <vector>

class Optimization {
public:
    Optimization(class System* system);
    virtual int  getNumberOfBatches() = 0;
    virtual std::string     getLabel() = 0;
    virtual Eigen::MatrixXd updateParameters() = 0;
    virtual ~Optimization() = 0;

    Eigen::MatrixXd getEnergyGradient();

protected:
    class System* m_system = nullptr;
    Eigen::VectorXd m_positions;
    int m_numberOfFreeDimensions = 0;
    int m_numberOfStepsAfterEquilibrium = 0;
    int m_numberOfElements = 0;
    int m_maxParameters = 0;
    double m_eta = 1;
    std::vector<class WaveFunction*> m_waveFunctionVector;
};
