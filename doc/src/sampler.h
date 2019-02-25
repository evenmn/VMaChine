#pragma once
#include <Eigen/Dense>
#include <fstream>

class Sampler {
public:
    Sampler(class System* system);
    void sample(const bool acceptedStep, const int stepNumber);
    void printOutputToTerminal(int iter, const int maxIter, const double time);
    void openOutputFiles(const std::string path);
    void printOutputToFile();
    void closeOutputFiles();
    void printImmediatelyToFile(const Eigen::VectorXd positions);
    void computeAverages();
    double getEnergy()          { return m_energy; }
    Eigen::MatrixXd getdE()     { return m_dE; }
    Eigen::MatrixXd getdEE()    { return m_dEE; }

private:
    int     m_totalNumberOfSteps = 0;
    int     m_maxNumberOfParametersPerElement = 0;
    int     m_numberOfMetropolisSteps = 0;
    int     m_numberOfParticles = 0;
    int     m_numberOfDimensions = 0;
    int     m_numberOfElements = 0;
    int     m_acceptenceRatio = 0;
    double  m_energy = 0;
    double  m_variance = 0;
    double  m_cumulativeEnergy = 0;
    double  m_equilibriumFraction = 0;
    double  m_SqrdE = 0;
    double  m_instantEnergy = 0;
    Eigen::MatrixXd  m_dE;
    Eigen::MatrixXd  m_dEE;
    std::ofstream m_energyFile;
    std::ofstream m_cumulativeFile;
    class System* m_system = nullptr;

    // One-body density related stuff
    bool m_calculateOneBody = true;
    int m_numberOfBins = 500;
    double m_maxRadius = 5;
    double m_radialStep = m_maxRadius/m_numberOfBins;
    Eigen::VectorXd m_binLinSpace = Eigen::VectorXd::LinSpaced(m_numberOfBins, 0, m_maxRadius);
    Eigen::VectorXd m_particlesPerBin = Eigen::VectorXd::Zero(m_numberOfBins);
    std::ofstream m_oneBodyFile;
};
