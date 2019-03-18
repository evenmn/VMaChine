#pragma once
#include <Eigen/Dense>
#include <fstream>

class Sampler {
public:
    Sampler(class System* system);
    void            sample(const bool acceptedStep, const int stepNumber);
    void            printOutputToTerminal(const int maxIter, const double time);
    void            openOutputFiles(const std::string path);
    void            printOutputToFile();
    void            closeOutputFiles();
    void            printInstantValuesToFile(const Eigen::VectorXd positions);
    void            computeAverages();
    double          getAverageEnergy()        { return m_averageEnergy; }
    Eigen::MatrixXd getAverageGradients()     { return m_averageGradients; }
    Eigen::MatrixXd getAverageGradientsE()    { return m_averageGradientsE; }
    std::string     generateFileName(const std::string name, const std::string extension);

private:
    int     m_totalNumberOfSteps = 0;
    int     m_maxNumberOfParametersPerElement = 0;
    int     m_numberOfMetropolisSteps = 0;
    int     m_numberOfParticles = 0;
    int     m_numberOfDimensions = 0;
    int     m_numberOfElements = 0;
    int     m_numberOfBatches = 0;
    int     m_numberOfStepsPerBatch = 0;
    int     m_acceptenceRatio = 0;
    int     m_iter = 0;
    bool    m_interaction = 0;
    double  m_variance = 0;
    double  m_equilibriumFraction = 0;
    double  m_omega = 0;

    double  m_averageEnergy = 0;
    double  m_averageEnergySqrd = 0;
    double  m_cumulativeEnergy = 0;
    double  m_cumulativeEnergySqrd = 0;
    double  m_instantEnergy = 0;

    Eigen::MatrixXd  m_averageGradients;
    Eigen::MatrixXd  m_averageGradientsE;
    Eigen::MatrixXd  m_cumulativeGradients;
    Eigen::MatrixXd  m_cumulativeGradientsE;
    Eigen::MatrixXd  m_instantGradients;

    std::ofstream m_averageEnergyFile;
    std::ofstream m_instantEnergyFile;
    std::string   m_filename = "Filename not generated yet";
    class System* m_system = nullptr;

    // One-body density related stuff
    bool            m_calculateOneBody = true;
    int             m_numberOfBins = 500;
    double          m_maxRadius = 5;
    double          m_radialStep = m_maxRadius/m_numberOfBins;
    Eigen::VectorXd m_binLinSpace = Eigen::VectorXd::LinSpaced(m_numberOfBins, 0, m_maxRadius);
    Eigen::VectorXd m_particlesPerBin = Eigen::VectorXd::Zero(m_numberOfBins);
    std::ofstream   m_oneBodyFile;
};
