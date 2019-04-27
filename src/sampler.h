#pragma once
#include <Eigen/Dense>
#include <fstream>

class Sampler {
public:
    Sampler(class System* system);
    void            sample(int numberOfSteps, int equilibriationSteps, const bool acceptedStep, const int stepNumber);
    void            printOutputToTerminal(const int maxIter, const double time);
    void            printFinalOutputToTerminal();
    void            openOutputFiles(const std::string path);
    void            printOutputToFile();
    void            closeOutputFiles();
    void            printInstantValuesToFile(const Eigen::VectorXd positions);
    void            computeAverages();
    double          getAverageEnergy()        { return m_averageEnergy; }
    Eigen::MatrixXd getAverageGradients()     { return m_averageGradients; }
    Eigen::MatrixXd getAverageGradientsE()    { return m_averageGradientsE; }
    std::string     generateFileName(std::string path, std::string name, std::string optimization, const std::string extension);

private:
    int              m_totalNumberOfSteps = 0;
    int              m_maxNumberOfParametersPerElement = 0;
    int              m_numberOfMetropolisSteps = 0;
    int              m_numberOfParticles = 0;
    int              m_numberOfDimensions = 0;
    int              m_numberOfElements = 0;
    int              m_numberOfBatches = 0;
    int              m_numberOfStepsPerBatch = 0;
    int              m_numberOfSteps = 0;
    int              m_equilibriationSteps = 0;
    int              m_acceptenceRatio = 0;
    int              m_iter = 0;
    bool             m_interaction = 0;
    double           m_variance = 0;
    double           m_equilibrationFraction = 0;
    double           m_omega = 0;

    double           m_averageEnergy = 0;
    double           m_averageEnergySqrd = 0;
    double           m_cumulativeEnergy = 0;
    double           m_cumulativeEnergySqrd = 0;
    double           m_instantEnergy = 0;

    Eigen::MatrixXd  m_averageGradients;
    Eigen::MatrixXd  m_averageGradientsE;
    Eigen::MatrixXd  m_cumulativeGradients;
    Eigen::MatrixXd  m_cumulativeGradientsE;
    Eigen::MatrixXd  m_instantGradients;

    std::ofstream    m_averageEnergyFile;
    std::ofstream    m_instantEnergyFile;
    std::string      m_averageEnergyFileName = "Filename not generated yet";
    std::string      m_instantEnergyFileName = "Filename not generated yet";
    class System*    m_system = nullptr;

    bool             m_printEnergyToFile = true;
    bool             m_printInstantEnergyToFile = true;

    // One-body density related stuff
    bool             m_calculateOneBody  = true;
    int              m_numberOfBins      = 100;
    double           m_maxRadius         = 10;
    double           m_radialStep        = 0.1;
    Eigen::VectorXd  m_binLinSpace;
    Eigen::VectorXd  m_particlesPerBin;
    std::ofstream    m_oneBodyFile;
};
