#pragma once
#include <Eigen/Dense>
#include <fstream>

class Sampler {
public:
    Sampler(class System* system);
    void            sample(const bool acceptedStep, const int stepNumber);
    void            printOutputToTerminal(const int maxIter, const double time);
    void            printFinalOutputToTerminal();
    void            openOutputFiles();
    void            printEnergyToFile();
    void            printParametersToFile();
    void            printOneBodyDensityToFile();
    void            printTwoBodyDensityToFile();
    void            closeOutputFiles();
    void            printInstantValuesToFile();
    void            computeOneBodyDensity(const Eigen::VectorXd positions);
    void            computeTwoBodyDensity(const Eigen::VectorXd positions);
    void            computeAverages();
    void            computeTotals();
    void            doResampling();
    void            appendInstantFiles();
    void            mergeOneBodyFiles();
    void            setNumberOfSteps(int numberOfStepsWOEqui, int totalNumberOfStepsWOEqui, int totalNumberOfStepsWEqui);
    double          getAverageEnergy()        { return m_averageEnergy; }
    Eigen::MatrixXd getAverageGradients()     { return m_averageGradients; }
    Eigen::MatrixXd getAverageGradientsE()    { return m_averageGradientsE; }
    std::string     generateFileName(std::string name, std::string extension);

private:
    int              m_numberOfStepsWOEqui = 0;
    int              m_totalNumberOfStepsWOEqui = 0;
    int              m_totalNumberOfStepsWEqui = 0;
    int              m_numberOfEquilibriationSteps = 0;
    int              m_initialTotalNumberOfStepsWOEqui = 0;

    int              m_maxNumberOfParametersPerElement = 0;
    int              m_numberOfProcesses = 0;
    int              m_numberOfParticles = 0;
    int              m_numberOfDimensions = 0;
    int              m_numberOfElements = 0;
    int              m_numberOfBatches = 0;
    int              m_numberOfStepsPerBatch = 0;
    int              m_totalAcceptence = 0;
    int              m_acceptence = 0;
    int              m_iter = 0;
    int              m_rank = 0;
    int              m_instantNumber = 0;
    bool             m_interaction = 0;
    double           m_variance = 0;
    double           m_equilibrationFraction = 0;
    double           m_omega = 0;

    double           m_stdError                  = 0;
    double           m_mseEnergy                 = 0;
    double           m_mseSTD                    = 0;
    double           m_mseVariance               = 0;
    double           m_averageEnergy = 0;
    double           m_averageEnergySqrd = 0;
    double           m_cumulativeEnergy = 0;
    double           m_cumulativeEnergySqrd = 0;
    double           m_instantEnergy = 0;

    double           m_totalCumulativeEnergy = 0;
    double           m_totalCumulativeEnergySqrd = 0;
    Eigen::MatrixXd  m_totalCumulativeGradients;
    Eigen::MatrixXd  m_totalCumulativeGradientsE;

    Eigen::MatrixXd  m_averageGradients;
    Eigen::MatrixXd  m_averageGradientsE;
    Eigen::MatrixXd  m_cumulativeGradients;
    Eigen::MatrixXd  m_cumulativeGradientsE;
    Eigen::MatrixXd  m_instantGradients;

    std::ofstream    m_averageEnergyFile;
    std::ofstream    m_instantEnergyFile;
    std::ofstream    m_parameterFile;
    std::string      m_averageEnergyFileName = "Filename not generated yet";
    std::string      m_instantEnergyFileName = "Filename not generated yet";
    std::string      m_path                  = "Path not specified";
    std::string      m_waveFunction          = "Wave function not found";

    bool             m_printEnergyToFile = true;
    bool             m_printInstantEnergyToFile = true;
    bool             m_printParametersToFile = true;

    // Electron density related stuff
    bool             m_computeOneBodyDensity = true;
    bool             m_computeTwoBodyDensity = true;
    int              m_numberOfBins      = 100;
    double           m_maxRadius         = 10;
    double           m_radialStep        = 0.1;
    Eigen::VectorXd  m_binLinSpace;
    Eigen::VectorXi  m_particlesPerBin;
    Eigen::VectorXi  m_totalParticlesPerBin;
    Eigen::MatrixXi  m_particlesPerBinPairwise;
    Eigen::MatrixXi  m_totalParticlesPerBinPairwise;
    std::ofstream    m_oneBodyFile;
    std::ofstream    m_twoBodyFile;

    class System*    m_system = nullptr;
};
