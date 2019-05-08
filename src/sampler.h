#pragma once
#include <Eigen/Dense>
#include <fstream>

class Sampler {
public:
    Sampler(class System* system);
    Eigen::MatrixXd getAverageGradients         ()  { return m_averageGradients; }
    Eigen::MatrixXd getAverageGradientsE        ()  { return m_averageGradientsE; }
    double          getAverageEnergy            ()  { return m_averageEnergy; }

    void            sample                      (const bool acceptedStep, const unsigned long stepNumber);
    void            printOutputToTerminal       (const unsigned int maxIter, const double time);
    void            printFinalOutputToTerminal  ();
    void            printParametersToFile       ();
    void            printEnergyToFile           ();
    void            printOneBodyDensityToFile   ();
    void            printTwoBodyDensityToFile   ();
    void            printInstantValuesToFile    ();
    void            openOutputFiles             ();
    void            closeOutputFiles            ();
    void            appendInstantFiles          ();
    void            doResampling                ();
    void            computeAverages             ();
    void            computeTotals               ();
    void            computeOneBodyDensity       (const Eigen::VectorXd positions);
    void            computeTwoBodyDensity       (const Eigen::VectorXd positions);
    void            setNumberOfSteps            (const unsigned long numberOfStepsWOEqui, const unsigned long totalNumberOfStepsWOEqui, const unsigned long totalNumberOfStepsWEqui);

    std::string     generateFileName(const std::string name, const std::string extension);

private:
    int             m_rank                      = 0;
    int             m_numberOfProcesses         = 0;
    int             m_numberOfBins              = 100;

    unsigned short  m_numberOfElements          = 0;
    unsigned short  m_numberOfDimensions        = 0;

    unsigned int    m_numberOfParticles         = 0;
    unsigned int    m_maxNumberOfParametersPerElement = 0;
    unsigned int    m_numberOfBatches           = 0;
    unsigned int    m_numberOfStepsPerBatch     = 0;
    unsigned int    m_totalAcceptence           = 0;
    unsigned int    m_acceptence                = 0;
    unsigned int    m_iter                      = 0;
    unsigned int    m_numberOfEquilibriationSteps = 0;

    unsigned long   m_instantNumber             = 0;
    unsigned long   m_numberOfStepsWOEqui       = 0;
    unsigned long   m_totalNumberOfStepsWOEqui  = 0;
    unsigned long   m_totalNumberOfStepsWEqui   = 0;
    unsigned long   m_initialTotalNumberOfStepsWOEqui = 0;

    double          m_stdError                  = 0;
    double          m_variance                  = 0;
    double          m_mseEnergy                 = 0;
    double          m_mseSTD                    = 0;
    double          m_mseVariance               = 0;
    double          m_equilibrationFraction     = 0;
    double          m_omega                     = 0;
    double          m_averageEnergy             = 0;
    double          m_averageEnergySqrd         = 0;
    double          m_cumulativeEnergy          = 0;
    double          m_cumulativeEnergySqrd      = 0;
    double          m_instantEnergy             = 0;
    double          m_totalCumulativeEnergy     = 0;
    double          m_totalCumulativeEnergySqrd = 0;
    double          m_maxRadius                 = 10;
    double          m_radialStep                = 0.1;

    bool            m_interaction              = true;
    bool            m_printParametersToFile    = true;
    bool            m_printEnergyToFile        = true;
    bool            m_printInstantEnergyToFile = true;
    bool            m_computeOneBodyDensity    = true;
    bool            m_computeTwoBodyDensity    = true;


    Eigen::MatrixXd m_totalCumulativeGradients;
    Eigen::MatrixXd m_totalCumulativeGradientsE;
    Eigen::MatrixXd m_averageGradients;
    Eigen::MatrixXd m_averageGradientsE;
    Eigen::MatrixXd m_cumulativeGradients;
    Eigen::MatrixXd m_cumulativeGradientsE;
    Eigen::MatrixXd m_instantGradients;
    Eigen::MatrixXi m_particlesPerBinPairwise;
    Eigen::MatrixXi m_totalParticlesPerBinPairwise;
    Eigen::VectorXi m_particlesPerBin;
    Eigen::VectorXi m_totalParticlesPerBin;
    Eigen::VectorXd m_binLinSpace;

    std::ofstream   m_parameterFile;
    std::ofstream   m_averageEnergyFile;
    std::ofstream   m_instantEnergyFile;
    std::ofstream   m_oneBodyFile;
    std::ofstream   m_twoBodyFile;

    std::string     m_instantEnergyFileName = "Filename not generated yet";
    std::string     m_path                  = "Path not specified";
    std::string     m_waveFunction          = "WaveFunction configuration now given";

    class System*   m_system                = nullptr;
};
