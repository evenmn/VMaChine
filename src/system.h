#pragma once
#include <Eigen/Dense>
#include <vector>
#include <string>

class System {
public:
    int  dynamicSteps               ();
    void runMetropolisCycles        ();
    void checkingConvergence        ();
    void setNumberOfFreeDimensions  ();
    void setMaxParameters           ();
    void setGradients               ();
    void setGlobalArraysToCalculate ();
    void resetAllArrays             ();
    void setAllConstants            ();
    void collectAllLabels           ();
    void setAtomicNumber            (const int Z);
    void runIterations              (const int numberOfIterations);
    void printToTerminal            (const int numberOfIterations);
    void setNumberOfParticles       (const int numberOfParticles);
    void setNumberOfDimensions      (const int numberOfDimensions);
    void setNumberOfHiddenNodes     (const int numberOfHiddenNodes);
    void setNumberOfMetropolisSteps (const int steps);
    void setNumberOfElements        (const unsigned long numberOfElements);
    void setStepLength              (const double stepLength);
    void setEquilibrationFraction   (const double equilibrationFraction);
    void setFrequency               (const double omega);
    void setWidth                   (const double sigma);
    void setLearningRate            (const double eta);
    void setPath                    (const std::string path);
    void parser                     (const std::string configFile, int &numberOfIterations);

    void setInteraction             (const bool interaction);
    void setParameterPrintingTools  (const bool printParametersToFile);
    void setConvergenceTools        (const bool checkConvergence,
                                     const int numberOfEnergies,
                                     const double tolerance);
    void setDynamicStepTools        (const bool applyDynamicSteps,
                                     const int rangeOfDynamicSteps,
                                     const int additionalSteps,
                                     const int additionalStepsLastIteration);
    void setDensityTools            (const bool computeDensity,
                                     const bool computeTwoBodyDensity,
                                     const int numberOfBins,
                                     const double maxRadius);
    void setEnergyPrintingTools     (const bool printEnergyFile,
                                     const bool printInstantEnergyFile);
    void setMPITools                (const int myRank, int numberOfProcesses);

    void updateAllParameters        (const Eigen::MatrixXd parameters);
    void initializeAllArrays        (const Eigen::VectorXd positions,
                                     const Eigen::VectorXd radialVector,
                                     const Eigen::MatrixXd distanceMatrix);
    void updateAllArrays            (const Eigen::VectorXd positions,
                                     const Eigen::VectorXd radialVector,
                                     const Eigen::MatrixXd distanceMatrix,
                                     const int changedCoord);

    double                          evaluateWaveFunctionRatio    ();
    double                          getKineticEnergy             ();
    Eigen::MatrixXd                 getAllInstantGradients       ();

    void setHamiltonian             (class Hamiltonian* hamiltonian);
    void setBasis                   (class Basis* basis);
    void setInitialState            (class InitialState* initialState);
    void setInitialWeights          (class InitialWeights* initialWeights);
    void setMetropolis              (class Metropolis* metropolis);
    void setOptimization            (class Optimization* optimization);
    void setRandomNumberGenerator   (class RandomNumberGenerator* randomNumberGenerator);
    void setWaveFunctionElements    (std::vector<class WaveFunction*> waveFunctionElements);

    class WaveFunction*             getWaveFunction()            { return m_waveFunction; }
    class Hamiltonian*              getHamiltonian()             { return m_hamiltonian; }
    class Basis*                    getBasis()                   { return m_basis; }
    class Sampler*                  getSampler()                 { return m_sampler; }
    class Optimization*             getOptimization()            { return m_optimization; }
    class InitialWeights*           getInitialWeights()          { return m_initialWeights; }
    class InitialState*             getInitialState()            { return m_initialState; }
    class RandomNumberGenerator*    getRandomNumberGenerator()   { return m_randomNumberGenerator; }

    unsigned long                   getNumberOfElements()        { return m_numberOfElements; }
    int                             getNumberOfProcesses()       { return m_numberOfProcesses; }
    int                             getNumberOfParticles()       { return m_numberOfParticles; }
    int                             getNumberOfDimensions()      { return m_numberOfDimensions; }
    int                             getNumberOfHiddenNodes()     { return m_numberOfHiddenNodes; }
    int                             getNumberOfFreeDimensions()  { return m_numberOfFreeDimensions; }
    int                             getTotalNumberOfParameters() { return m_totalNumberOfParameters; }
    int                             getMaxParameters()           { return m_maxParameters; }
    int                             getAtomicNumber()            { return m_Z; }
    int                             getRank()                    { return m_myRank; }
    int                             getNumberOfBins()            { return m_numberOfBins; }

    int                             getTotalStepsWOEqui()        { return m_totalStepsWOEqui; }
    int                             getTotalStepsWEqui()         { return m_totalStepsWEqui; }
    int                             getTotalEquilibriationSteps(){ return m_totalEquilibriationSteps; }
    int                             getStepsWOEqui()             { return m_stepsWOEqui; }
    int                             getStepsWEqui()              { return m_stepsWEqui; }
    int                             getInitialTotalStepsWOEqui() { return m_initialTotalStepsWOEqui; }
    int                             getEquilibriationSteps()     { return m_equilibriationSteps; }

    double                          getMaxRadius()               { return m_maxRadius; }
    double                          getEquilibrationFraction()   { return m_equilibrationFraction; }
    double                          getFrequency()               { return m_omega; }
    double                          getWidth()                   { return m_sigma; }
    double                          getLearningRate()            { return m_eta; }
    double                          getStepLength()              { return m_stepLength; }
    bool                            getInteraction()             { return m_interaction; }
    bool                            getDensity()                 { return m_computeDensity; }
    bool                            computeTwoBodyDensity()      { return m_computeTwoBodyDensity; }
    bool                            getPrintEnergy()             { return m_printEnergyFile; }
    bool                            getPrintParametersToFile()   { return m_printParametersToFile; }
    bool                            getPrintInstantEnergy()      { return m_printInstantEnergyFile; }
    bool                            getCalculateDistanceMatrix() { return m_calculateDistanceMatrix; }
    bool                            getCalculateRadialVector()   { return m_calculateRadialVector; }
    Eigen::VectorXd                 getPositions()               { return m_positions; }
    Eigen::VectorXd                 getRadialVector()            { return m_radialVector; }
    Eigen::MatrixXd                 getDistanceMatrix()          { return m_distanceMatrix; }
    Eigen::MatrixXd                 getWeights()                 { return m_parameters; }
    std::string                     getPath()                    { return m_path; }
    std::string                     getTrialWaveFunction()       { return m_trialWaveFunction; }
    std::vector<class WaveFunction*>getWaveFunctionElements()    { return m_waveFunctionElements; }


private:
    unsigned long                       m_numberOfElements          = 0;
    int                                 m_numberOfHiddenNodes       = 0;
    int                                 m_numberOfParticles         = 0;
    int                                 m_numberOfDimensions        = 0;
    int                                 m_numberOfFreeDimensions    = 0;
    int                                 m_maxParameters             = 0;
    int                                 m_totalNumberOfParameters   = 0;
    int                                 m_Z                         = 1;
    int                                 m_rangeOfDynamicSteps       = 10;
    int                                 m_additionalSteps           = 4;
    int                                 m_additionalStepsLastIteration = 8;
    int                                 m_lastIteration             = 1;
    int                                 m_numberOfEnergies          = 0;
    int                                 m_numberOfBins              = 1;
    int                                 m_numberOfProcesses         = 1;
    int                                 m_myRank                    = 0;
    int                                 m_iter                      = 0;

    int                                 m_totalStepsWOEqui          = 0;
    int                                 m_totalStepsWEqui           = 0;
    int                                 m_totalEquilibriationSteps  = 0;
    int                                 m_stepsWOEqui               = 0;
    int                                 m_stepsWEqui                = 0;
    int                                 m_equilibriationSteps       = 0;
    int                                 m_initialStepsWOEqui        = 0;
    int                                 m_initialTotalStepsWOEqui   = 0;

    double                              m_equilibrationFraction     = 0.0;
    double                              m_stepLength                = 0.1;
    double                              m_omega                     = 1.0;
    double                              m_sigma                     = 1.0;
    double                              m_eta                       = 0.1;
    double                              m_tolerance                 = 1e-7;
    double                              m_maxRadius                 = 10;
    double                              m_totalTime                 = 0;

    bool                                m_interaction               = true;
    bool                                m_checkConvergence          = true;
    bool                                m_applyDynamicSteps         = true;
    bool                                m_computeDensity            = true;
    bool                                m_computeTwoBodyDensity     = true;
    bool                                m_printEnergyFile           = true;
    bool                                m_printInstantEnergyFile    = true;
    bool                                m_printParametersToFile     = true;
    bool                                m_calculateDistanceMatrix   = true;
    bool                                m_calculateRadialVector     = true;

    class WaveFunction*                 m_waveFunction              = nullptr;
    class Hamiltonian*                  m_hamiltonian               = nullptr;
    class Basis*                        m_basis                     = nullptr;
    class InitialState*                 m_initialState              = nullptr;
    class InitialWeights*               m_initialWeights            = nullptr;
    class Sampler*                      m_sampler                   = nullptr;
    class Metropolis*                   m_metropolis                = nullptr;
    class Optimization*                 m_optimization              = nullptr;
    class RandomNumberGenerator*        m_randomNumberGenerator     = nullptr;
    std::vector<class WaveFunction*>    m_waveFunctionElements;

    std::string                         m_path                      = "../data/";
    std::string                         m_trialWaveFunction         = "Not defined";

    Eigen::VectorXd                     m_positions;
    Eigen::MatrixXd                     m_parameters;
    Eigen::MatrixXd                     m_distanceMatrix;
    Eigen::VectorXd                     m_radialVector;
    Eigen::MatrixXd                     m_gradients;
    Eigen::VectorXd                     m_energies;
};

