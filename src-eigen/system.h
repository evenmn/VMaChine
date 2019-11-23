#pragma once
#include <cassert>
#include <ctime>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

#include "Eigen/Dense"
#include "main.h"

class System
{
public:
    int adaptiveSteps();
    void initializeSystem();
    void runMetropolisCycles();
    void checkingConvergence();
    void setNumberOfFreeDimensions();
    void setMaxParameters();
    void setGradients();
    void setGlobalArraysToCalculate();
    void resetAllArrays();
    void setAllConstants();
    void collectAllLabels();
    void setAtomicNumber(const int Z);
    void runSimulation();
    void printToTerminal();
    void setNumberOfParticles(const int numberOfParticles);
    void setNumberOfDimensions(const int numberOfDimensions);
    void setNumberOfHiddenUnits(const int numberOfHiddenUnits);
    void setNumberOfMetropolisSteps(const int steps);
    void setNumberOfIterations(const int numberOfIterations);
    void setNumberOfElements(const unsigned long numberOfElements);
    void setStepLength(const double stepLength);
    void setEquilibrationFraction(const double equilibrationFraction);
    void setFrequency(const double omega);
    void setWidth(const double sigma);
    void setTotalSpin(const double totalSpin);
    void setLearningRate(const double eta);
    void setPath(const std::string path);
    void initializeFromConfig(int argc, char** argv);
    void parser(const std::string configFile);
    void searchShortning(const std::vector<std::string> labels,
                         const std::string newLabel,
                         std::string &allLabels);

    void setInteraction(const bool interaction);
    void setParameterPrintingTools(const bool printParametersToFile = true);
    void setConvergenceTools(const int numberOfEnergies = 5,
                             const double tolerance = 1e-6);
    void setAdaptiveStepTools(const int rangeOfAdaptiveSteps = 10,
                              const int additionalSteps = 4,
                              const int additionalStepsLastIteration = 8);
    void setDensityTools(const bool computeOneBodyDensity = false,
                         const bool computeSpatialOneBodyDensity = false,
                         const bool computeTwoBodyDensity = false,
                         const int numberOfBins = 1000,
                         const double maxRadius = 50);
    void setScreeningTools(const double screeningStrength, const double dsl);
    void dumpEnergyToFile(bool printEnergyFile = true);
    void checkResampling(bool printInstantEnergyFile = true);
    void initializeMPI();

    void updateAllParameters(const Eigen::MatrixXd parameters);
    void initializeAllArrays(const Eigen::VectorXd positions,
                             const Eigen::VectorXd radialVector,
                             const Eigen::MatrixXd distanceMatrix);
    void updateAllArrays(const Eigen::VectorXd positions,
                         const Eigen::VectorXd radialVector,
                         const Eigen::MatrixXd distanceMatrix,
                         const int changedCoord);

    double evaluateProbabilityRatio();
    double getKineticEnergy();
    Eigen::MatrixXd getAllParameterGradients();

    void setHamiltonian(class Hamiltonian *hamiltonian);
    void setBasis(class Basis *basis);
    void setInitialState(class InitialState *initialState);
    void setInitialWeights(class InitialWeights *initialWeights);
    void setMetropolis(class Metropolis *metropolis);
    void setOptimization(class Optimization *optimization);
    void setRandomNumberGenerator(class RandomNumberGenerator *randomNumberGenerator);
    void setWaveFunctionElements(std::vector<class WaveFunction *> waveFunctionElements);
    void setWaveFunctionElement(WaveFunction *waveFunction);
    void setInputLayer(int numberOfUnits);
    void addDenseLayer(int numberOfUnits, Activation *activation);

    class WaveFunction *getWaveFunction() { return m_waveFunction; }
    class Hamiltonian *getHamiltonian() { return m_hamiltonian; }
    class Basis *getBasis() { return m_basis; }
    class Sampler *getSampler() { return m_sampler; }
    class Optimization *getOptimization() { return m_optimization; }
    class InitialWeights *getInitialWeights() { return m_initialWeights; }
    class InitialState *getInitialState() { return m_initialState; }
    class RandomNumberGenerator *getRandomNumberGenerator() { return m_randomNumberGenerator; }

    int getNumberOfElements() { return m_numberOfElements; }
    int getNumberOfProcesses() { return m_numberOfProcesses; }
    int getNumberOfParticles() { return m_numberOfParticles; }
    int getNumberOfDimensions() { return m_numberOfDimensions; }
    int getNumberOfHiddenUnits() { return m_numberOfHiddenUnits; }
    int getNumberOfFreeDimensions() { return m_degreesOfFreedom; }
    int getTotalNumberOfParameters() { return m_totalNumberOfParameters; }
    int getMaxParameters() { return m_maxParameters; }
    int getAtomicNumber() { return m_Z; }
    int getRank() { return m_myRank; }
    int getNumberOfBins() { return m_numberOfBins; }

    int getTotalStepsWOEqui() { return m_totalStepsWOEqui; }
    int getTotalStepsWEqui() { return m_totalStepsWEqui; }
    int getTotalEquilibriationSteps() { return m_totalEquilibriationSteps; }
    int getStepsWOEqui() { return m_stepsWOEqui; }
    int getStepsWEqui() { return m_stepsWEqui; }
    int getInitialTotalStepsWOEqui() { return m_initialTotalStepsWOEqui; }
    int getEquilibriationSteps() { return m_equilibriationSteps; }

    double getMaxRadius() { return m_maxRadius; }
    double getEquilibrationFraction() { return m_equilibrationFraction; }
    double getFrequency() { return m_omega; }
    double getWidth() { return m_sigma; }
    double getLearningRate() { return m_eta; }
    double getStepLength() { return m_stepLength; }
    double getTotalSpin() { return m_totalSpin; }
    double getScreeningStrength() { return m_screeningStrength; }
    double getDSL() { return m_dsl; }

    bool getScreening() { return m_screening; }
    bool getInteraction() { return m_interaction; }
    bool computeOneBodyDensity() { return m_computeOneBodyDensity; }
    bool computeOneBodyDensity2() { return m_computeOneBodyDensity2; }
    bool computeTwoBodyDensity() { return m_computeTwoBodyDensity; }
    bool printEnergyToFile() { return m_printEnergyToFile; }
    bool printParametersToFile() { return m_printParametersToFile; }
    bool doResampling() { return m_doResampling; }
    bool getCalculateDistanceMatrix() { return m_calculateDistanceMatrix; }
    bool getCalculateRadialVector() { return m_calculateRadialVector; }

    Eigen::VectorXd getPositions() { return m_positions; }
    Eigen::VectorXd getRadialVector() { return m_radialVector; }
    Eigen::MatrixXd getDistanceMatrix() { return m_distanceMatrix; }
    Eigen::MatrixXd getWeights() { return m_parameters; }

    std::string getPath() { return m_path; }
    std::string getTrialWaveFunction() { return m_trialWaveFunction; }

    std::vector<class WaveFunction *> getWaveFunctionElements() { return m_waveFunctionElements; }
    std::vector<class Layer *> getLayers() { return m_layers; }
    std::vector<int> getHiddenUnits() { return m_hiddenUnits; }

private:
    int m_numberOfParticles;
    int m_numberOfDimensions;
    int m_degreesOfFreedom;
    int m_numberOfIterations;
    int m_numberOfElements;
    int m_numberOfHiddenUnits;
    int m_maxParameters;
    int m_totalNumberOfParameters;
    int m_Z;
    int m_numberOfProcesses;
    int m_myRank;

    int m_rangeOfAdaptiveSteps = 10;
    int m_additionalSteps = 4;
    int m_additionalStepsLastIter = 8;
    int m_lastIteration = 1;
    int m_numberOfEnergies = 5;
    int m_numberOfBins = 1000;
    int m_iter = 0;
    int m_args = 1;

    int m_totalStepsWOEqui = int(pow(2,19));
    int m_totalStepsWEqui = int(pow(2,19));
    int m_totalEquilibriationSteps = 0;
    int m_stepsWOEqui = int(pow(2, 17));
    int m_stepsWEqui = int(pow(2, 17));
    int m_equilibriationSteps = 0;
    int m_initialStepsWOEqui = int(pow(2, 17));
    int m_initialTotalStepsWOEqui = int(pow(2, 19));

    double m_equilibrationFraction = 0;
    double m_tolerance = 1e-7;
    double m_maxRadius = 5;
    double m_totalTime = 0;
    double m_totalSpin = 0;
    double m_screeningStrength = 100;
    double m_dsl = 100;
    double m_globalTime = 0;

    double m_omega = 1.0;
    double m_sigma = 1.0;
    double m_stepLength = 0.05;
    double m_eta = 0.1;

    bool m_interaction = true;
    bool m_checkConvergence = false;
    bool m_applyAdaptiveSteps = false;
    bool m_computeOneBodyDensity = false;
    bool m_computeOneBodyDensity2 = false;
    bool m_computeTwoBodyDensity = false;
    bool m_printEnergyToFile = false;
    bool m_doResampling = true;
    bool m_printParametersToFile = false;
    bool m_calculateDistanceMatrix = false;
    bool m_calculateRadialVector = false;
    bool m_screening = false;

    class WaveFunction *m_waveFunction = nullptr;
    class Hamiltonian *m_hamiltonian = new HarmonicOscillator(this);
    class Basis *m_basis = new Hermite(this);
    class InitialState *m_initialState = new RandomNormal(this);
    class InitialWeights *m_initialWeights = new Automatize(this);
    class Sampler *m_sampler = nullptr;
    class Metropolis *m_metropolis = new ImportanceSampling(this);
    class Optimization *m_optimization = new ADAM(this);
    class RandomNumberGenerator *m_randomNumberGenerator = new MersenneTwister();

    std::vector<class WaveFunction *> m_waveFunctionElements;
    std::vector<class Layer *> m_layers;
    std::vector<int> m_hiddenUnits;

    std::string m_path = "../data/";
    std::string m_trialWaveFunction;
    std::string m_configFile = "none";

    Eigen::VectorXd m_positions;
    Eigen::MatrixXd m_parameters;
    Eigen::MatrixXd m_distanceMatrix;
    Eigen::VectorXd m_radialVector;
    Eigen::MatrixXd m_gradients;
    Eigen::VectorXd m_energies;
};
