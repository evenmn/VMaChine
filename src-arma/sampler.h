#pragma once
#include "Hamiltonians/hamiltonian.h"
#include "Optimization/optimization.h"
#include "RNG/rng.h"
#include "WaveFunctions/wavefunction.h"
#include "block/c++/blocker.h"
#include "system.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <string>
#include <armadillo>

class Sampler
{
public:
    Sampler(class System *system);
    void sample(const bool acceptedStep, const arma::uword stepNumber);
    void printOutputToTerminal(const arma::uword maxIter, const double time);
    void printFinalOutputToTerminal();
    void openOutputFiles();
    void printEnergyToFile();
    void printParametersToFile();
    void printOneBodyDensityToFile();
    void printOneBodyDensity2ToFile();
    void printTwoBodyDensityToFile();
    void closeOutputFiles();
    void printInstantValuesToFile();
    void computeOneBodyDensity(const arma::vec radialVector);
    void computeTwoBodyDensity(const arma::vec radialVector);
    void computeOneBodyDensity2(const arma::vec positions);
    void computeAverages();
    void computeTotals();
    void doResampling();
    void appendInstantFiles(const std::string extension);
    void mergeOneBodyFiles();
    void setNumberOfSteps(arma::uword numberOfStepsWOEqui,
                          arma::uword totalNumberOfStepsWOEqui,
                          arma::uword totalNumberOfStepsWEqui);
    double getAverageEnergy() { return m_averageEnergy; }
    arma::mat getAverageGradients() { return m_averageGradients; }
    arma::mat getAverageGradientsE() { return m_averageGradientsE; }
    std::string getParameterFileName() { return m_parameterFileName; }
    std::string generateFileName(std::string name, std::string extension);

private:
    arma::uword m_stepsWOEqui = 0;
    arma::uword m_totalStepsWOEqui = 0;
    arma::uword m_totalStepsWEqui = 0;
    arma::uword m_equilibriationSteps = 0;
    arma::uword m_initialTotalStepsWOEqui = 0;

    arma::uword m_maxParameters = 0;
    arma::uword m_numberOfProcesses = 0;
    arma::uword m_numberOfParticles = 0;
    arma::uword m_numberOfDimensions = 0;
    arma::uword m_numberOfElements = 0;
    arma::uword m_numberOfBatches = 0;
    arma::uword m_numberOfStepsPerBatch = 0;
    arma::uword m_totalAcceptence = 0;
    arma::uword m_acceptence = 0;
    arma::uword m_iter = 0;
    arma::uword m_rank = 0;
    arma::uword m_instantNumber = 0;
    bool m_interaction = 0;
    double m_equilibrationFraction = 0;
    double m_omega = 0;

    double m_variance = 0;
    double m_stdError = 0;
    double m_mseEnergy = 0;
    double m_mseSTD = 0;
    double m_mseVariance = 0;
    double m_varianceKinetic = 0;
    double m_stdErrorKinetic = 0;
    double m_mseEnergyKinetic = 0;
    double m_mseSTDKinetic = 0;
    double m_mseVarianceKinetic = 0;
    double m_varianceExternal = 0;
    double m_stdErrorExternal = 0;
    double m_mseEnergyExternal = 0;
    double m_mseSTDExternal = 0;
    double m_mseVarianceExternal = 0;
    double m_varianceInteraction = 0;
    double m_stdErrorInteraction = 0;
    double m_mseEnergyInteraction = 0;
    double m_mseSTDInteraction = 0;
    double m_mseVarianceInteraction = 0;
    double m_averageKineticEnergy = 0;
    double m_averageExternalEnergy = 0;
    double m_averageInteractionEnergy = 0;
    double m_averageEnergy = 0;
    double m_averageEnergySqrd = 0;
    double m_cumulativeKineticEnergy = 0;
    double m_cumulativeExternalEnergy = 0;
    double m_cumulativeInteractionEnergy = 0;
    double m_cumulativeEnergy = 0;
    double m_cumulativeEnergySqrd = 0;
    double m_kineticEnergy = 0;
    double m_externalEnergy = 0;
    double m_interactionEnergy = 0;
    double m_instantEnergy = 0;

    double m_totalCumulativeKineticEnergy = 0;
    double m_totalCumulativeExternalEnergy = 0;
    double m_totalCumulativeInteractionEnergy = 0;
    double m_totalCumulativeEnergy = 0;
    double m_totalCumulativeEnergySqrd = 0;
    arma::mat m_totalCumulativeGradients;
    arma::mat m_totalCumulativeGradientsE;

    arma::mat m_averageGradients;
    arma::mat m_averageGradientsE;
    arma::mat m_cumulativeGradients;
    arma::mat m_cumulativeGradientsE;
    arma::mat m_instantGradients;

    std::ofstream m_averageEnergyFile;
    std::ofstream m_averageKineticEnergyFile;
    std::ofstream m_averageExternalEnergyFile;
    std::ofstream m_averageInteractionEnergyFile;
    std::ofstream m_instantEnergyFile;
    std::ofstream m_instantKineticEnergyFile;
    std::ofstream m_instantExternalEnergyFile;
    std::ofstream m_instantInteractionEnergyFile;
    std::ofstream m_parameterFile;
    std::string m_parameterFileName = "Filename not generated yet";
    std::string m_instantEnergyFileName = "Filename not generated yet";
    std::string m_instantKineticEnergyFileName = "Filename not generated yet";
    std::string m_instantExternalEnergyFileName = "Filename not generated yet";
    std::string m_instantInteractionEnergyFileName = "Filename not generated yet";
    std::string m_path = "Path not specified";
    std::string m_trialWaveFunction = "Wave function not found";
    std::string m_hamiltonian = "Hamiltonian not found";

    bool m_printEnergyToFile = true;
    bool m_printInstantEnergyToFile = true;
    bool m_printParametersToFile = true;

    // Electron density related stuff
    bool m_computeOneBodyDensity = true;
    bool m_computeOneBodyDensity2 = true;
    bool m_computeTwoBodyDensity = true;
    arma::uword m_numberOfBins = 1000;
    arma::uword m_numberOfBinsHalf = m_numberOfBins / 2;
    double m_maxRadius = 10.;
    double m_radialStep = 0.1;

    arma::vec m_particlesPerBin;
    arma::vec m_totalParticlesPerBin;
    arma::mat m_densityGrid;
    arma::mat m_totalDensityGrid;
    arma::mat m_particlesPerBinPairwise;
    arma::mat m_totalParticlesPerBinPairwise;
    std::ofstream m_oneBodyFile;
    std::ofstream m_oneBodyFile2;
    std::ofstream m_twoBodyFile;

    class System *m_system = nullptr;
};
