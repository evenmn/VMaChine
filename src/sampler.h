#pragma once
#include "Eigen/Dense"
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

class Sampler
{
public:
    Sampler(class System *system);
    void sample(const bool acceptedStep, const int stepNumber);
    void printOutputToTerminal(const int maxIter, const double time);
    void printFinalOutputToTerminal();
    void openOutputFiles();
    void printEnergyToFile();
    void printParametersToFile();
    void printOneBodyDensityToFile();
    void printTwoBodyDensityToFile();
    void closeOutputFiles();
    void printInstantValuesToFile();
    void computeOneBodyDensity(const Eigen::VectorXd radialVector);
    void computeTwoBodyDensity(const Eigen::VectorXd radialVector);
    void computeOneBodyDensity2(const Eigen::VectorXd positions);
    void computeAverages();
    void computeTotals();
    void doResampling();
    void appendInstantFiles(const std::string extension);
    void mergeOneBodyFiles();
    void setNumberOfSteps(int numberOfStepsWOEqui,
                          int totalNumberOfStepsWOEqui,
                          int totalNumberOfStepsWEqui);
    double getAverageEnergy() { return m_averageEnergy; }
    Eigen::MatrixXd getAverageGradients() { return m_averageGradients; }
    Eigen::MatrixXd getAverageGradientsE() { return m_averageGradientsE; }
    std::string getParameterFileName() { return m_parameterFileName; }
    std::string generateFileName(std::string name, std::string extension);

private:
    int m_stepsWOEqui = 0;
    int m_totalStepsWOEqui = 0;
    int m_totalStepsWEqui = 0;
    int m_equilibriationSteps = 0;
    int m_initialTotalStepsWOEqui = 0;

    int m_maxParameters = 0;
    int m_numberOfProcesses = 0;
    int m_numberOfParticles = 0;
    int m_numberOfDimensions = 0;
    int m_numberOfElements = 0;
    int m_numberOfBatches = 0;
    int m_numberOfStepsPerBatch = 0;
    int m_totalAcceptence = 0;
    int m_acceptence = 0;
    int m_iter = 0;
    int m_rank = 0;
    int m_instantNumber = 0;
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
    Eigen::MatrixXd m_totalCumulativeGradients;
    Eigen::MatrixXd m_totalCumulativeGradientsE;

    Eigen::MatrixXd m_averageGradients;
    Eigen::MatrixXd m_averageGradientsE;
    Eigen::MatrixXd m_cumulativeGradients;
    Eigen::MatrixXd m_cumulativeGradientsE;
    Eigen::MatrixXd m_instantGradients;

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
    bool m_computeTwoBodyDensity = true;
    int m_numberOfBins = 1000;
    int m_numberOfBinsHalf = m_numberOfBins / 2;
    double m_maxRadius = 10;
    double m_radialStep = 0.1;

    Eigen::VectorXi m_particlesPerBin;
    Eigen::VectorXi m_totalParticlesPerBin;
    Eigen::MatrixXi m_densityGrid;
    Eigen::MatrixXi m_totalDensityGrid;
    Eigen::MatrixXi m_particlesPerBinPairwise;
    Eigen::MatrixXi m_totalParticlesPerBinPairwise;
    std::ofstream m_oneBodyFile;
    std::ofstream m_twoBodyFile;

    class System *m_system = nullptr;
};
