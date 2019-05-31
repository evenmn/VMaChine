#include "automatize.h"
#include "constant.h"
#include "randomize.h"
#include <iostream>
#include <fstream>
#include "../Optimization/optimization.h"
#include "../system.h"

Automatize::Automatize(System* system)  :  InitialWeights(system) {
    m_system             = system;
    m_path               = m_system->getPath();
    m_numberOfParticles  = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_numberOfElements                = m_system->getNumberOfElements();
    m_maxNumberOfParametersPerElement = m_system->getMaxParameters();
    m_interaction        = m_system->getInteraction();
    m_omega              = m_system->getFrequency();
    m_initialTotalStepsWOEqui = m_system->getInitialTotalStepsWOEqui();
    m_trialWaveFunction  = m_system->getTrialWaveFunction();
    setupInitialWeights();
}

std::string Automatize::generateFileName(std::string name, std::string extension) {
    std::string fileName = m_path;
    fileName += "int" + std::to_string(m_interaction) + "/";
    fileName += name                                  + "/";
    fileName += m_trialWaveFunction                   + "/";
    fileName += std::to_string(m_numberOfDimensions)  + "D/";
    fileName += std::to_string(m_numberOfParticles)   + "P/";
    fileName += std::to_string(m_omega)               + "w/";
    fileName += m_system->getOptimization()->getLabel();
    fileName += "_MC" + std::to_string(1048575);
    fileName += extension;
    return fileName;
}

void writeFileContentIntoEigenMatrix(std::ifstream infile, Eigen::MatrixXd &matrix) {
    std::string line;
    double value;
    char delimiter;
    int i = 0;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        int j = 0;
        while (iss >> value) {
            matrix(i,j) = value;
            iss >> delimiter;
            j++;
        }
        i++;
    }
}

void Automatize::setupInitialWeights() {
    std::ifstream infile(generateFileName("weights", ".dat"));
    if(infile.is_open()) {
        m_parameters = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxNumberOfParametersPerElement);
        std::string line;
        double value;
        char delimiter;
        int i = 0;
        while(std::getline(infile, line)) {
            std::istringstream iss(line);
            int j = 0;
            while (iss >> value) {
                m_parameters(i,j) = value;
                iss >> delimiter;
                j++;
            }
            i++;
        }
        m_system->updateAllParameters(m_parameters);
    }
    else {
        if(m_trialWaveFunction == "VMC") {
            Constant(m_system, 1.0);
        }
        else if(m_trialWaveFunction == "RBM") {
            Randomize(m_system, 0.5);
            //m_parameters = Randomize::getParameters();
        }
        else if(m_trialWaveFunction == "RBMPJ") {
            Randomize(m_system, 0.1);
        }
        else if(m_trialWaveFunction == "PRBM") {
            Constant(m_system, 0.0);
        }
        else if(m_trialWaveFunction == "DRBM") {
            Constant(m_system, 0.0);
        }
        else {
            Randomize(m_system, 0.01);
        }
    }
}
