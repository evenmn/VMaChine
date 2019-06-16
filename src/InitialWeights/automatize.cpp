#include "automatize.h"
#include "constant.h"
#include "randomize.h"
#include <iostream>
#include <fstream>
#include "../Optimization/optimization.h"
#include "../system.h"

Automatize::Automatize(System* system) :
    InitialWeights(system) {
    m_system                    = system;
    m_path                      = m_system->getPath();
    m_numberOfParticles         = m_system->getNumberOfParticles();
    m_numberOfDimensions        = m_system->getNumberOfDimensions();
    m_numberOfElements          = m_system->getNumberOfElements();
    m_maxParameters             = m_system->getMaxParameters();
    m_interaction               = m_system->getInteraction();
    m_omega                     = m_system->getFrequency();
    m_initialTotalStepsWOEqui   = m_system->getInitialTotalStepsWOEqui();
    m_trialWaveFunction         = m_system->getTrialWaveFunction();
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
    fileName += "_MC" + std::to_string(m_initialTotalStepsWOEqui);
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
        m_parameters = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxParameters);
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
    }
    else {
        if(m_trialWaveFunction == "VMC") {
            //Constant(m_system, 1.0);
            m_method = new Constant(m_system, 1.0);
        }
        else if(m_trialWaveFunction == "RBM") {
            m_method = new Randomize(m_system, 0.5);
        }
        else if(m_trialWaveFunction == "RBMSJ") {
            m_method = new Randomize(m_system, 0.1);
        }
        else if(m_trialWaveFunction == "RBMPJ") {
            m_method = new Randomize(m_system, 0.1);
        }
        else if(m_trialWaveFunction == "PRBM") {
            m_method = new Constant(m_system, 0.0);
        }
        else {
            m_method = new Randomize(m_system, 0.01);
        }
        m_parameters = m_method->getParameters();
    }
    m_system->updateAllParameters(m_parameters);
}

Eigen::MatrixXd Automatize::getParameters() {
    return m_parameters;
}
