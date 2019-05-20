#include "automatize.h"
#include "constant.h"
#include "randomize.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include "../system.h"

Automatize::Automatize(System* system)  :  InitialWeights(system) {
    m_system             = system;
    m_trialWaveFunction  = m_system->getTrialWaveFunction();
    setupInitialWeights();
}

void Automatize::setupInitialWeights() {
    std::vector<std::vector<int>> v;       // declare vector of vectors
    std::ifstream infile("myfile.txt");       // open the file
    if(infile.is_open()) {
        std::string tempstr;                   // declare a temporary string
        int tempint;                           // declare a temporary integer
        char delimiter;                        // declare a temporary delimiter
        while (std::getline(infile, tempstr)) {   // read line by line from a file into a string
            std::istringstream iss(tempstr);   // initialize the stringstream with that string
            std::vector<int> tempv;            // declare a temporary vector for the row
            while (iss >> tempint) {           // extract the numbers from a stringstream
                tempv.push_back(tempint);      // push it onto our temporary vector
                iss >> delimiter;              // read the , delimiter
            }
            v.push_back(tempv);                // push the vector onto vector of vectors
        }
    }
    else {
        if(m_trialWaveFunction == "VMC") {
            Constant(m_system, 1.0);
        }
        else if(m_trialWaveFunction == "RBM") {
            Randomize(m_system, 0.5);
        }
        else if(m_trialWaveFunction == "RBMPJ") {
            Randomize(m_system, 0.1);
        }
        else {
            Randomize(m_system, 0.01);
        }
    }
}
