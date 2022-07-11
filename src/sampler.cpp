#include "sampler.h"

using std::cout;
using std::endl;

/* ----------------------------------------------------------------------------
  Sampler constructor. Initialize everything that is needed in class.
---------------------------------------------------------------------------- */

Sampler::Sampler(System *system)
{
    m_system = system;
    m_numberOfProcesses = m_system->getNumberOfProcesses();
    m_numberOfParticles = m_system->getNumberOfParticles();
    m_numberOfDimensions = m_system->getNumberOfDimensions();
    m_numberOfElements = m_system->getNumberOfElements();
    m_maxParameters = m_system->getMaxParameters();

    m_totalStepsWOEqui = m_system->getTotalStepsWOEqui();
    m_totalStepsWEqui = m_system->getTotalStepsWEqui();
    m_stepsWOEqui = m_system->getStepsWOEqui();
    m_initialTotalStepsWOEqui = m_system->getInitialTotalStepsWOEqui();
    m_equilibriationSteps = m_system->getEquilibriationSteps();

    m_omega = m_system->getFrequency();
    m_numberOfBatches = m_system->getOptimization()->getNumberOfBatches();
    m_numberOfStepsPerBatch = int(m_stepsWOEqui / m_numberOfBatches);
    m_computeOneBodyDensity = m_system->radialOneBodyDensity();
    m_computeOneBodyDensity2 = m_system->spatialOneBodyDensity();
    m_computeTwoBodyDensity = m_system->radialTwoBodyDensity();
    m_printEnergyToFile = m_system->printEnergyToFile();
    m_printInstantEnergyToFile = m_system->doResampling();
    m_printParametersToFile = m_system->printParametersToFile();
    m_numberOfBins = m_system->getNumberOfBins();
    m_numberOfBinsHalf = m_numberOfBins / 2;
    m_maxRadius = m_system->getMaxRadius();
    m_rank = m_system->getRank();
    m_path = m_system->getPath();
    m_trialWaveFunction = m_system->getTrialWaveFunction();
    m_hamiltonian = m_system->getHamiltonian()->getLabel();
    m_radialStep = m_maxRadius / m_numberOfBins;
    m_densityGrid = Eigen::MatrixXi::Zero(m_numberOfBins, m_numberOfBins);
    m_particlesPerBin = Eigen::VectorXi::Zero(m_numberOfBins);
    m_particlesPerBinPairwise = Eigen::MatrixXi::Zero(m_numberOfBins, m_numberOfBins);
}


/* ----------------------------------------------------------------------------
  Sample particles
---------------------------------------------------------------------------- */

void Sampler::sample(const bool acceptedStep, const int stepNumber)
{
    if (stepNumber == m_equilibriationSteps) {
        m_acceptence = 0;
        m_cumulativeKineticEnergy = 0;
        m_cumulativeExternalEnergy = 0;
        m_cumulativeInteractionEnergy = 0;
        m_cumulativeEnergy = 0;
        m_cumulativeEnergySqrd = 0;
        m_cumulativeGradients = Eigen::MatrixXd::Zero(Eigen::Index(m_numberOfElements),
                                                      m_maxParameters);
        m_cumulativeGradientsE = Eigen::MatrixXd::Zero(Eigen::Index(m_numberOfElements),
                                                       m_maxParameters);
    }
    m_kineticEnergy = m_system->getKineticEnergy();
    m_externalEnergy = m_system->getHamiltonian()->getExternalEnergy();
    m_interactionEnergy = m_system->getInteractionStyle()->getInteractionEnergy();
    m_instantEnergy = m_kineticEnergy + m_externalEnergy + m_interactionEnergy;
    m_instantGradients = m_system->getAllParameterGradients();

    m_cumulativeEnergy += m_instantEnergy;
    m_cumulativeKineticEnergy += m_kineticEnergy;
    m_cumulativeExternalEnergy += m_externalEnergy;
    m_cumulativeInteractionEnergy += m_interactionEnergy;

    m_cumulativeEnergySqrd += m_instantEnergy * m_instantEnergy;
    m_cumulativeKineticEnergySqrd += m_kineticEnergy * m_kineticEnergy;
    m_cumulativeExternalEnergySqrd += m_externalEnergy * m_externalEnergy;
    m_cumulativeInteractionEnergySqrd += m_interactionEnergy * m_interactionEnergy;

    if (stepNumber > (m_numberOfStepsPerBatch * (m_iter % m_numberOfBatches))
        && stepNumber < (m_numberOfStepsPerBatch * (m_iter % m_numberOfBatches + 1))) {
        m_cumulativeGradients += m_instantGradients;
        m_cumulativeGradientsE += m_instantGradients * m_instantEnergy;
    }
    if (acceptedStep) {
        m_acceptence += 1;
    }
}


/* ----------------------------------------------------------------------------
  Set number of cycles for each iteration
---------------------------------------------------------------------------- */

void Sampler::setNumberOfSteps(int numberOfStepsWOEqui,
                               int totalNumberOfStepsWOEqui,
                               int totalNumberOfStepsWEqui)
{
    m_stepsWOEqui = numberOfStepsWOEqui;
    m_totalStepsWOEqui = totalNumberOfStepsWOEqui;
    m_totalStepsWEqui = totalNumberOfStepsWEqui;
    m_numberOfStepsPerBatch = int(m_totalStepsWOEqui / m_numberOfBatches);
}


/* ----------------------------------------------------------------------------
  Sum the expectation values from all processes
---------------------------------------------------------------------------- */
void Sampler::computeTotals()
{
    int parameterSlots = int(m_numberOfElements * m_maxParameters);
    m_totalCumulativeGradients = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxParameters);
    m_totalCumulativeGradientsE = Eigen::MatrixXd::Zero(m_numberOfElements, m_maxParameters);
    MPI_Reduce(&m_acceptence, &m_totalAcceptence, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&m_cumulativeKineticEnergy,
               &m_totalCumulativeKineticEnergy,
               1,
               MPI_DOUBLE,
               MPI_SUM,
               0,
               MPI_COMM_WORLD);
    MPI_Reduce(&m_cumulativeExternalEnergy,
               &m_totalCumulativeExternalEnergy,
               1,
               MPI_DOUBLE,
               MPI_SUM,
               0,
               MPI_COMM_WORLD);
    MPI_Reduce(&m_cumulativeInteractionEnergy,
               &m_totalCumulativeInteractionEnergy,
               1,
               MPI_DOUBLE,
               MPI_SUM,
               0,
               MPI_COMM_WORLD);
    MPI_Reduce(&m_cumulativeEnergy,
               &m_totalCumulativeEnergy,
               1,
               MPI_DOUBLE,
               MPI_SUM,
               0,
               MPI_COMM_WORLD);
    MPI_Reduce(&m_cumulativeEnergySqrd,
               &m_totalCumulativeEnergySqrd,
               1,
               MPI_DOUBLE,
               MPI_SUM,
               0,
               MPI_COMM_WORLD);
    MPI_Reduce(&m_cumulativeKineticEnergySqrd,
               &m_totalCumulativeKineticEnergySqrd,
               1,
               MPI_DOUBLE,
               MPI_SUM,
               0,
               MPI_COMM_WORLD);
    MPI_Reduce(&m_cumulativeExternalEnergySqrd,
               &m_totalCumulativeExternalEnergySqrd,
               1,
               MPI_DOUBLE,
               MPI_SUM,
               0,
               MPI_COMM_WORLD);
    MPI_Reduce(&m_cumulativeInteractionEnergySqrd,
               &m_totalCumulativeInteractionEnergySqrd,
               1,
               MPI_DOUBLE,
               MPI_SUM,
               0,
               MPI_COMM_WORLD);
    MPI_Reduce(m_cumulativeGradients.data(),
               m_totalCumulativeGradients.data(),
               parameterSlots,
               MPI_DOUBLE,
               MPI_SUM,
               0,
               MPI_COMM_WORLD);
    MPI_Reduce(m_cumulativeGradientsE.data(),
               m_totalCumulativeGradientsE.data(),
               parameterSlots,
               MPI_DOUBLE,
               MPI_SUM,
               0,
               MPI_COMM_WORLD);
}


/* ----------------------------------------------------------------------------
  Compute the average of values over all processes
---------------------------------------------------------------------------- */

void Sampler::computeAverages()
{
    m_averageKineticEnergy = m_totalCumulativeKineticEnergy / m_totalStepsWOEqui;
    m_averageExternalEnergy = m_totalCumulativeExternalEnergy / m_totalStepsWOEqui;
    m_averageInteractionEnergy = m_totalCumulativeInteractionEnergy / m_totalStepsWOEqui;
    m_averageEnergy = m_totalCumulativeEnergy / m_totalStepsWOEqui;

    m_averageKineticEnergySqrd = m_totalCumulativeKineticEnergySqrd / m_totalStepsWOEqui;
    m_averageExternalEnergySqrd = m_totalCumulativeExternalEnergySqrd / m_totalStepsWOEqui;
    m_averageInteractionEnergySqrd = m_totalCumulativeInteractionEnergySqrd / m_totalStepsWOEqui;
    m_averageEnergySqrd = m_totalCumulativeEnergySqrd / m_totalStepsWOEqui;

    m_variance = (m_averageEnergySqrd - m_averageEnergy * m_averageEnergy) / m_totalStepsWOEqui;
    m_varianceKinetic = (m_averageKineticEnergySqrd - m_averageKineticEnergy * m_averageKineticEnergy) / m_totalStepsWOEqui;
    m_varianceExternal = (m_averageExternalEnergySqrd - m_averageExternalEnergy * m_averageExternalEnergy) / m_totalStepsWOEqui;
    m_varianceInteraction = (m_averageInteractionEnergySqrd - m_averageInteractionEnergy * m_averageInteractionEnergy) / m_totalStepsWOEqui;

    m_stdError = sqrt(m_variance);
    m_stdErrorKinetic = sqrt(m_varianceKinetic);
    m_stdErrorExternal = sqrt(m_varianceExternal);
    m_stdErrorInteraction = sqrt(m_varianceInteraction);

    m_averageGradients = m_totalCumulativeGradients / m_numberOfStepsPerBatch;
    m_averageGradientsE = m_totalCumulativeGradientsE / m_numberOfStepsPerBatch;

    if (std::isnan(m_averageKineticEnergy)) {
        perror("Energy exploded, please decrease the learning rate");
        cout << " Kinetic energy        : " << m_averageKineticEnergy << endl;
        cout << " External energy       : " << m_averageExternalEnergy << endl;
        cout << " Interaction energy    : " << m_averageInteractionEnergy << endl;
        MPI_Abort(MPI_COMM_WORLD, 143);
    }
}


/* ----------------------------------------------------------------------------
  Print output to the terminal for an ordinary iteration
---------------------------------------------------------------------------- */
void Sampler::printOutputToTerminal(const int maxIter, const double time)
{
    m_iter += 1;
    //cout << endl;
    cout << std::fixed;
    cout << std::setprecision(8);
    cout << m_iter << "  "
         << m_averageEnergy << "  "
         << m_stdError << "  "
         << m_averageKineticEnergy << "  "
         << m_stdErrorKinetic << "  "
         << m_averageExternalEnergy << "  "
         << m_stdErrorExternal << "  "
         << m_averageInteractionEnergy << "  "
         << m_stdErrorInteraction << "  "
         << double(m_totalAcceptence) / m_totalStepsWOEqui << "  "
         << time << std::endl;

    //cout << "  -- System info: " << " -- " << endl;
    //cout << " Iteration progress      : " << m_iter << "/" << maxIter << endl;
    //cout << " Number of particles     : " << m_numberOfParticles << endl;
    //cout << " Number of dimensions    : " << m_numberOfDimensions << endl;
    //cout << " Number of processes     : " << m_numberOfProcesses << endl;
    //cout << " Number of parameters    : " << m_system->getTotalNumberOfParameters() << endl;
    //cout << " Oscillator frequency    : " << m_omega << endl;
    //cout << " Trial wave function     : " << m_trialWaveFunction << endl;
    //cout << " Hamiltonian             : " << m_hamiltonian << endl;
    //cout << " # Monte Carlo Cycles    : " << m_totalStepsWEqui << " (" << m_totalStepsWOEqui
    //     << " after burn-in)" << endl;
    //cout << endl;
    //cout << "  -- Results -- " << endl;
    //cout << " Total energy            : " << m_averageEnergy;
    //cout << "        (with STD = "           << m_stdError << ")" << endl;
    //cout << " Kinetic energy          : " << m_averageKineticEnergy;
    //cout << "        (with STD = "           << m_stdErrorKinetic << ")" << endl;
    //cout << " External energy         : " << m_averageExternalEnergy;
    //cout << "        (with STD = "           << m_stdErrorExternal << ")" << endl;
    //cout << " Interaction energy      : " << m_averageInteractionEnergy;
    //cout << "        (with STD = "           << m_stdErrorInteraction << ")" << endl;
    //cout << " Variance                : " << m_variance << endl;
    //cout << " STD                     : " << sqrt(m_variance) << endl;
    //cout << " Acceptence Ratio        : " << double(m_totalAcceptence) / m_totalStepsWOEqui << endl;
    //cout << " CPU Time                : " << time << endl;
    //cout << endl;
}


/* ----------------------------------------------------------------------------
  Print output to terminal for final iteration
---------------------------------------------------------------------------- */

void Sampler::printFinalOutputToTerminal(std::chrono::system_clock::time_point start)
{
    cout << endl;
    cout << std::fixed;
    cout << std::setprecision(10);
    cout << "               Final results" << endl;
    cout << "==============================================" << endl;
    cout << " Total energy            : " << m_averageEnergy;
    cout << " (with STD = " << m_stdError << ")" << endl;
    cout << " Kinetic energy          : " << m_averageKineticEnergy
         << " (with STD = " << m_stdErrorKinetic << ")" << endl;
    cout << " External energy         : " << m_averageExternalEnergy
         << " (with STD = " << m_stdErrorExternal << ")" << endl;
    cout << " Interaction energy      : " << m_averageInteractionEnergy
         << " (with STD = " << m_stdErrorInteraction << ")" << endl;
    //cout << " Variance                : " << m_variance << " (with MSE = " << m_mseVariance << ")"
    //     << endl;
    //cout << " STD                     : " << m_stdError << " (with MSE = " << m_mseSTD << ")" << endl;
    cout << " Acceptence Ratio        : " << double(m_totalAcceptence) / m_totalStepsWOEqui << endl;
    cout << endl;

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    cout << "Finished computation at " << std::ctime(&end_time)
         << "Elapsed time: " << elapsed_seconds.count() << "s\n";
}


/* ----------------------------------------------------------------------------
  Perform blocking, based on the method and code implemented by Marius Jonsson.
  Please see github.com/computative
---------------------------------------------------------------------------- */

void Sampler::doResampling()
{
    if (m_printInstantEnergyToFile) {
        appendInstantFiles(".dat");
        appendInstantFiles("_kin.dat");
        appendInstantFiles("_ext.dat");
        appendInstantFiles("_int.dat");
        std::vector<std::string> files;
        files.push_back(m_instantEnergyFileName);
        files.push_back(m_instantKineticEnergyFileName);
        files.push_back(m_instantExternalEnergyFileName);
        files.push_back(m_instantInteractionEnergyFileName);

        std::ifstream infile(m_instantEnergyFileName.c_str());
        std::vector<double> instantEnergies;
        std::string line;
        while (std::getline(infile, line)) {
            instantEnergies.push_back(strtod(line.c_str(), nullptr));
        }
        Blocker block(instantEnergies);
        m_averageEnergy = block.mean;
        m_stdError = block.stdErr;
        m_variance = m_stdError * m_stdError;
        m_mseEnergy = block.mse_mean;
        m_mseSTD = block.mse_stdErr;
        m_mseVariance = m_mseSTD * m_mseSTD;
        //if(remove(m_instantEnergyFileName.c_str()) != 0) {
        //    perror( "Could not remove blocking file" );
        //}

        std::ifstream infile2(m_instantKineticEnergyFileName.c_str());
        std::vector<double> instantEnergies2;
        std::string line2;
        while (std::getline(infile2, line2)) {
            instantEnergies2.push_back(strtod(line2.c_str(), nullptr));
        }
        Blocker block2(instantEnergies2);
        m_averageKineticEnergy = block2.mean;
        m_stdErrorKinetic = block2.stdErr;
        m_varianceKinetic = m_stdErrorKinetic * m_stdErrorKinetic;
        m_mseEnergyKinetic = block2.mse_mean;
        m_mseSTDKinetic = block2.mse_stdErr;
        m_mseVarianceKinetic = m_mseSTDKinetic * m_mseSTDKinetic;
        if (remove(m_instantKineticEnergyFileName.c_str()) != 0) {
            perror("Could not remove blocking file");
        }

        std::ifstream infile3(m_instantExternalEnergyFileName.c_str());
        std::vector<double> instantEnergies3;
        std::string line3;
        while (std::getline(infile3, line3)) {
            instantEnergies3.push_back(strtod(line3.c_str(), nullptr));
        }
        Blocker block3(instantEnergies3);
        m_averageExternalEnergy = block3.mean;
        m_stdErrorExternal = block3.stdErr;
        m_varianceExternal = m_stdErrorExternal * m_stdErrorExternal;
        m_mseEnergyExternal = block3.mse_mean;
        m_mseSTDExternal = block3.mse_stdErr;
        m_mseVarianceExternal = m_mseSTDExternal * m_mseSTDExternal;
        if (remove(m_instantExternalEnergyFileName.c_str()) != 0) {
            perror("Could not remove blocking file");
        }

        std::ifstream infile4(m_instantInteractionEnergyFileName.c_str());
        std::vector<double> instantEnergies4;
        std::string line4;
        while (std::getline(infile4, line4)) {
            instantEnergies4.push_back(strtod(line4.c_str(), nullptr));
        }
        Blocker block4(instantEnergies4);
        m_averageInteractionEnergy = block4.mean;
        m_stdErrorInteraction = block4.stdErr;
        m_varianceInteraction = m_stdErrorInteraction * m_stdErrorInteraction;
        m_mseEnergyInteraction = block4.mse_mean;
        m_mseSTDInteraction = block4.mse_stdErr;
        m_mseVarianceInteraction = m_mseSTDInteraction * m_mseSTDInteraction;
        if (remove(m_instantInteractionEnergyFileName.c_str()) != 0) {
            perror("Could not remove blocking file");
        }
    }
}


/* ----------------------------------------------------------------------------
  Append the instant file from each process to a total instant file. To be used
  for blocking
---------------------------------------------------------------------------- */

void Sampler::appendInstantFiles(const std::string extension)
{
    std::string outfileName = m_path + std::to_string(m_instantNumber) + "_"
                              + std::to_string(m_rank) + extension;
    std::ofstream outfile(outfileName.c_str(), std::ios::out | std::ios::app);
    for (int i = 1; i < m_numberOfProcesses; i++) {
        std::string name = m_path + std::to_string(m_instantNumber) + "_" + std::to_string(i)
                           + extension;
        std::ifstream infile(name.c_str(), std::ios::in);
        if (!infile.is_open()) {
            perror("File not found");
        } else {
            outfile << infile.rdbuf();
        }
        if (remove(name.c_str()) != 0) {
            perror(" Could not remove blocking file");
        }
    }
}


/* ----------------------------------------------------------------------------
  Open the generated files to store local energy, electron densirt, weights
  and so on
---------------------------------------------------------------------------- */

void Sampler::openOutputFiles()
{
    //if (m_rank == 0) {
    //    m_instantNumber = m_system->getRandomNumberGenerator()->nextInt(1e6);
    //}
    MPI_Bcast(&m_instantNumber, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (m_printParametersToFile && m_rank == 0) {
        m_parameterFileName = m_path + "weights.dat";
        m_parameterFile.open(m_parameterFileName);
    }
    if (m_printEnergyToFile && m_rank == 0) {
        std::string averageEnergyFileName = m_path + "energy.dat";
        m_averageEnergyFile.open(averageEnergyFileName);
        std::string averageKineticEnergyFileName = m_path + "kinetic.dat";
        m_averageKineticEnergyFile.open(averageKineticEnergyFileName);
        std::string averageExternalEnergyFileName = m_path + "external.dat";
        m_averageExternalEnergyFile.open(averageExternalEnergyFileName);
        std::string averageInteractionEnergyFileName = m_path + "interaction.dat";
        m_averageInteractionEnergyFile.open(averageInteractionEnergyFileName);
    }
    if (m_computeOneBodyDensity && m_rank == 0) {
        std::string oneBodyFileName = m_path + "onebody.dat";
        m_oneBodyFile.open(oneBodyFileName);
    }
    if (m_computeOneBodyDensity2 && m_rank == 0) {
        std::string oneBodyFileName2 = m_path + "onebody2.dat";
        m_oneBodyFile2.open(oneBodyFileName2);
    }
    if (m_computeTwoBodyDensity && m_rank == 0) {
        std::string twoBodyFileName = m_path + "twobody.dat";
        m_twoBodyFile.open(twoBodyFileName);
    }
    if (m_printInstantEnergyToFile) {
        if (m_rank == 0) {
            m_instantNumber = m_system->getRandomNumberGenerator()->nextInt(unsigned(1e6));
        }
        MPI_Bcast(&m_instantNumber, 1, MPI_INT, 0, MPI_COMM_WORLD);
        m_instantEnergyFileName = m_path + std::to_string(m_instantNumber) + "_"
                                  + std::to_string(m_rank) + ".dat";
        m_instantKineticEnergyFileName = m_path + std::to_string(m_instantNumber) + "_"
                                         + std::to_string(m_rank) + "_kin.dat";
        m_instantExternalEnergyFileName = m_path + std::to_string(m_instantNumber) + "_"
                                          + std::to_string(m_rank) + "_ext.dat";
        m_instantInteractionEnergyFileName = m_path + std::to_string(m_instantNumber) + "_"
                                             + std::to_string(m_rank) + "_int.dat";
        m_instantEnergyFile.open(m_instantEnergyFileName);
        m_instantKineticEnergyFile.open(m_instantKineticEnergyFileName);
        m_instantExternalEnergyFile.open(m_instantExternalEnergyFileName);
        m_instantInteractionEnergyFile.open(m_instantInteractionEnergyFileName);
    }
}


/* ----------------------------------------------------------------------------
  Print total energy, kinetic energy, external potential and interaction energy
  to files
---------------------------------------------------------------------------- */

void Sampler::printEnergyToFile()
{
    if (m_printEnergyToFile && m_rank == 0) {
        m_averageEnergyFile << m_averageEnergy << endl;
        m_averageKineticEnergyFile << m_averageKineticEnergy << endl;
        m_averageExternalEnergyFile << m_averageExternalEnergy << endl;
        m_averageInteractionEnergyFile << m_averageInteractionEnergy << endl;
    }
}


/* ----------------------------------------------------------------------------
  Print parameters to file
---------------------------------------------------------------------------- */

void Sampler::printParametersToFile()
{
    if (m_printParametersToFile && m_rank == 0) {
        std::string parameterFileName = m_path + "weights.dat";
        m_parameterFile.open(parameterFileName);
        m_parameterFile << m_system->getWeights() << endl;
        m_parameterFile.close();
    }
}


/* ----------------------------------------------------------------------------
  Print radial one-body density to file
---------------------------------------------------------------------------- */

void Sampler::printOneBodyDensityToFile()
{
    if (m_computeOneBodyDensity) {
        m_totalParticlesPerBin = Eigen::VectorXi::Zero(m_numberOfBins);
        MPI_Reduce(m_particlesPerBin.data(),
                   m_totalParticlesPerBin.data(),
                   m_numberOfBins,
                   MPI_INT,
                   MPI_SUM,
                   0,
                   MPI_COMM_WORLD);
        if (m_rank == 0) {
            m_oneBodyFile << m_totalParticlesPerBin << endl;
        }
    }
}


/* ----------------------------------------------------------------------------
  Print spatial one-body density to file
---------------------------------------------------------------------------- */

void Sampler::printOneBodyDensity2ToFile()
{
    if (m_computeOneBodyDensity2) {
        m_totalDensityGrid = Eigen::MatrixXi::Zero(m_numberOfBins, m_numberOfBins);
        MPI_Reduce(m_densityGrid.data(),
                   m_totalDensityGrid.data(),
                   m_numberOfBins * m_numberOfBins,
                   MPI_INT,
                   MPI_SUM,
                   0,
                   MPI_COMM_WORLD);
        if (m_rank == 0) {
            m_oneBodyFile2 << m_totalDensityGrid << endl;
        }
    }
}


/* ----------------------------------------------------------------------------
  Print radial two-body density to file
----------------------------------------------------------------------------- */

void Sampler::printTwoBodyDensityToFile()
{
    if (m_computeTwoBodyDensity) {
        m_totalParticlesPerBinPairwise = Eigen::MatrixXi::Zero(m_numberOfBins, m_numberOfBins);
        MPI_Reduce(m_particlesPerBinPairwise.data(),
                   m_totalParticlesPerBinPairwise.data(),
                   m_numberOfBins * m_numberOfBins,
                   MPI_INT,
                   MPI_SUM,
                   0,
                   MPI_COMM_WORLD);
        if (m_rank == 0) {
            m_twoBodyFile << m_totalParticlesPerBinPairwise << endl;
        }
    }
}


/* ----------------------------------------------------------------------------
  Close all output files after printing
---------------------------------------------------------------------------- */

void Sampler::closeOutputFiles()
{
    if (m_averageEnergyFile.is_open()) {
        m_averageEnergyFile.close();
    }
    if (m_averageKineticEnergyFile.is_open()) {
        m_averageKineticEnergyFile.close();
    }
    if (m_averageExternalEnergyFile.is_open()) {
        m_averageExternalEnergyFile.close();
    }
    if (m_averageInteractionEnergyFile.is_open()) {
        m_averageInteractionEnergyFile.close();
    }
    if (m_oneBodyFile.is_open()) {
        m_oneBodyFile.close();
    }
    if (m_oneBodyFile2.is_open()) {
        m_oneBodyFile2.close();
    }
    if (m_twoBodyFile.is_open()) {
        m_twoBodyFile.close();
    }
    if (m_instantEnergyFile.is_open()) {
        m_instantEnergyFile.close();
    }
    if (m_parameterFile.is_open()) {
        m_parameterFile.close();
    }
}


/* ----------------------------------------------------------------------------
  Print instant energies to file for blocking
---------------------------------------------------------------------------- */

void Sampler::printInstantValuesToFile()
{
    if (m_printInstantEnergyToFile) {
        m_instantEnergyFile << m_instantEnergy << endl;
        m_instantKineticEnergyFile << m_kineticEnergy << endl;
        m_instantExternalEnergyFile << m_externalEnergy << endl;
        m_instantInteractionEnergyFile << m_interactionEnergy << endl;
    }
}


/* ----------------------------------------------------------------------------
  Compute radial one-body density
---------------------------------------------------------------------------- */

void Sampler::computeOneBodyDensity(const Eigen::VectorXd radialVector)
{
    if (m_computeOneBodyDensity) {
        for (int i_p = 0; i_p < m_numberOfParticles; i_p++) {
            if(radialVector(i_p) < m_maxRadius) {
                int bin = int(radialVector(i_p) / m_radialStep);
                if (bin < m_numberOfBins) {
                    m_particlesPerBin(bin)++;
                }
            }
        }
    }
}


/* ----------------------------------------------------------------------------
  Compute spatial one-body density
---------------------------------------------------------------------------- */

void Sampler::computeOneBodyDensity2(const Eigen::VectorXd positions)
{
    if (m_computeOneBodyDensity2) {
        for (int i_p = 0; i_p < m_numberOfParticles; i_p++) {
            Eigen::VectorXi indices = Eigen::VectorXi::Zero(m_numberOfDimensions);
            for (int i_d = 0; i_d < m_numberOfDimensions; i_d++) {
                int i = i_p * m_numberOfDimensions + i_d;
                indices(i_d) = int(positions(i) / (2 * m_radialStep));
            }
            if (abs(indices(0)) < m_numberOfBinsHalf && abs(indices(1)) < m_numberOfBinsHalf) {
                m_densityGrid(indices(1) + m_numberOfBinsHalf, indices(0) + m_numberOfBinsHalf)++;
            }
        }
    }
}


/* ----------------------------------------------------------------------------
  Compute radial two-body density
---------------------------------------------------------------------------- */

void Sampler::computeTwoBodyDensity(const Eigen::VectorXd radialVector)
{
    if (m_computeTwoBodyDensity) {
        for (int i_p = 0; i_p < m_numberOfParticles; i_p++) {
            int bin_i = int(radialVector(i_p) / m_radialStep) + 1;
            for (int j_p = i_p + 1; j_p < m_numberOfParticles; j_p++) {
                int bin_j = int(radialVector(j_p) / m_radialStep) + 1;
                if (bin_i < m_numberOfBins && bin_j < m_numberOfBins) {
                    m_particlesPerBinPairwise(bin_i, bin_j)++;
                }
            }
        }
    }
}
