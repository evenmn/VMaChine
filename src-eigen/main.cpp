#include "main.h"

int main(int argc, char *argv[])
{
    // Define system
    System *QD = new System();

    QD->setNumberOfParticles(6);
    QD->setNumberOfDimensions(2);
    QD->setFrequency(1.0);
    QD->setNumberOfMetropolisSteps(int(pow(2, 18)));

    QD->setLearningRate(0.05);
    QD->setStepLength(0.05);

    QD->checkResampling(true);
    QD->computeRadialOneBodyDensity();

    QD->setHamiltonian(new HarmonicOscillator(QD));

    // Define trial wave function ansatz
    //QD->setWaveFunctionElement(new Gaussian(QD));
    QD->setWaveFunctionElement(new SlaterDeterminant(QD));
    QD->setWaveFunctionElement(new SimpleJastrow(QD));
    QD->setWaveFunctionElement(new RBMProduct(QD));
    QD->setWaveFunctionElement(new RBMGaussian(QD));

    QD->setBasis(new Hermite(QD));

    QD->setNumberOfIterations(30);
    QD->runSimulation();
    return 0;
}
