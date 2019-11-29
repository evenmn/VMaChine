#include "main.h"

int main(int argc, char *argv[])
{
    // Define system
    System *QD = new System();

    QD->setNumberOfParticles(2);
    QD->setNumberOfDimensions(2);
    QD->setFrequency(1.0);
    QD->setNumberOfMetropolisSteps(int(pow(2, 21)));

    QD->setLearningRate(0.5);
    QD->setStepLength(0.01);

    QD->setHamiltonian(new HarmonicOscillator(QD));

    // Define trial wave function ansatz
    QD->setWaveFunctionElement(new Gaussian(QD));
    QD->setWaveFunctionElement(new SlaterDeterminant(QD));
    QD->setWaveFunctionElement(new SimpleJastrow(QD));
    //QD->setWaveFunctionElement(new RBMProduct(QD));
    //QD->setWaveFunctionElement(new RBMGaussian(QD));

    QD->setBasis(new Hermite(QD));

    QD->setNumberOfIterations(1000);
    QD->runSimulation();
    return 0;
}
