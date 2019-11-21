#include "main.h"

int main(int argc, char *argv[])
{
    // Define system
    System *QD = new System();

    QD->setNumberOfParticles(2);
    QD->setNumberOfDimensions(2);
    QD->setFrequency(1.0);
    QD->setInteraction(true);

    QD->setLearningRate(0.05);
    QD->setStepLength(0.005);
    QD->setNumberOfMetropolisSteps(int(pow(2, 20)));

    QD->setAdaptiveStepTools(true, 100, 4, 4);

    QD->setHamiltonian(new HarmonicOscillator(QD));

    // Define trial wave function
    QD->setBasis(new Hermite(QD));
    //QD->setWaveFunctionElement(new RBMGaussian(QD));
    QD->setWaveFunctionElement(new Gaussian(QD));
    QD->setWaveFunctionElement(new SlaterDeterminant(QD));
    QD->setWaveFunctionElement(new PadeJastrow(QD));
    //QD->setWaveFunctionElement(new DoubleProduct(QD));

    //QD->addDenseLayer(5, new ReLU(QD));

    QD->setNumberOfIterations(1000);
    QD->runSimulation();
    return 0;
}
