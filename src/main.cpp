#include "main.h"

int main(int argc, char **argv)
{
    // Define system
    System *QD = new System();

    QD->setNumberOfParticles(2);
    QD->setNumberOfDimensions(2);
    QD->setFrequency(1.0);
    QD->setNumberOfMetropolisCycles(int(pow(2, 19)));

    QD->setInteraction(false);
    QD->setLearningRate(0.01);
    QD->setStepLength(0.1);

    QD->dumpEnergyToFile();

    QD->setHamiltonian(new HarmonicOscillator(QD));

    // Define trial wave function ansatz
    QD->setWaveFunctionElement(new SlaterDeterminant(QD));
    //QD->setWaveFunctionElement(new RBMProduct(QD));
    //QD->setWaveFunctionElement(new RBMGaussian(QD));
    //QD->setWaveFunctionElement(new FNN(QD));
    QD->setWaveFunctionElement(new Gaussian(QD));

    //QD->addDenseLayer(25, new ReLU(QD));

    //QD->setInitialWeights(new Customized(QD));
    //QD->setOptimization(new GradientDescent(QD));

    QD->setConvergenceTools(4, 1e-4);
    QD->setAdaptiveStepTools(10, 2, 2);
    QD->setEquilibrationFraction(0.001);

    QD->setNumberOfIterations(100);

    QD->initializeFromConfig(argc, argv);
    QD->runSimulation();
    return 0;
}
