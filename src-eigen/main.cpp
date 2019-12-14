#include "main.h"

int main(int argc, char **argv)
{
    // Define system
    System *QD = new System();

    QD->setNumberOfParticles(2);
    QD->setNumberOfDimensions(2);
    QD->setFrequency(1.0);
    QD->setNumberOfMetropolisCycles(int(pow(2, 19)));

    QD->setLearningRate(0.3);
    QD->setStepLength(0.1);

    QD->doBlocking(false);

    QD->setHamiltonian(new HarmonicOscillator(QD));

    // Define trial wave function ansatz
    QD->setWaveFunctionElement(new PadeJastrow(QD));
    QD->setWaveFunctionElement(new RBMProduct(QD));
    QD->setWaveFunctionElement(new RBMGaussian(QD));

    QD->setConvergenceTools(4, 1e-4);
    QD->setAdaptiveStepTools(10, 3, 3);
    QD->setEquilibrationFraction(0.001);

    QD->setNumberOfIterations(100);
    QD->runSimulation();
    return 0;
}
