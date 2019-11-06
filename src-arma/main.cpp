#include "main.h"

int main(int argc, char *argv[])
{
    // Define system
    System *QD = new System();

    QD->setNumberOfParticles(2);
    QD->setNumberOfDimensions(2);
    QD->setFrequency(1.0);
    QD->setInteraction(true);

    QD->setEquilibrationFraction(0.01);
    QD->setLearningRate(0.1);
    QD->setStepLength(0.05);
    QD->setNumberOfMetropolisSteps(arma::uword(pow(2, 20)));

    QD->setHamiltonian(new HarmonicOscillator(QD));

    // Define trial wave function
    QD->setBasis(new Hermite(QD));
    QD->setWaveFunctionElement(new Gaussian(QD));
    QD->setWaveFunctionElement(new SlaterDeterminant(QD));
    QD->setWaveFunctionElement(new PadeJastrow(QD));

    QD->setOptimization(new GradientDescent(QD, 0, 0.1));

    QD->runSimulation(1000);
    return 0;
}
