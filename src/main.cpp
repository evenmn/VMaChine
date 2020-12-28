#include "main.h"

int main(int argc, char **argv)
{
    // Define system
    System *QD = new System();

    QD->setNumberOfParticles(2);
    QD->setNumberOfDimensions(2);
    QD->setEquilibrationFraction(0.001);

    QD->initializeFromConfig(argc, argv);
    QD->runSimulation();
    return 0;
}
