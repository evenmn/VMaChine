#include "system.h"

int main(int argc, char **argv)
{
    // Define system
    System base;

    base.setNumberOfParticles(2);
    base.setNumberOfDimensions(2);
    base.setBurnInSteps(100);

    base.initializeFromConfig(argc, argv);
    base.runSimulation();

    return 0;
}
