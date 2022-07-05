#include <memory>

#include "main.h"

int main(int argc, char **argv)
{
    // Define system
    System base;

    //base.setNumberOfParticles(2);
    //base.setNumberOfDimensions(2);
    //base.setEquilibrationFraction(0.001);

    base.initializeFromConfig(argc, argv);
    base.runSimulation();

    return 0;
}
