About VMaChine
================
:code:`VMaChine` is a package for running variational Monte Carlo (VMC) simulations optimized for machine learning-inspired trial wave functions. 

A configuration file has to be included when executing the code. A minimal configuration file might look like this:

.. code-block::

    # system
    numParticles: 6
    numDimensions: 2
    hamiltonian: harmonicOscillator
    omega: 1.0
    interactionStyle: coulomb

    # wave function
    waveFunctionElement: slaterDeterminant
    waveFunctionElement: gaussian
    waveFunctionElement: padeJastrow
    basis: hermite

    # simulation
    numIterations: 5000
    numSteps: 10000
    learningRate: 0.1
    stepLength: 0.05

If we name the configuration script :code:`input.in`, :code:`VMaChine` can be executed using

.. code-block::

   vmachine input.in
