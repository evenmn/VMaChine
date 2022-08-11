Random Number Generator
=======================

The code relies on the Mersenne-Twister random number generator, which has a fairly long period and a relatively efficient implementation in the C++ standard library. Other random number generators can be added, but currently Mersenne-Twister is the only one available. However, if you want an additional line in your configuration file, you can add

.. code-block::

    rng: MersenneTwister
