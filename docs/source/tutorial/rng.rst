Random Number Generator
=======================

The code relies on the Mersenne-Twister random number generator, which has a fairly long period and a relatively efficient implementation in the C++ standard library. Other random number generators can be added, but currently Mersenne-Twister is the only one available. However, if you want an additional line in your configuration file, you can add

.. code-block::

    rng: MersenneTwister

Setting seed
------------

Sometimes it cna be convenient to set the RNG seed manually. The seed can be set to 12345 by

.. code-block::

   rng: MersenneTwister 12345

If seed is not manually set, it is picked by a random device. 
