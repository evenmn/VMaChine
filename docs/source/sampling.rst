Sampling
=========

For sampling the system we use Metropolis sampling, where the state probability is determined by :math:`P(\boldsymbol{r})=|\Psi(\boldsymbol{r})|^2`. There are two sampling algorithms available: brute force Metropolis (:code:`bruteForce`) and importance sampling (:code:`importanceSampling`). Brute force sampling moves particles in random directions, while importance sampling moves particles in the direction of the quantum force. One can therefore expect the acceptance ratio to be higher for importance sampling, which is the preferred and default sampling algorithm. Brute force sampling can be set using

.. code-block::

   sampling: bruteForce

The sampling algorithms are also associated with a maximum step length, which is set by

.. code-block::

   stepLength: 0.1

The number of steps per iteration is set using:

.. code-block::

   numSteps: 10000

If doing blocking resampling, make sure that the number of steps is an integer power of 2.

Burn-in
--------

Since the initial particle configuration might be an unlikely state, it is normal to run equilibriation steps before estimating expectation values. This period is called the burn-in time.
