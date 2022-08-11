Optimization
=============

Supported optimization algorithms are the Adam optimizer (:code:`adam`), gradient descent (:code:`gd`) and stochastic gradient descent (:code:`sgd`). They are all based on the gradient only. Algorithms based on the Hessian matrix are not yet supported. Set optimization algorithms by

.. code-block::

    optimization: adam

Gradient descent and stochastic gradient descent both take the first momentum and decay rate as arguments. 

Learning rate
--------------
The learning rate is set with

.. code-block::

   learningRate: 0.01


Adaptive number of steps
-------------------------

The number of steps determine how accurate the estimation of expectation values is. This includes the ground state energy as well as the parameter gradients and the electron densities. To update the parameters correctly, the parameter gradients should be fairly accurate, but in the begining they do not need to be spot on. Therefore, it is often advantagous having an adaptive number of steps during the optimization. The last iteration is usually the production run, and should have a large number of steps to provide good estimates. Set adaptive steps by

.. code-block::

   applyAdaptiveSteps: true

and then set the range of the adaptive steps using

.. code-block::

   rangeOfAdaptiveSteps: 10

The number of additional steps is given by

.. code-block::

   additionalSteps: 2

where the value corresponds to the 2-exponent (since resampling usually is combined with additional steps). Setting 2 means increasing the number of steps by a factor 2^2. Similarly, one can set the additional steps for the very last iteration by

.. code-block::

   additionalStepsLastIter: 4

which is again compared to the original number of steps.
