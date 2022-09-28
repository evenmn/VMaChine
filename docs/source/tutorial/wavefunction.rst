Wave function elements
=========================

The code supports a range of wave function elements that can be combined arbitrarily. The final trial wave function is the product of all the wave function elements specified. A wave function element is added using the :code:`waveFunctionElement`-keyword, for example:

.. code-block::

   waveFunctionElement: gaussian

The available wave function elements are :code:`gaussian`, :code:`slaterDeterminant`, :code:`padeJastrow`, :code:`simpleJastrow`, :code:`RBMGaussian`, :code:`RBMProduct`, :code:`hydrogenLike` and :code:`hardCoreJastrow`. 

Restricted Boltzmann machine (RBM)
-----------------------------------

The restricted Boltzmann machine (RBM) consists of the two wave function elements :code:`RBMGaussian` and :code:`RBMProduct`. They can be set independently using:

.. code-block::

   waveFunctionElement: RBMGaussian
   waveFunctionElement: RBMProduct

or using a predefined trial wave function:

.. code-block::

    waveFunction: RBM

This latter method will also include the Slater determinant. The number of hidden nodes is set using 

.. code-block::

   numHiddenNodes: 2

Sorting of inputs is a common method to ensure that the artificial neural network is symmetric (such that the entire wave function becomes antisymmetric). Apply this by setting

.. code-block::

   sorting: true


Slater determinant and basis set
----------------------------------

The Slater determinant should be included to keep the wave function anti-symmetric for many-fermion systems. It is added in the natural way,

.. code-block::

   waveFunctionElement: slaterDeterminant

The Slater determinant relies on a basis set, which has to be chosen carefully for every single system. For circular quantum dots, the Hermite polynomials is usually a wise choice as it is the analytical basis set for a non-interacting system. They can be set using

.. code-block::

   basis: hermite

The Hermite polynomials are hard-coded up to 17th degree for computational efficiency. Quantum dots with more than 110 particles may rely on recursive formulae for the Hermite polynomials. For atoms, the hydrogen-like orbitals are implemented. They can be set using

.. code-block::

   basis: hydrogenOrbital

The Slater determinant is usually among the most computationally intensive parts of the code to evaluate. Finding the gradient and the Laplacian of the Slater determinant with respect to the particle coordinates is also costly, as the closed forms depend on the inverse of the Slater matrix. We find these iteratively to keep the cost to a minimum. 

Parameter initialisation
-------------------------

The parameters in the wave function elements can be initialised randomly from a uniform distribution, normal distribution or all set to the same value. Parameters can also be initialised according to Xavier initialisation or from a parameter file (usually a checkpoint file). The default parameters are chosen according to the :code:`Automatize` class, which initialises the various wave function elements in a way that is known to work. To overwrite this, use the :code:`initialWeights`-keyword:

.. code-block::

   initialWeights: automatize

The available initializations are :code:`automatize`, :code:`randomuniform`, :code:`randomnormal`, :code:`constant`, :code:`xavier` and :code:`fromfile`.

The :code:`constant` initialisation takes the value as argument:

.. code-block::

   initialWeights: constant 1.0

The :code:`randomuniform` initialisation takes the maximum value of the range as argument:

.. code-block::

   initialWeights: randomuniform 0.2

The :code:`fromfile` initialisation takes the parameter file as argument:

.. code-block::

   initialWeights: fromfile weights.dat
