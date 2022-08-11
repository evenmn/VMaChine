Hamiltionans
=============

Usually, the Hamiltionian is what we first want to specify when deciding on system. The Hamiltionian can be set using the rather intuitive :code:`hamiltionian`-keyword, ex:

.. code-block::

   hamiltonian: harmonicOscillator

The available Hamiltonians include :code:`harmonicOscillator`, :code:`doubleWell`, :code:`atomicNucleus` and :code:`ellipticalHarmonicOscillator`. More Hamiltonians can easily be added.

Parameters
-----------

Often, the various Hamiltonians depend on parameters. The preferred way of setting parameters is to specify them on the same line as the Hamiltonian itself is specified. For instance, the double well takes a parameter defining the distance between the two wells. It is specified in the following way:

.. code-block::

   hamiltonian: doubleWell 0.1

Similarly, the elliptical harmonic oscillator takes the ellipticity as argument:

.. code-block::

   hamiltonian: ellipticalHarmonicOscillator 0.5

However, for some Hamiltonians the parameter is crucial for other parts of the code than only the Hamiltonian. In those cases, the parameter is specified on a separate line to make it available globally. An example on this is the harmonic osciallator frequency,

.. code-block::

   omega: 1.0

which is used both by the circular and the elliptical harmonic oscillator. Another example is the atomic number that is used when simulating atoms,

.. code-block::

    hamiltonian: atomicNucleus
    atomicNumber: 2

Degrees of freedom
-------------------

When working with atoms, the degrees of freedom is well-defined by the atomic number as the system defaults to three dimensions. For quantum dot systems, on the other hand, the number of dimensions and particles have to be explicitly given by 

.. code-block::

   numParticles: 2
   numDimensions: 3

For closed-shell systems, make sure that the number of particles belong to the magic numbers :math:`N`:

.. math::

   N=s\binom{n+d}{d}

where :math:`s` is the number of spin configurations (2), :math:`d` is the number of dimensions and :math:`n` is an arbitrary integer. In two dimensions, the magic numbers are 2, 6, 12, 20 ...

Interaction
------------

Currently, the only available interaction styles are Coulomb interaction and no interaction. They can be specified using

.. code-block::

   interactionStyle: coulomb

or

.. code-block::

   interactionStyle: noInteraction

respectively. A screened Coulombic interaction may be implemented shortly.

Initialise
----------

The particles are usually initialised randomly after either a uniform or normal distribution. Set the initial state by either

.. code-block::

   initialState: randomUniform

or

.. code-block::

   initialState: randomNormal

The latter is the default.
