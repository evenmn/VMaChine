# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'VMaChine'
copyright = '2022, Even Marius Nordhagen'
author = 'Even Marius Nordhagen'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [ "breathe" ]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Breathe Configuration
breathe_projects = {"VMaChine": "../doxygen_out/xml/"}

breathe_projects_source = {
    "system": ("../../src", ["system.h"]),
    "sampler": ("../../src", ["sampler.h"]),
    "basis": ("../../src/Basis", ["basis.h", "hermite.h", "hydrogenorbital.h", "none.h"]),
    "hamiltonians": ("../../src/Hamiltonians", ["hamiltonian.h", "atomicnucleus.h", "doublewell.h", "ellipticalharmonicoscillator.h", "harmonicoscillator.h"]),
    "initialStates": ("../../src/InitialStates", ["initialstate.h", "randomuniform.h", "randomnormal.h"]),
    "initialWeights": ("../../src/InitialWeights", ["initialweights.h", "constant.h", "customized.h", "fromfile.h", "automatize.h", "randomnormal.h", "randomuniform.h", "xavier.h"]),
    "interaction": ("../../src/Interaction", ["interaction.h", "coulomb.h", "nointeraction.h"]),
    "metropolis": ("../../src/Metropolis", ["metropolis.h", "bruteforce.h", "importancesampling.h"]),
    "optimization": ("../../src/Optimization", ["optimization.h", "adam.h", "asgd.h", "barzilaiborwein.h", "gradientdescent.h", "sgd.h"]),
    "rng": ("../../src/RNG", ["rng.h", "mersennetwister.h"]),
    "wavefunctions": ("../../src/WaveFunctions", ["wavefunction.h", "gaussian.h", "hydrogenlike.h", "padejastrow.h", "rbmgaussian.h", "rbmproduct.h", "simplejastrow.h", "slaterdeterminant.h"]),
}

breathe_default_project = 'VMaChine'
