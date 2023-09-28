# -*- coding: utf-8 -*-
"""

@author: Dr. William A. Pisani

ElasticMD
========

A Python package to enable easy analysis of stress-strain curves produced
from LAMMPS MD simulations.

Provides
    1. A class to extract and export thermodynamic data from LAMMPS log files
    2. Functions to analyze bulk modulus simulations
    3. Functions to analyze stress-strain data and produce elastic moduli
       including Young's modulus, Poisson's ratio, and shear modulus
    4. GUIs to visually analyze stress-strain data, the recommended way
    5. Functions to generate LAMMPS input scripts for isotropic compression (bulk modulus), uniaxial tension (elastic modulus, Poisson's ratio), and shear (shear modulus)

"""

# Import modules
from . import extract
from . import analysis
from . import generate
