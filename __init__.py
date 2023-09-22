# -*- coding: utf-8 -*-
"""
ElasticMD
========

A Python package to enable easier analysis of stress-strain curves produced
from LAMMPS MD simulations.

Provides (work-in-progress)
    1. A class to extract thermodynamic data from LAMMPS log files
    2. A function to export user-specified data as CSV files
    3. Functions to analyze bulk modulus simulations
    4. Functions to analyze stress-strain data and produce elastic moduli
       including Young's modulus, Poisson's ratio, and shear modulus
    5. GUIs to visually analyze stress-strain data, the recommended way

@author: William A. Pisani, Ph.D.
"""

# Import modules
from . import extract
from . import export
from . import analysis
