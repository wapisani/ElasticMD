# -*- coding: utf-8 -*-
"""

@author: Dr. William A. Pisani

This script serves as an example of how to use ElasticMD to write LAMMPS 
scripts for predicting the elastic moduli (bulk, shear, elastic, Poisson).
"""


import os,sys
# Add the directory containing the ElasticMD folder to the path so it can be 
# found 
sys.path.insert(0,r'F:\Documents\Programming\Python')

import ElasticMD

directory = r'F:\Documents\Programming\Python\ElasticMD\examples\PA6'
os.chdir(directory)

# Let's try writing out a bulk modulus LAMMPS input script!
example_data_filename = r'PA6_10052023_IFFREq_180S1.dat.gz'
bulk_filename = example_data_filename.split('.')[0] + '_BK.in'

ElasticMD.generate.bulk(bulk_filename,example_data_filename)

# Now let's try writing out shear straining LAMMPS input scripts!
# Please note that "Sh1" "Sh2" "Sh3" will be added to the end of your filename
# so do not include the file extension in your filename. These are needed for 
# the shear_analysis_gui.py script to work properly.
shear_filename = example_data_filename.split('.')[0] + '_'
# We will reuse the data file from before
ElasticMD.generate.shear(shear_filename,example_data_filename,strain_rate=2e8,strain=0.2214028,directions=(1,2,3))

# Now let's try writing out uniaxial straining LAMMPS input scripts!
# Please note that "YM1" "YM2" "YM3" will be added to the end of your filename
# so do not include the file extension in your filename. These are needed for 
# the uniaxial_analysis_gui.py script to work properly.
elastic_filename = example_data_filename.split('.')[0] + '_'
# We will reuse the data file from before
ElasticMD.generate.elastic(elastic_filename,example_data_filename,strain_rate=2e8,strain=0.2214028,directions=(1,2,3))



