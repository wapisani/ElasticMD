# -*- coding: utf-8 -*-
"""

@author: Dr. William A. Pisani

This file contains a function to write out LAMMPS scripts for bulk modulus
simulations.

The files created by this function will not be ready to run because force field
parameters and a few other things need to be added or modified. This function
largely serves as a template that may be modified or incorporated into another
script generator.
"""


def write_bulk(script_filename,data_filename):
    """
    This function will write out a LAMMPS script that will run the system at
    standard temperature and pressure (300 K, 1.0 atm) for 500 picoseconds and 
    then ramp the pressure to 5000 atm while keeping the temperature the same.
    
    You will need to insert your force field parameters and ensure that the 
    units and timestep match your system.

    Parameters
    ----------
    script_filename : str
        Filename of bulk modulus LAMMPS script to be written.
    
    data_filename : str
        Filename of data file to be read into LAMMPS simulation.

    Returns
    -------
    None.

    """
    
    with open(script_filename,'w') as f:
        f.write("# This script was written with ElasticMD.\n")
        f.write("# This script will run at 300 K and 1.0 atm for 500 ps and then ramp to 5000 atm and 300 K for 500 ps. Both runs use NPT with aniso.\n")
        f.write("""#################################
### Initialization & Settings ###
#################################

# Debugging, I recommend having this for all simulations
echo both

units 		real
dimension	3
boundary	p p p
""")
        f.write(f"variable myid string {script_filename[:-3]}\n")
        f.write(f"variable data string {data_filename}\n")
        f.write("""log ${myid}.log.lammps

timestep        1 # 1 fs

##############################
### Force Field Parameters ###
##############################
# Insert your force field parameters here
# atom_style
# bond_style
# angle_style
# dihedral_style
# improper_style
# special_bonds
# pair_style
# pair_modify
# kspace_style
               
read_data ${data}

variable thermo_freq equal 10000
variable dump_freq equal 20*${thermo_freq}

##############################
########## Stresses ##########
##############################

compute         p all pressure thermo_temp
variable        sxx equal -0.101325*c_p[1] #in MPa
variable        syy equal -0.101325*c_p[2] #in MPa
variable        szz equal -0.101325*c_p[3] #in MPa
fix             sxx_ave all ave/time 1 ${thermo_freq} ${thermo_freq} v_sxx
fix             syy_ave all ave/time 1 ${thermo_freq} ${thermo_freq} v_syy
fix             szz_ave all ave/time 1 ${thermo_freq} ${thermo_freq} v_szz 

##############################
### Thermo & Dump Settings ###
##############################

thermo ${thermo_freq}
thermo_style custom step temp press etotal ke pe epair ebond eangle edihed eimp elong lx ly lz pxx pyy pzz v_sxx v_syy v_szz f_sxx_ave f_syy_ave f_szz_ave vol density

dump 1 all custom/gz ${dump_freq} ${myid}.lammpstrj.gz x y z type id 

####################
### Run Settings ###
####################
neigh_modify    delay 0 every 1 check yes one 2500 page 100000
fix	1 all npt temp 300 300 1000 aniso 1.0 1.0 10000
run 100000 # 100 ps
run 100000 # 100 ps
run 100000 # 100 ps
run 100000 # 100 ps
run 100000 # 100 ps
unfix 1

fix 1 all npt temp 300 300 1000 aniso 5000.0 5000.0 10000
run 100000 # 100 ps
run 100000 # 100 ps
run 100000 # 100 ps
run 100000 # 100 ps
run 100000 # 100 ps
""")
    print(f"{script_filename} written successfully!")
 
 
 
 