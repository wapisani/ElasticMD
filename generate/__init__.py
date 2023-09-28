# -*- coding: utf-8 -*-
"""
@author: Dr. William A. Pisani

This is the entry point for ElasticMD.generate
"""



def bulk(script_filename,data_filename):
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
atom_style
bond_style
angle_style
dihedral_style
improper_style
special_bonds
pair_style
pair_modify
kspace_style
               
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
 
 
def shear(script_filename,data_filename,strain_rate=2e8,strain=0.221402758160,directions = (1,2,3)):
    """
    This function will write out a LAMMPS script that will shear the system at 
    the specified strain rate until the specified level of engineering strain
    is achieved. Three scripts will be written out, each shearing a different 
    plane.

    Parameters
    ----------
    script_filename : str
        Filename of script, please note that '.in' should not be in filename.
    data_filename : str
        Data file to be read into LAMMPS.
    strain_rate : float, optional
        Strain rate in 1/s. The default is 2e8.
    strain : float, optional
        Engineering strain. The default is 0.221402758160. 
    directions : tuple, optional
        The directions to strain the system in. The default is (1,2,3).
        Can be specified to only generate one direction such as directions=(1)

    Returns
    -------
    None.

    """
    dir_dict = {1: 'xy', 2: 'xz', 3: 'yz'}
    
    for direction in directions:
        dir_str = dir_dict[direction]
        shear_filename = script_filename + f'Sh{direction}'
        with open(shear_filename+'.in','w') as f:
            f.write("# This script was written with ElasticMD.\n")
            f.write(f"# This script will strain the system in the {dir_str}-direction \
at a strain rate of {strain_rate} until an engineering strain of {strain} is achieved. \
Simulation will run at 300 K and 1.0 atm in the non-straining directions.\n")
            f.write("""
#################################
### Initialization & Settings ###
#################################

# Debugging
echo both

units 		real
dimension	3
boundary	p p p

""")
            f.write(f"variable myid string {shear_filename}\n")
            f.write(f"variable data string {data_filename}\n")
            f.write("log ${myid}.log.lammps\n")
            f.write("timestep 1 # 1 fs\n")
            f.write(f"variable strain equal {strain}\n")
            f.write(f"variable strain_rate_s equal {strain_rate}\n")
            f.write(f"variable dir equal {direction}\n")
            f.write("""# dir 1 means xy
# dir 2 means xz
# dir 3 means yz
##############################
### Force Field Parameters ###
##############################
# Insert your force field parameters here
atom_style	
bond_style 
angle_style 
dihedral_style 
improper_style 
special_bonds lj/coul 0 0 1
pair_style 
pair_modify mix sixthpower
kspace_style pppm 1e-6
read_data       ${data}

change_box all triclinic remap
kspace_style pppm 1e-6 # You MUST redefine kspace after changing the box to triclinic


variable        strain_rate_fs equal ${strain_rate_s}*1e-15 # in 1/fs
variable        totaltime equal ${strain}/${strain_rate_fs} # total time in femtoseconds
variable        steps equal $(round(v_totaltime/dt))
variable        quarter_steps equal $(floor(v_steps/4))
variable        eeng equal time*${strain_rate_fs}
variable        etrue equal ln(1+v_eeng)

variable thermo_freq equal $(round(v_steps/2000))
variable dump_freq equal $(round(v_steps/30))

##############################
########### Strains ##########
##############################
#cellgamma, cellbeta, and cellalpha refer to the crystallographic angles

variable	eengxy equal (PI/180*(90-cellgamma))
variable	eengxz equal (PI/180*(90-cellbeta))
variable	eengyz equal (PI/180*(90-cellalpha))
variable    etruexy equal ln(1+v_eengxy)
variable    etruexz equal ln(1+v_eengxz)
variable    etrueyz equal ln(1+v_eengyz)


##############################
########## Stresses ##########
##############################

compute         p all pressure thermo_temp
variable        sxx equal -0.101325*c_p[1] #in MPa
variable        syy equal -0.101325*c_p[2] #in MPa
variable        szz equal -0.101325*c_p[3] #in MPa
variable        sxy equal -0.101325*c_p[4] #in MPa
variable        sxz equal -0.101325*c_p[5] #in MPa
variable        syz equal -0.101325*c_p[6] #in MPa
fix             sxx_ave all ave/time 1 ${thermo_freq} ${thermo_freq} v_sxx
fix             syy_ave all ave/time 1 ${thermo_freq} ${thermo_freq} v_syy
fix             szz_ave all ave/time 1 ${thermo_freq} ${thermo_freq} v_szz
fix             sxy_ave all ave/time 1 ${thermo_freq} ${thermo_freq} v_sxy
fix             sxz_ave all ave/time 1 ${thermo_freq} ${thermo_freq} v_sxz
fix             syz_ave all ave/time 1 ${thermo_freq} ${thermo_freq} v_syz

##############################
### Thermo & Dump Settings ###
##############################

thermo ${thermo_freq}
thermo_style custom step temp press etotal ke pe epair ebond eangle edihed eimp elong lx ly lz vol density v_eengxy v_eengxz v_eengyz v_etruexy v_etruexz v_etrueyz v_sxx v_syy v_szz v_sxy v_sxz v_syz f_sxx_ave f_syy_ave f_szz_ave f_sxy_ave f_sxz_ave f_syz_ave v_dir
dump 1 all custom/gz ${dump_freq} ${myid}.lammpstrj.gz x y z type id
 
####################
### Run Settings ###
####################
neigh_modify    delay 0 every 1 check yes one 2500 page 100000
if "${dir} == 1" then &
  "fix          1 all npt temp 300 300 100 xz 1 1 1000 yz 1 1 1000 x 1 1 1000 y 1 1 1000 z 1 1 1000" &
  "fix          2 all deform 1 xy erate ${strain_rate_fs} flip no"
if  "${dir} == 2" then &
  "fix          1 all npt temp 300 300 100 xy 1 1 1000 yz 1 1 1000 x 1 1 1000 y 1 1 1000 z 1 1 1000" &
  "fix          2 all deform 1 xz erate ${strain_rate_fs} flip no"
if "${dir} == 3" then &
  "fix          1 all npt temp 300 300 100 xy 1 1 1000 xz 1 1 1000 x 1 1 1000 y 1 1 1000 z 1 1 1000" &
  "fix          2 all deform 1 yz erate ${strain_rate_fs} flip no"
# Breaking up the runs resets the pppm grid points which increases accuracy when box size changes
# Please note that the start and stop keywords are necessary so that the fix deform command deforms the box correctly over several runs
run		${quarter_steps} start 0 stop $(v_quarter_steps*4)
run		${quarter_steps} start 0 stop $(v_quarter_steps*4)
run		${quarter_steps} start 0 stop $(v_quarter_steps*4)
run		${quarter_steps} start 0 stop $(v_quarter_steps*4)
write_data ${myid}.dat
""")
            print(f"{shear_filename}.in written successfully!")
            
 
def elastic():
    pass