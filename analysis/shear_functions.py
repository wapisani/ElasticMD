# -*- coding: utf-8 -*-
"""

@author: Dr. William A. Pisani

This file contains function definitions needed to analyze
shear simulation data from LAMMPS. These functions
are intended to be used with PySimpleGUI.
"""

import piecewise_regression
import matplotlib.pyplot as plt
from ..extract import log

def load_data(log_file):
    """
    This function loads in the shear data from a LAMMPS simulation.
    Prior to using this function, os.chdir() should be used to change to the 
    appropriate directory.

    Parameters
    ----------
    log_file : str
        The file name of the LAMMPS log file.

    Raises
    ------
    Exception
        If the file name is not in the required format, an exception is raised.
        The log file name should have Sh1, Sh2, or Sh3 in it. Sh=shear

    Returns
    -------
    strain_dir : str
        Direction of shear strain.
    primary_strain : numpy 1D array
        The true strain values in the direction of strain.
    primary_stress : numpy 1D array
        The true stress values in the direction of strain.
    """
    thermo = log(log_file)
    
    if 'Sh1' in log_file:
        strain_dir = 'xy'
        primary_strain, primary_stress = thermo.get(('v_etruexy','f_sxy_ave'),0)
        
        
    elif 'Sh2' in log_file:
        strain_dir = 'xz'
        primary_strain, primary_stress = thermo.get(('v_etruexz','f_sxz_ave'),0)
        
        
    elif 'Sh3' in log_file:
        strain_dir = 'yz'
        primary_strain, primary_stress = thermo.get(('v_etrueyz','f_syz_ave'),0)
        
        
    else:
        raise Exception(f"Something is wrong. Sh1, Sh2, Sh3 not found in {log_file}.\
                        Are you sure this is a shear sim log file?")

    return strain_dir, primary_strain, primary_stress



def compute_shear_modulus(primary_strain,primary_stress,n_breakpoints,data_range=(0,0)):
    """
    This function will use the piecewise-regression library to compute the
    shear modulus of the stress-strain data.

    Parameters
    ----------
    primary_strain : numpy 1D array
        The true strain values in the direction of strain.
    primary_stress : numpy 1D array
        The true stress values in the direction of strain.
    n_breakpoints : int
        Number of breakpoints to use in the segmented regression.
    data_range : tuple, optional
        The data range to be used in the analysis. The default is (0,0).

    Yields
    ------
    fit : piecewise-regression object
        The fit of the data.
    shearMod : float
        The shear modulus.
    shearMod_stderror : float
        The standard error of the shear modulus.
    yield_strain : float
        The yield strain.
    yield_strength : float
        The yield strength.

    """
    data_start = data_range[0]
    data_end = data_range[1]
    if data_end == 0:
        strain = primary_strain[data_start:]
        stress = primary_stress[data_start:]
    else:
        strain = primary_strain[data_start:data_end]
        stress = primary_stress[data_start:data_end]
        
    shearModFlag = True
    count = 0

    
    while shearModFlag:
        
        fit = piecewise_regression.Fit(strain,stress,n_breakpoints=n_breakpoints)
        
        results = fit.get_results()
        
        if results['converged'] == True:
            # Shear modulus will be the alpha of the first linear regression (alpha1)
            shearMod = results['estimates']['alpha1']['estimate']/1000 # in GPa
            shearMod_stderror = results['estimates']['alpha1']['se']/1000
            
            shearModFlag = False
            print("Shear modulus fitting has converged!")
            
            # Need to get yield strength, yield strain out as well
            # Apparently, the yield strain will not always be the first breakpoint (breakpoint1)
            breakpoint_names = [key for key in results['estimates'].keys() if 'breakpoint' in key]
            breakpoint_xvalues = [results['estimates'][key]['estimate'] for key in breakpoint_names]
            # The yield strain will be the breakpoint with the smallest x-value
            yield_breakpoint = breakpoint_names[breakpoint_xvalues.index(min(breakpoint_xvalues))]
            yield_strain = results['estimates'][yield_breakpoint]['estimate']
            c = results['estimates']['const']['estimate']
            yield_strength = shearMod*1000 * yield_strain + c
            
            return fit, shearMod, shearMod_stderror, yield_strain, yield_strength
            
        else:
            print("Shear modulus fitting did not converge. Trying again....")
            if count > 4:
                print("Try a different number of breakpoints or adjust the data range")
                shearMod = 0
                shearMod_stderror = 0
                yield_strain = 0
                yield_strength = 0
                return fit, shearMod, shearMod_stderror, yield_strain, yield_strength
            count += 1
            continue
    
    


def plot_shear(fit, strain_dir,shearMod, yield_strain, yield_strength):
    """
    This function will plot the shear modulus stress-strain data and the fit
    of that data.

    Parameters
    ----------
    fit : piecewise-regression object
        The fit of the elastic modulus data.
    strain_dir : str
        Direction of uniaxial strain.
    shearMod : float
        The shear modulus.
    yield_strain : float
        The yield strain.
    yield_strength : float
        The yield strength.

    Returns
    -------
    fig : matplotlib.pyplot figure object 
        The shear modulus figure object.
    ax : matplotlib.pyplot axis object
        The axes of the shear modulus figure object.

    """
    fit.plot_data(color="grey",s=20)
    fit.plot_fit(color="red",linewidth=4)
    fit.plot_breakpoints()
    fit.plot_breakpoint_confidence_intervals()
    plt.scatter(yield_strain,yield_strength,s=50,color='black',marker='^',zorder=4,label='$\sigma_{yield}$ = '+f'{yield_strength:0.2f} MPa' + '\n$\epsilon_{yield}$ = '+f'{yield_strain:0.4f}')
    plt.xlabel(f'True Strain in {strain_dir}')
    plt.ylabel(f'True Stress in {strain_dir}, MPa')
    plt.title(f"Shear modulus in {strain_dir} = {shearMod:0.2f} GPa")
    plt.legend(loc='lower right')
    fig = plt.gcf()
    ax = plt.gca()
    
    return fig, ax
    

def write_summaries(log_filename,strain_dir,shearMod,yield_strain,\
                    yield_strength,shear_fit, shearMod_stderror):
    """
    This function writes out two summary files of the shear modulus fit.
    The first is a text file of the shear modulus as well as a summary of the 
    fit. The second is an HTML file of the shear modulus as well as plots of 
    each of the fits. The W3.css file is needed in the same directory as the 
    HTML file for optimal viewing. The CSS file may be obtained here: 
    https://www.w3schools.com/w3css/4/w3.css.

    Parameters
    ----------
    log_filename : str
        Filename of the log file being analyzed.
    strain_dir : str
        The direction of strain.
    shearMod : float
        The predicted shear modulus.
    yield_strain : float
        The yield strain.
    yield_strength : float
        The yield strength.
    shear_fit : segmented-regression fit object
        The fit of the shear modulus stress-strain data.
    shearMod_stderror : float
        The standard error of the shear modulus fit.

    Yields
    ------
    None.

    """
    shear_plot_filename = log_filename.split('.log')[0]+f'_G{strain_dir}.pdf'
    
    
    # write out summaries of the fits to a file
    output_file = log_filename.split('.log')[0] + "_regression_summaries.txt"
    with open(output_file,'w') as handle:
        handle.write(f"This file contains the segmented regression results from a shear\
deformation simulation for {log_filename}.\n")
        handle.write("Shear modulus segmented regression results:\n")
        handle.write(f"Shear modulus in {strain_dir} = {shearMod:0.4f} GPa\n")
        handle.write(f"A check for shear modulus: yield strength divided by yield strain\
will probably be close to the value of shear modulus: {yield_strength/yield_strain:0.5f}.\n")
        handle.write("Shear modulus fit summary:\n")
        handle.write(shear_fit.summary())
       
    print(f"{output_file} has been saved!")  
    
    html_output = log_filename.split('.log')[0] + f'G{strain_dir}.html'
    with open(html_output,'w') as handle:
        handle.write('<html>\n')
        handle.write('<head><link rel="stylesheet" href="./w3.css"></head>\n')
        handle.write('<body>\n')
        handle.write(f'<h2>Analysis of {log_filename}</h2>\n')
        handle.write('<div class="w3-row">\n')
        handle.write(f"E_{strain_dir} = {shearMod:0.4f} +- {shearMod_stderror:0.4f} GPa\n")
        handle.write(f"<br />Yield strain = {yield_strain:0.6f}\n")
        handle.write(f"<br />Yield stress = {yield_strength:0.4f}\n")
        handle.write(f"<br /><br /><iframe height='100%' width='100%' src='./{shear_plot_filename}'></iframe>\n")
        handle.write("</div>\n</body>\n</html>")

    print(f"{html_output} has been saved!")
