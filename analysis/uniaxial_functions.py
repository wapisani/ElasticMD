# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 12:04:12 2022

@author: Dr. Will Pisani

This file contains function definitions needed to analyze
uniaxial tensile simulation data from LAMMPS. These functions
are intended to be used with PySimpleGUI.

NOTE: Yield strength for compression is currently bugged. Will fix eventually.
"""

import piecewise_regression
import matplotlib.pyplot as plt
from ..extract import log
from sklearn.linear_model import LinearRegression

def load_data(log_file):
    """
    This function loads in the uniaxial tensile data from a LAMMPS simulation.
    Prior to using this function, os.chdir() should be used to change to the 
    appropriate directory. A CSV file can be used instead of a LAMMPS log file.

    Parameters
    ----------
    log_file : str
        The file name of the LAMMPS log file.

    Raises
    ------
    Exception
        If the file name is not in the required format, an exception is raised.
        The log file name should have YM1, YM2, or YM3 in it. YM=Young's Modulus

    Returns
    -------
    strain_dir : str
        Direction of uniaxial strain.
    primary_strain : numpy 1D array
        The true strain values in the direction of strain.
    primary_stress : numpy 1D array
        The true stress values in the direction of strain.
    secondary_strain1 : numpy 1D array
        The true strain values in the first of the two directions perpendicular
        to the primary strain.
    secondary_strain2 : numpy 1D array
        The true strain values in the second of the two directions perpendicular
        to the primary strain.
    Nu1_label : str
        The direction of the first Poisson's ratio. Used in plotting
    Nu2_label : str
        The direction of the second Poisson's ratio. Used in plotting

    """
    if '.csv' in log_file:
        import pandas as pd
        pd_flag = True
        thermo_pd = pd.read_csv(log_file)
    else:
        pd_flag = False
        thermo = log(log_file)
    
    if 'YM1' in log_file or 'Comp1' in log_file:
        strain_dir = 'x'
        Nu1_label = 'xy'
        Nu2_label = 'xz'
        if pd_flag:
            primary_strain = thermo_pd.v_etruex
            primary_stress = thermo_pd.f_sxx_ave
            secondary_strain1 = thermo_pd.v_etruey
            secondary_strain2 = thermo_pd.v_etruez
        else:
            primary_strain, primary_stress = thermo.get(('v_etruex','f_sxx_ave'),0)
            secondary_strain1, secondary_strain2 = thermo.get(('v_etruey','v_etruez'),0)
            
    elif 'YM2' in log_file or 'Comp2' in log_file:
        strain_dir = 'y'
        Nu1_label = 'yx'
        Nu2_label = 'yz'
        if pd_flag:
            primary_strain = thermo_pd.v_etruey
            primary_stress = thermo_pd.f_syy_ave
            secondary_strain1 = thermo_pd.v_etruex
            secondary_strain2 = thermo_pd.v_etruez
        else:
            primary_strain, primary_stress = thermo.get(('v_etruey','f_syy_ave'),0)
            secondary_strain1, secondary_strain2 = thermo.get(('v_etruex','v_etruez'),0)
        
        
    elif 'YM3' in log_file or 'Comp3' in log_file:
        strain_dir = 'z'
        Nu1_label = 'zx'
        Nu2_label = 'zy'
        if pd_flag:
            primary_strain = thermo_pd.v_etruez
            primary_stress = thermo_pd.f_szz_ave
            secondary_strain1 = thermo_pd.v_etruex
            secondary_strain2 = thermo_pd.v_etruey
        else:
            primary_strain, primary_stress = thermo.get(('v_etruez','f_szz_ave'),0)
            secondary_strain1, secondary_strain2 = thermo.get(('v_etruex','v_etruey'),0)
        
        
    else:
        raise Exception(f"Something is wrong. YM1, YM2, YM3 or Comp1, Comp2, Comp3 not found in {log_file}")

    return strain_dir, primary_strain, primary_stress, secondary_strain1, secondary_strain2, Nu1_label, Nu2_label



def compute_elastic_modulus(primary_strain,primary_stress,n_breakpoints,data_range=(0,0),compression_flag=False):
    """
    This function will use the piecewise-regression library to compute the
    elastic modulus of the stress-strain data.

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
    compression_flag : bool
        Flag used to determine if the simulation data being analyzed
        was created from uniaxial tension or uniaxial compression

    Yields
    ------
    fit : piecewise-regression object
        The fit of the data.
    elasticMod : float
        The elastic modulus.
    elasticMod_stderror : float
        The standard error of the elastic modulus.
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
        
    youngsModFlag = True
    ym_count = 0

    
    while youngsModFlag:
        
        if n_breakpoints == 0:
            model = LinearRegression()
            model.fit(strain.reshape((-1,1)),stress)
            r_squared = model.score(strain.reshape((-1,1)),stress)
            elasticMod = model.coef_[0]/1000
            yield_strain = 0
            yield_strength = 0
            return model, elasticMod, r_squared, yield_strain, yield_strength
        
        fit = piecewise_regression.Fit(strain,stress,n_breakpoints=n_breakpoints)
        
        results = fit.get_results()
        
        if results['converged'] == True:
            # Young's modulus will be the alpha of the first linear regression (alpha1)
            elasticMod = results['estimates']['alpha1']['estimate']/1000 # in GPa
            elasticMod_stderror = results['estimates']['alpha1']['se']/1000
            
            youngsModFlag = False
            print("Elastic modulus fitting has converged!")
            
            # Need to get yield strength, yield strain out as well
            # Apparently, the yield strain will not always be the first breakpoint (breakpoint1)
            breakpoint_names = [key for key in results['estimates'].keys() if 'breakpoint' in key]
            breakpoint_xvalues = [results['estimates'][key]['estimate'] for key in breakpoint_names]
            # The yield strain will be the breakpoint with the smallest x-value
            if compression_flag:
                yield_breakpoint = breakpoint_names[breakpoint_xvalues.index(max(breakpoint_xvalues))]    
                beta = f'beta{yield_breakpoint[-1]}' # Get number of breakpoint
                yield_strain = results['estimates'][yield_breakpoint]['estimate']
                c = results['estimates'][beta]['estimate']
                yield_strength = yield_strain * c # Yield strain times the beta seems to give the correct yield strength for compression, but not for tension
            else:
                yield_breakpoint = breakpoint_names[breakpoint_xvalues.index(min(breakpoint_xvalues))]
                yield_strain = results['estimates'][yield_breakpoint]['estimate']
                c = results['estimates']['const']['estimate']
                yield_strength = elasticMod*1000 * yield_strain + c
            
            return fit, elasticMod, elasticMod_stderror, yield_strain, yield_strength
            
        else:
            print("Elastic modulus fitting did not converge. Trying again....")
            if ym_count > 4:
                print("Try a different number of breakpoints or adjust the data range")
                elasticMod = 0
                elasticMod_stderror = 0
                yield_strain = 0
                yield_strength = 0
                return fit, elasticMod, elasticMod_stderror, yield_strain, yield_strength
            ym_count += 1
            continue
    
    


def compute_poisson(primary_strain, secondary_strain, n_breakpoints,Nu_label,data_range=(0,0)):
    """
    This function computes the Poisson's ratio 

    Parameters
    ----------
    primary_strain : numpy 1D array
        The true strain values in the direction of strain.
    secondary_strain : numpy 1D array
        The true strain values in the direction of the secondary strain.
    n_breakpoints : int
        Number of breakpoints to use in the segmented regression.
    Nu_label : str
        The direction of the Poisson's ratio. Used in plotting
    data_range : tuple, optional
        The data range to be used in the analysis. The default is (0,0).

    Returns
    -------
    fit : piecewise-regression object
        The fit of the data.
    Nu : float
        The Poisson's ratio.
    Nu_stderror : float
        The standard error of the fit of the Poisson's ratio.

    """
    data_start = data_range[0]
    data_end = data_range[1]
    if data_end == 0:
        p_strain = primary_strain[data_start:]
        s_strain = secondary_strain[data_start:]
    else:
        p_strain = primary_strain[data_start:data_end]
        s_strain = secondary_strain[data_start:data_end]
    
    if n_breakpoints == 0:
        model = LinearRegression()
        model.fit(p_strain.reshape((-1,1)),s_strain)
        r_squared = model.score(p_strain.reshape((-1,1)),s_strain)
        Nu = -1*model.coef_[0]
        return model, Nu, r_squared
    
    n_count = 0
    
    NuFlag = True
    while NuFlag:
        fit = piecewise_regression.Fit(p_strain,s_strain,n_breakpoints=n_breakpoints)
        
        results = fit.get_results()
        if results['converged'] == True:
            NuFlag = False
            # Nu will be the alpha of the first linear regression (alpha1),
            Nu = -results['estimates']['alpha1']['estimate']
            Nu_stderror = results['estimates']['alpha1']['se']
            print(f'Poisson ratio in {Nu_label} converged!')
            
            return fit, Nu, Nu_stderror
            
        else:
            if n_breakpoints == 0:
                print("Trying linear regression on specified data range...")
                model = LinearRegression()
                model.fit(p_strain.reshape((-1,1)),s_strain)
                r_squared = model.score(p_strain.reshape((-1,1)),s_strain)
                Nu = -1*model.coef_[0]
                return model, Nu, r_squared
            
            n_count += 1
            if n_count > 4:
                print("Try a different number of breakpoints or adjust the data range")
                Nu = 0
                Nu_stderror = 0
                return fit, Nu, Nu_stderror
            print(f'Poisson ratio in {Nu_label} did not converge. Trying again...')
        
    
    


def plot_elastic(fit, strain_dir, elasticMod, yield_strain, yield_strength,strain=None, stress=None):
    """
    This function will plot the elastic modulus stress-strain data and the fit
    of that data.

    Parameters
    ----------
    fit : piecewise-regression object
        The fit of the elastic modulus data.
    strain_dir : str
        Direction of uniaxial strain.
    elasticMod : float
        The elastic modulus.
    yield_strain : float
        The yield strain.
    yield_strength : float
        The yield strength.

    Returns
    -------
    fig : matplotlib.pyplot figure object 
        The elastic modulus figure object.
    ax : matplotlib.pyplot axis object
        The axes of the elastic modulus figure object.

    """
    if 'sklearn' in str(type(fit)):
        fig, ax = plt.subplots(figsize=(6.4,5),tight_layout=True)
        m = fit.coef_
        b = fit.intercept_
        y = m*strain + b
        ax.scatter(strain,stress,color="grey",s=20)
        ax.plot(strain,y,color="red",linewidth=4)
        ax.set_xlabel(f'True Strain in {strain_dir}')
        ax.set_ylabel(f'True Stress in {strain_dir}, MPa')
        ax.set_title(f"Young's modulus in {strain_dir} = {elasticMod:0.2f} GPa")
    else:

        fit.plot_data(color="grey",s=20)
        fit.plot_fit(color="red",linewidth=4)
        fit.plot_breakpoints()
        fit.plot_breakpoint_confidence_intervals()
        plt.scatter(yield_strain,yield_strength,s=50,color='black',marker='^',zorder=4,label='$\sigma_{yield}$ = '+f'{yield_strength:0.2f} MPa' + '\n$\epsilon_{yield}$ = '+f'{yield_strain:0.4f}')
        plt.xlabel(f'True Strain in {strain_dir}')
        plt.ylabel(f'True Stress in {strain_dir}, MPa')
        plt.title(f"Young's modulus in {strain_dir} = {elasticMod:0.2f} GPa")
        plt.legend(loc='lower right')
        fig = plt.gcf()
        ax = plt.gca()
    
    return fig, ax
    

def plot_poisson(fit, Nu_label,Nu,primary_strain=None,secondary_strain=None):
    """
    This function will plot the Poisson's ratio data along with the fit.

    Parameters
    ----------
    fit : piecewise-regression object
        The fit of the Poisson's ratio data.
    Nu1_label : str
        The direction of the Poisson's ratio. Used in plotting
    Nu1 : float
        The Poisson's ratio.
    primary_strain : 1D NumPy array
        Primary strain. Used only when the fit is a scikit-learn linear regression
    secondary_strain : 1D NumPy array
        Secondary strain. Used only when the fit is a scikit-learn linear regression

    Returns
    -------
    fig : matplotlib.pyplot figure object 
        The Poisson's ratio figure object.
    ax : matplotlib.pyplot axis object
        The axes of the Poisson's ratio figure object.

    """
    if 'sklearn' in str(type(fit)):
        fig, ax = plt.subplots(figsize=(6.4,5),tight_layout=True)
        m = fit.coef_
        b = fit.intercept_
        y = m*primary_strain + b
        ax.scatter(primary_strain,secondary_strain,color="grey",s=20)
        ax.plot(primary_strain,y,color="red",linewidth=4)
        ax.set_xlabel(f'True Strain in {Nu_label[0]}')
        ax.set_ylabel(f'True Strain in {Nu_label[1]}')
        ax.set_title(f"Nu_{Nu_label} = {Nu:0.5f}")
    else:
        fit.plot_data(color="grey",s=20)
        fit.plot_fit(color="red",linewidth=4)
        fit.plot_breakpoints()
        fit.plot_breakpoint_confidence_intervals()
        plt.xlabel(f'True Strain in {Nu_label[0]}')
        plt.ylabel(f'True Strain in {Nu_label[1]}')
        plt.title(f"Nu_{Nu_label} = {Nu:0.5f}")
        fig = plt.gcf()
        ax = plt.gca()
    
    return fig, ax
    

def write_summaries(log_filename,strain_dir,elasticMod,Nu1_label,Nu1,Nu2_label,\
                    Nu2,yield_strain,yield_strength,elastic_fit,Nu1_fit,Nu2_fit,\
                    ym_stderror,Nu1_stderror,Nu2_stderror):
    """
    This function writes out two summary files of the elastic modulus and 
    Poisson's ratio fits. The first is a text file of the elastic modulus and 
    Poisson's ratios as well as a summary of the fits of each. The second is 
    an HTML file of the elastic modulus and Poisson's ratio values as well 
    as plots of each of the fits. The W3.css file is needed in the same 
    directory as the HTML file for optimal viewing. The CSS file may be obtained
    here: https://www.w3schools.com/w3css/4/w3.css.

    Parameters
    ----------
    log_filename : str
        Filename of the log file being analyzed.
    strain_dir : str
        The direction of strain.
    elasticMod : float
        The predicted elastic modulus.
    Nu1_label : str
        The direction of the first Poisson's ratio.
    Nu1 : float
        The first Poisson's ratio.
    Nu2_label : str
        The direction of the second Poisson's ratio.
    Nu2 : float
        The second Poisson's ratio.
    yield_strain : float
        The yield strain.
    yield_strength : float
        The yield strength.
    elastic_fit : segmented-regression fit object
        The fit of the elastic modulus stress-strain data.
    Nu1_fit : segmented-regression fit object
        The fit of the primary strain and first secondary strain data.
    Nu2_fit : segmented-regression fit object
        The fit of the primary strain and second secondary strain data.
    ym_stderror : float
        The standard error of the elastic modulus fit.
    Nu1_stderror : float
        The standard error of the first Poisson's ratio fit.
    Nu2_stderror : float
        The standard error of the second Poisson's ratio fit.

    Yields
    ------
    None.

    """
    elastic_plot_filename = log_filename.split('.log')[0]+f'_stress_strain_{strain_dir}.pdf'
    Nu1_plot_filename = log_filename.split('.log')[0]+f'Nu{Nu1_label}.pdf'
    Nu2_plot_filename = log_filename.split('.log')[0]+f'Nu{Nu2_label}.pdf'
    
    # write out summaries of the fits to a file
    output_file = log_filename.split('.log')[0] + "_regression_summaries.txt"
    with open(output_file,'w') as handle:
        handle.write(f"This file contains the segmented regression results from a uniaxial\
                     tensile deformation simulation for {log_filename}.\n")
        handle.write("Elastic modulus segmented regression results:\n")
        handle.write(f"Elastic modulus in {strain_dir} = {elasticMod:0.4f} GPa\n")
        handle.write(f"Poisson's ratio in {Nu1_label} = {Nu1:0.7f}\n")
        handle.write(f"Poisson's ratio in {Nu2_label} = {Nu2:0.7f}\n\n")
        if yield_strain > 0:
            handle.write(f"A check for elastic modulus: yield strength divided by yield strain\
                         will probably be close to the value of elastic modulus: {yield_strength/yield_strain:0.5f}.\n")
        handle.write("Elastic modulus fit summary:\n")
        if 'sklearn' in str(type(elastic_fit)):
            handle.write("Linear regression from scikit-learn\n")
            handle.write(f"Rsquared = {ym_stderror:0.4f}\n")
            handle.write(f"y = {elastic_fit.coef_[0]/1000:0.4f} * x + {elastic_fit.intercept_:0.4f}\n")
        else:
            handle.write(elastic_fit.summary())
        handle.write(f"\nPoisson's ratio in {Nu1_label} fit summary:\n")
        if 'sklearn' in str(type(Nu1_fit)):
            handle.write("Linear regression from scikit-learn\n")
            handle.write(f"Rsquared = {Nu1_stderror:0.6f}\n")
            handle.write(f"y = {Nu1_fit.coef_[0]:0.6f} * x + {Nu1_fit.intercept_:0.6f}\n")
        else:
            handle.write(Nu1_fit.summary())
        handle.write(f"\nPoisson's ratio in {Nu2_label} fit summary:\n")
        if 'sklearn' in str(type(Nu2_fit)):
            handle.write("Linear regression from scikit-learn\n")
            handle.write(f"Rsquared = {Nu2_stderror:0.6f}\n")
            handle.write(f"y = {Nu2_fit.coef_[0]:0.6f} * x + {Nu2_fit.intercept_:0.6f}\n")
        else:
            handle.write(Nu2_fit.summary())
        
        
    print(f"{output_file} has been saved!")  
    
    html_output = log_filename.split('.log')[0] + f'E{strain_dir}_Nu{Nu1_label}_Nu{Nu2_label}.html'
    with open(html_output,'w') as handle:
        handle.write('<html>\n')
        handle.write('<head><link rel="stylesheet" href="./w3.css"></head>\n')
        handle.write('<body>\n')
        handle.write(f'<h2>Analysis of {log_filename}</h2>\n')
        handle.write('<div class="w3-row">\n')
        handle.write('<div class="w3-third">\n')
        if 'sklearn' in str(type(elastic_fit)):
            handle.write(f"E_{strain_dir} = {elasticMod:0.4f} GPa<br />Rsquared = {ym_stderror:0.4f} <br />\n")
        else:
            handle.write(f"E_{strain_dir} = {elasticMod:0.4f} +- {ym_stderror:0.4f} GPa\n")
            handle.write(f"<br />Yield strain = {yield_strain:0.6f}\n")
            handle.write(f"<br />Yield stress = {yield_strength:0.4f}\n")
        handle.write(f"<br /><br /><iframe height='100%' width='100%' src='./{elastic_plot_filename}'></iframe>\n")
        handle.write("</div>\n")
        handle.write('<div class="w3-third">\n')
        if 'sklearn' in str(type(Nu1_fit)):
            handle.write(f"Nu_{Nu1_label} = {Nu1:0.6f}<br />Rsquared = {Nu1_stderror:0.6f} <br />\n")
        else:
            handle.write(f"Nu_{Nu1_label} = {Nu1:0.6f} +- {Nu1_stderror:0.6f} <br /><br />\n")
        handle.write(f"<br /><br /><iframe height='100%' width='100%' src='./{Nu1_plot_filename}'></iframe>\n")
        handle.write("</div>\n")
        handle.write('<div class="w3-third">\n')
        if 'sklearn' in str(type(Nu2_fit)):
            handle.write(f"Nu_{Nu2_label} = {Nu2:0.6f}<br />Rsquared = {Nu2_stderror:0.6f} <br />\n")
        else:
            handle.write(f"Nu_{Nu2_label} = {Nu2:0.6f} +- {Nu2_stderror:0.6f} <br /><br />\n")
        handle.write(f"<br /><br /><iframe height='100%' width='100%' src='./{Nu2_plot_filename}'></iframe>\n")
        handle.write("</div>\n</div>\n</body>\n</html>")

    print(f"{html_output} has been saved!")






