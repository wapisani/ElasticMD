# -*- coding: utf-8 -*-
"""
@author: Dr. William A. Pisani

"""

def bulk_modulus(steps,pressure,volume,plot_flag=False,plot_filename=''):
    """
    This function takes in three numpy lists, steps, pressure and volume, and returns
    the predicted bulk modulus according to the equation 
    K = - (Delta P)/(Delta V/V). I assume that the first half of the data is 
    from an NPT run at 1.0 atm and that the second half of the data is from an 
    NPT run at a much higher pressure (such as 5000 atm). I also assume that 
    the two runs are equal in length.

    Parameters
    ----------
    steps : numpy list
        Time step values from a bulk modulus simulation
    pressure : numpy list
        Pressure values in atmospheres (atm) from a bulk modulus simulation. 
    volume : numpy list
        Volume values from a bulk modulus simulation. Volume units cancel out 
        in the calculation so units do not matter.
    plot_flag : boolean
        If True, a plot of the pressure and volume will be created
    plot_filename : str
        If set, the pressure/volume plot will be saved to this filename.
        The user must set the filename correctly. pyplot supports a number of 
        image formats.
        
    Returns
    -------
    bulkMod : float
        The predicted bulk modulus in GPa.
    """
    from collections import defaultdict
    import numpy as np
    
    if len(pressure) != len(volume) != len(steps):
        raise Exception("Error: steps, pressure, and volume lists must be the same length")
    
    
    # Filter out all data from timesteps that repeat, e.g. a run ends on 100000
    # but starts at 100000 will have similar values
    # PaulMcG's solution from https://stackoverflow.com/questions/5419204/index-of-duplicates-items-in-a-python-list
    def listDuplicates(seq):
        tally = defaultdict(list)
        for i, item in enumerate(seq):
            tally[item].append(i)
        return [[key,locs] for key,locs in tally.items() if len(locs)>1]
        
    duplicates = listDuplicates(steps)
    indices_to_delete = [index[1][1] for index in duplicates]
    # Delete duplicate entries and have new arrays start at first non-zero step
    steps = np.delete(steps,indices_to_delete)[1:]
    pressure = np.delete(pressure,indices_to_delete)[1:]
    volume = np.delete(volume,indices_to_delete)[1:]
            
    midpoint = int(len(steps)/2)+1
    
    pressure1 = pressure[:midpoint]
    volume1 = volume[:midpoint]
    pressure2 = pressure[midpoint:]
    volume2 = volume[midpoint:]
    
    p1_avg = np.average(pressure1)
    v1_avg = np.average(volume1)
    p2_avg = np.average(pressure2)
    v2_avg = np.average(volume2)
    
    bulk_mod = -1*((p2_avg-p1_avg)/((v2_avg-v1_avg)/v1_avg)) *101325*10**-9 # in GPa 
    
    if plot_flag:
        import matplotlib.pyplot as plt
        # Plot pressure-volume curve
        fig, ax = plt.subplots()
        ax.scatter(volume,pressure,marker='o',facecolor='none',edgecolors='#244395')
        ax.scatter([v1_avg,v2_avg],[p1_avg,p2_avg],marker='^',facecolor='red',label="Avg press/vol")
        ax.set_xlabel('Volume, Angstrom$^{3}$')
        ax.set_ylabel('Pressure, atm')
        ax.annotate(str(np.round(p1_avg,1))+" atm",xy=(v1_avg-0.01*v1_avg,p1_avg+1000))
        ax.annotate(str(np.round(p2_avg,1))+" atm",xy=(v2_avg+0.005*v2_avg,p2_avg+500))
        ax.legend()
        ax.set_title("Bulk modulus = " + str(np.round(bulk_mod,4)) + " GPa")
        if plot_filename != '':
            fig.savefig(plot_filename)
        
    return bulk_mod 