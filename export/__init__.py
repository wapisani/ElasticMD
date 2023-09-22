# -*- coding: utf-8 -*-

"""
Created on Wed Nov 25 16:30:25 2020

@author: Dr. William A. Pisani
"""

def export2Csv(lists2Export,csvname,csvdir,header=''):
    """
    Takes a list of lists and exports them to a CSV file called "csvname"
    located in csvdir.

    Returns
    -------
    None.

    """
    import os
    import numpy as np
    os.chdir(csvdir)
    ncolumns = len(lists2Export) 
    nrows = len(lists2Export[0])
    
    # Add lists to NumPy 2D array
    data = np.zeros((nrows,ncolumns))
    for col in range(ncolumns):
        data[:,col] = lists2Export[col]
    
    # Save NumPy data
    np.savetxt(csvname,data,delimiter=",",header=header)
    print(f"{csvname} saved to {csvdir}")