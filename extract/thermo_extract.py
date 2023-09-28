# -*- coding: utf-8 -*-

"""

@author: Dr. William A. Pisani 

This module will read in a LAMMPS log file and extract all thermodynamic data.
"""

import numpy as np

class log:
    
    def __init__(self,logname,suppressOutput=False):
        self.keywords = [] # List of lists of thermo keywords (Step, Temp, Press, etc.), if different runs have different numbers of keywords
        self.nkeywords = [] # List of numbers of keywords
        self.keywordIndices = [] # List of dictionaries of thermo keyword-index pairs corresponding to each run
        self.headerIndices = [] # List of indices where thermo keywords are
        self.endRunIndices = [] # List of indices where "Loop time string is, marking the end of each run
        self.data = [] # List of 2d numpy arrays
        self.logname = logname # Name of log file
        self.extractKeywordsAndData()
        
        if suppressOutput == False:
            self.printHeaders()
    
    def __repr__(self):
        return f"{self.__class__.__name__} object containing thermodynamic data from {self.logname}"
    
    def printHeaders(self):
        """
        Prints thermodynamic headers in order.
        
        Returns
        -------
        None.

        """
        print(f"Data extraction from {self.logname} was successful!")
        print(f"{len(self.keywords)} thermodynamic section(s) found")
        print("\nThermodynamic keyword headers are as follows (index, header):")
        print("-------------------------------------------------------------")
        for index,header in enumerate(self.keywords):
            print(f"{index}\t{header}\n")
        if len(self.keywords) > 1:
            print("Please note that since more than one (1) section of non-identical thermodynamic data was found, you will need to specify which section of data you wish to extract.")
            print("For example, thermo.get(('Step','Temp','Press'),0) to get the step, temperature, and pressure from the first section of data")
    
    
    def extractKeywordsAndData(self):
        """
        Get all keyword/data sections from all LAMMPS runs in log file 

        Returns
        -------
        None.

        """
        with open(self.logname,'r') as logfile:
            logContents = logfile.read()
            
        splitLogContents = logContents.split('\n')
        
        for index,line in enumerate(splitLogContents):
            if line.find("Per MPI ") > -1 or line.find("Memory ") > -1: # Per MPI is for modern versions of LAMMPS, Memory is for older versions
                self.headerIndices.append(index+1) # Per MPI rank memory always occurs one line before the thermo keywords line
            elif line.find("Loop time ") > -1:
                self.endRunIndices.append(index)
        
        for index,headerIndex in enumerate(self.headerIndices):
            # Thermo keywords
            line = splitLogContents[headerIndex]
            headerLine = " ".join(line.split()).split(' ')
            
            # Raw Data
            try:
                start, stop = headerIndex+1, self.endRunIndices[index]
                rawData = splitLogContents[start:stop]
            except IndexError: # IndexError in self.endRunIndices[index]
                start = headerIndex+1
                rawData = splitLogContents[start:]
                print("From the log file, it appears that the simulation did not finish.\
                      This may impact your data analysis.")
            
            if headerLine not in self.keywords:
                self.keywords.append(headerLine)
                self.nkeywords.append(len(headerLine))
                
                keywordPairs = {}
                for index,keyword in enumerate(headerLine):
                    keywordPairs.update({keyword:index})
                self.keywordIndices.append(keywordPairs)
                
                # Convert raw data to numpy array
                npData = np.zeros((len(rawData),len(headerLine)))
                for i,dataLine in enumerate(rawData):
                    dataLine = " ".join(dataLine.split())
                    for j,value in enumerate(dataLine.split(' ')):
                        npData[i,j] = value
                self.data.append(npData)
                
            else: # If thermo header is identical to one already stored
                # Get index of first occurence of thermo header that is identical to the next header/data set to be stored
                
                firstOccurenceIndex = self.keywords.index(headerLine)
                # Convert raw data to numpy array
                npData = np.zeros((len(rawData),len(headerLine)))
                for i,dataLine in enumerate(rawData):
                    dataLine = " ".join(dataLine.split())
                    if 'WARNING' in dataLine:
                        continue
                    for j,value in enumerate(dataLine.split(' ')):
                        npData[i,j] = value
                # Add numpy array to numpy array of first occurence
                self.data[firstOccurenceIndex] = np.concatenate((self.data[firstOccurenceIndex],npData))
            
 
    def get(self,keys,index=0):
        """
        Retrieves the specified columns from the parsed data
        
        Parameters
        ----------
        keys : tuple
            Tuple of thermodynamic keywords. Common examples are "Step", "Temp", and "Press".
            Fixes, computes, and variables can also be extracted if given the appropriate term (e.g. "f_sxx_ave").
        index : int, optional
            Index for which section of data you wish to pull from. The default is 0.

        Raises
        ------
        Exception
            If no keywords are specified, an exception will be raised. At least one keyword must be specified.

        Returns
        -------
        properties : list
            List of 1d numpy arrays corresponding to the input tuple.

        """
        if len(keys) == 0:
            raise Exception("no keywords specified, you must specify at least one keyword (e.g. Step, Temp, etc)")
        
        if type(keys) == tuple:
            properties = []
            for key in keys:
                keyValue = self.keywordIndices[index][key]
                data = self.data[index][:,keyValue]
                properties.append(data)
        elif type(keys) == str:
            keyValue = self.keywordIndices[index][keys]
            data = self.data[index][:,keyValue]
            properties = data
        
        return properties
            
    def get_df(self, index=0):
        """
        Retrieves the thermodynamic data at the specified index and returns
        as a pandas dataframe

        Parameters
        ----------
        index : int, optional
            Index for which section of data you wish to pull from. The default is 0.

        Returns
        -------
        properties : pandas dataframe
            pandas dataframe containing the thermodynamic data of the specified
            section of data

        """
        
        import pandas as pd
        
        properties = pd.DataFrame(self.data[index],columns=self.keywords[0])
        
        return properties
    
    def export2Csv(lists2Export,csvname,csvdir,header=''):
        """
        Takes a list of lists and exports them to a CSV file called "csvname"
        located in csvdir.

        Parameters
        ----------
        lists2Export : list of lists
            List of lists to export to CSV.
        csvname : str
            Filename of CSV file to be written.
        csvdir : str
            Path to CSV file.
        header : str, optional
            If you want to add a comment to the CSV file. The default is ''.

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
    