"""
This module contains functions to process free energy
files produced by the PLUMED sum_hills command.

List of functions
rename_profiles
shift_to_zero
"""
import glob
import os

def _read_in_folder(folder, datatype):
    """
    Helper function that constructs a list of all files present in a directory
    
    Parameters
    ----------
    folder : str
        Path to the desired directory from which files
        should be checked out
    
    Returns
    -------
    dirList : str
        List containing path+file names
    """
    dirList = []
    for files in os.listdir(folder):
        if files.endswith(datatype):
            dirList.append(files)
    dirList.sort()
    if dirList == []:
        return None
    else:
        return dirList

def _get_CV(LINES):
    """
    Read the collective variable column of the .dat file
    
    Parameters
    ----------
    LINES : float
        Lines read in from a file
    
    Returns
    -------
    CV : float
        List of collective variables
    """
    CV = []
    for line in LINES:
        words = line.split()
        CV.append(float(words[0]))
    return CV

def _get_FE(LINES):
    """
    Read the free energy column of the .dat file
    
    Parameters
    ----------
    LINES : float
        Lines read in from a file
    
    Returns
    -------
    FE : float
        List of free energy values along the collective variable
    """
    FE = []
    for line in LINES:
        words = line.split()
        FE.append(float(words[1]))
    return FE

def _get_min(FILES):
    """
    Searches for the minimum value in the free energy columns 
    
    Parameters
    ----------
    FILES : str
        list of path+files that should be checked for the 
        minimum energy value
    
    Returns
    -------
    MIN : float
        lowest free energy value
    """
    FILE = open(FILES[-1], "r")
    LINES = FILE.readlines()[5:]
    FE = _get_FE(LINES)
    MIN = min(FE)
    return MIN


def _shift_FE(FE, MIN):
    """
    Shifts the energy column by the minimum free energy to align it to 0
    
    Parameters
    ----------
    FE : float
        List of free energy values along the collective variable
    MIN : float
        lowest free energy value
    
    Returns
    -------
    FE_NEW : float
        List of updated free energy values along the collective variable
    """
    FE_NEW = []
    for fe in FE:
        fe_shifted = fe - MIN
        FE_NEW.append(fe_shifted)
    return FE_NEW

    
def rename_profiles(path):
    """
    Function that reformats fes file names in order to get them read sorted.
    
    Instead of label like *1.dat, *2.dat ..., the code converts it 
    to *00001.dat, *00002.dat and so on.
    
    Parameters
    ----------
    path : str
        directory path including the folder with the fes files
    
    Returns
    -------
    """
    temp_dirList = _read_in_folder(path,'')
    dirList = []
    for folder in temp_dirList:
        dirList.append(path+folder)
    del temp_dirList
    for ii in range(len(dirList)):
        os.system('mv %sfes_%d.dat %sfes_%05d.dat' %(path, ii, path, ii))

def shift_to_zero(path):
    """
    Function that shifts all fes profiles in the path to zero.
    
    In order to get the free energy profiles aligned to zero,
    the code checks for the minimum free energy value and substracts
    it from all free energy profiles.
    
    Parameters
    ----------
    path : str
        directory path including the folder with the fes files
    
    Returns
    -------
    """
    FILES = glob.glob(path + "/fes*.dat")
    FILES.sort()
    MIN = _get_min(FILES)

    for item in FILES:
        FILE = open(item, "r")
        LINES = FILE.readlines()[5:]
        FILE.close()
        CV = _get_CV(LINES)
        FE = _get_FE(LINES)
        FE_NEW = _shift_FE(FE, MIN)
        OUTPUT = open(item, "w")
        for num in range(len(CV)):
            OUTPUT.write("{:<14.9F}{:<14.9F}\n".format(CV[num], FE_NEW[num]))
        OUTPUT.close()
