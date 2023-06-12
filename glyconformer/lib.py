__author__ = "Isabell Louise Grothaus"
__license__ = ""
__version__ = ""
__email__ = "grothaus@uni-bremen.de"
__status__ = "Development"

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plumed
from scipy.signal import argrelextrema
import json
import matplotlib.ticker as ticker
import importlib.resources
import csv


class glyconformer:

    def __init__(self, inputfile, outputdir="./", length=None, glycantype=None ,angles=None, omega_angles=None, separator_index=None, separator=None, fepdir = None, order_max = None, order_min = None):
        
        #set vars
        self.inputfile = inputfile
        self.outputdir = outputdir 
        
        if glycantype is None:
            self.angles = angles
            self.omega_angles = omega_angles
            self.separator_index = separator_index
            self.separator = separator
            self.fepdir = fepdir
            self.order_max = order_max
            self.order_min = order_min
            self.maxima,self.minima = self.find_min_max(self.fepdir, self.angles, self.order_max, self.order_min)    
        else:
            self.glycantype = glycantype
            self.angles = self._readfeature("LIBRARY_GLYCANS.{}".format(self.glycantype),"angles.dat")
            self.omega_angles = self._readfeature("LIBRARY_GLYCANS.{}".format(self.glycantype),"omega_angles.dat")
            self.separator_index, self.separator = self._readseparator("LIBRARY_GLYCANS.{}".format(self.glycantype),"separator.dat")
            self.fepdir = "LIBRARY_GLYCANS/{}".format(self.glycantype)
            self.minima = self._readdict("LIBRARY_GLYCANS.{}".format(self.glycantype),"minima.dat")
            self.maxima = self._readdict("LIBRARY_GLYCANS.{}".format(self.glycantype),"maxima.dat")
            
        self.label = self.label_min(self.minima, self.angles, self.omega_angles)
        
        #read inputfile
        
        if length is None:
            self.colvar = self._readinputfile()
            self.length = len(self.colvar) 
        else:
            self.length = length 
            colvar = self._readinputfile()
            self.colvar = colvar.iloc[::round(colvar.shape[0]/self.length), :]
            
    def _readinputfile(self):
        try:
            colvar = plumed.read_as_pandas(self.inputfile)
        except:
            colvar = pd.read_csv(self.inputfile, delim_whitespace=True)
        finally:
            colvar = colvar[self.angles]
            
        return colvar
    
    def _readfeature(self, path, file):
        feature = importlib.resources.read_text(path,file)
        feature = feature.split()

        return feature

    def _readdict(self, path, file):
        """
        Function loading a dictionary from file.

        Parameters
        ----------
        path, file : str
            Name of the path + file to read in,
            including file extension 

        Returns
        -------
        dict
            Variable holding the dictionary
        """
    
        with importlib.resources.open_text(path,file) as f:
            data = f.read()
        dict = json.loads(data)
    
        return dict 
    
    def _readseparator(self, path, file):
        
        index = []
        sep = []
        with importlib.resources.open_text(path,file) as f:
            reader = csv.reader(f, delimiter=' ')
            for row in reader:
                index.append(int(row[0]))
                sep.append(str(row[1]))
        return index, sep
    
    def find_min_max(self, fepdir, features, order_max, order_min):
        """ 
        Finds minima and maxima of a 1D array read from file.

        Function reading in the free energy profiles (as arrays 
        from fes_*.dat files) of each feature (torsion angle) and 
        determining the maxima and minima of the array. The minima 
        array is altered in a way that an array with 2 maxima 
        corresponds to 3 minima, considering periodicity.

        Parameters
        ----------
        features : str
            List of torsion angle names (= column names in the 
            COLVAR file) in the correct order
        fepdir : str
            Directory from which to read the free energy files 
            (fes_*.dat) for the different torsion angles
        order_max: int
            How many points on each side of a maxima to use for 
            the comparison to consider comparator(n, n+x) to be True.
        order_min: int
            How many points on each side of a minima to use for 
            the comparison to consider comparator(n, n+x) to be True.

        Raises
        ------
        ValueError
            If free energy has more than 4 minimas, the algorithm can
            not handle it anymore.

        Returns
        -------
        maxmima_dict : dict
            Dictionary with features as keys and maxima points as values
        minima_dict : dict
            Dictionary with features as keys and minima points as values
        """

        maxima_dict = {}
        minima_dict = {}

        for f in features:
            profile = pd.read_csv("{}/fes_{}.dat".format(fepdir, f), delim_whitespace=True, names = ["x","y"])
            profmin = profile.index[profile['y'] == 0.0]
            profmin_value = profile.iloc[profmin[0],0]
            p1 = profile.iloc[:profmin[0],:]
            p2 = profile.iloc[profmin[0]:,:]
            p = p2.append([p1])
            p = p.to_numpy()

            maxima = argrelextrema(p[:,1], np.greater, order=order_max)
            minima = argrelextrema(p[:,1], np.less, order=order_min)

            if len(maxima[0]) == 0:
                maxima_dict["{}".format(f)] = [0]
            elif len(maxima[0]) == 1:
                maxima_dict["{}".format(f)] = [p[maxima[0][0],[0]].item()]
            elif len(maxima[0]) == 2:
                x = [p[maxima[0][0],[0]].item(), p[maxima[0][1],[0]].item()]
                x.sort()
                maxima_dict["{}".format(f)] = x
            elif len(maxima[0]) == 3:
                x = [p[maxima[0][0],[0]].item(), p[maxima[0][1],[0]].item(), p[maxima[0][2],[0]].item()]
                x.sort()
                maxima_dict["{}".format(f)] = x



            if len(minima[0]) == 0:
                minima_dict["{}".format(f)] = [profmin_value]

            elif len(minima[0]) == 1:   
                x = [profmin_value, p[minima[0][0],[0]].item()]
                x.sort()

                test = -3.5 <= x[0] <= maxima_dict["{}".format(f)][0]

                if test == True:                                  
                    x = [x[0],x[1],x[0]] 
                    minima_dict["{}".format(f)] = x

                else:
                    x = [x[1],x[0],x[1]] 
                    minima_dict["{}".format(f)] = x    


            elif len(minima[0]) == 2:
                x = [profmin_value, p[minima[0][0],[0]].item(), p[minima[0][1],[0]].item()]
                x.sort()

                test = -3.5 <= x[0] <= maxima_dict["{}".format(f)][0]

                if test == True:                                  
                    x = [x[0],x[1],x[2],x[0]] 
                    minima_dict["{}".format(f)] = x

                else:
                    x = [x[2],x[0],x[1],x[2]] 
                    minima_dict["{}".format(f)] = x                     

            elif len(minima[0]) >= 3:
                raise ValueError('More than 3 minima detected')

        return maxima_dict, minima_dict

    def label_min(self, minima, angles, omega_angles):
        """
        Function that labels minima by their corresponding IUPAC name.

        Reads in the created minima_dict, containing the position of 
        minimas along the torsion angles. These values are replaced by 
        capital letters originating from the IUPAC nomenclature for 
        torsion angles. 

        Parameters
        ----------
        cd : dict
            Minima_dict
        features : str
            List of torsion angle names (= column names in COLVAR file) 
            in the correct order
        f_omega : str
            List of omega torsion angle names in the correct order

        Returns
        -------
        label_dict : dict
            Dictionary with feature names as keys and sorted labels for each minima as values
        """

        tangle_dict = {
        " C ": (-0.5236, 0),
        " c ": (0, 0.5236),
        " G₊": (0.5236, 1.5708),
        " A₊": (1.5708, 2.6180),
        " T ": (2.6180, 3.5),
        " t ": (-3.5, -2.6180),
        " A₋": (-2.6180, -1.5708),
        " G₋": (-1.5708, -0.5236)
        }

        oangle_dict = {
        " gg": (-2.6180, 0),
        " gt": (0, 2.6180),
        " tg": (-3.5, -2.6180),
        " TG": (2.6180, 3.5)
        }

        name_dict = {}
        name_dict1 = {}
        name_dict2 = {}
        name_dict3 = {}


        for f in angles:
            if f in omega_angles:
                for key, value in oangle_dict.items():
                    if value[0] <= minima["{}".format(f)][0] < value[1]:
                        name_dict["{}".format(f)] = key
                    if minima["{}".format(f)][0] == 4:
                        name_dict["{}".format(f)] = 'none'      
            else:
                for key, value in tangle_dict.items():
                    if value[0] <= minima["{}".format(f)][0] < value[1]:
                        name_dict["{}".format(f)] = key
                    if minima["{}".format(f)][0] == 4:
                        name_dict["{}".format(f)] = 'none'

        for f in angles:
            if len(minima["{}".format(f)]) > 1:
                if f in omega_angles:
                    for key, value in oangle_dict.items():
                        if value[0] <= minima["{}".format(f)][1] < value[1]:                    
                            name_dict1["{}".format(f)] = key
                        if minima["{}".format(f)][1] == 4:
                            name_dict1["{}".format(f)] = 'none'        
                else:
                    for key,value in tangle_dict.items():
                        if value[0] <= minima["{}".format(f)][1] < value[1]:                    
                            name_dict1["{}".format(f)] = key
                        if minima["{}".format(f)][1] == 4:
                            name_dict1["{}".format(f)] = 'none'

        for f in angles:                
            if len(minima["{}".format(f)]) > 2:
                if f in omega_angles:
                    for key, value in oangle_dict.items():
                        if value[0] <= minima["{}".format(f)][2] < value[1]:                    
                            name_dict2["{}".format(f)] = key
                        if minima["{}".format(f)][2] == 4:
                            name_dict2["{}".format(f)] = 'none'                    
                else:
                    for key,value in tangle_dict.items():
                        if value[0] <= minima["{}".format(f)][2] < value[1]:                    
                            name_dict2["{}".format(f)] = key
                        if minima["{}".format(f)][2] == 4:
                            name_dict2["{}".format(f)] = 'none'

        for f in angles:                
            if len(minima["{}".format(f)]) > 3:
                if f in omega_angles:
                    for key, value in oangle_dict.items():
                        if value[0] <= minima["{}".format(f)][3] < value[1]:                    
                            name_dict3["{}".format(f)] = key
                else:
                    for key,value in tangle_dict.items():
                        if value[0] <= minima["{}".format(f)][3] < value[1]:                    
                            name_dict3["{}".format(f)] = key

        label_dict = self._combine_dict(name_dict,name_dict1,name_dict2,name_dict3)
        return label_dict
    
    def validate_fep(self):

        """
        Function that plots free energy profiles with annotated minima and maxima.

        Evaluation of how good the identification of minima and maxima worked.
        -------
        """

        y0 = 0
        y1 = 40
        C = np.arange(-0.5236, 0 + 0.001, 0.001)
        c = np.arange(0, 0.5236 + 0.001, 0.001)
        g = np.arange(0.5236, 1.5708 + 0.001, 0.001)
        a = np.arange(1.5708, 2.6180 + 0.001, 0.001)
        T = np.arange(2.6180, 3.5 + 0.001, 0.001)
        t = np.arange(-3.5, -2.6180 + 0.001, 0.001)
        A = np.arange(-2.6180, -1.5708 + 0.001, 0.001)
        G = np.arange(-1.5708, -0.5236 + 0.001, 0.001)

        co = "darkslategrey"

        fig, axs = plt.subplots(len(self.angles), figsize = (6,60), constrained_layout=True, )

        for i, f in zip(range(len(self.angles)),self.angles):
            profile = pd.read_csv("{}/fes_{}.dat".format(self.fepdir, f), delim_whitespace=True, names = ["x","y"])
            axs[i].set_title('{}'.format(f),fontsize = 15)
            axs[i].set_xlabel("Torsion angle [rad]", fontsize = 12)
            axs[i].set_ylabel("Free energy [kJ/mol]", fontsize = 12)
            axs[i].plot(profile.loc[:,"x"], profile.loc[:,"y"], c = "black", linewidth=2)
            axs[i].fill_between(C, y0, y1, alpha = 0.6, color = "paleturquoise", edgecolor = "none")
            axs[i].fill_between(c, y0, y1, alpha = 0.6, color = "paleturquoise", edgecolor = "none")
            axs[i].fill_between(a, y0, y1, alpha = 0.6, color = "darkturquoise", edgecolor = "none")
            axs[i].fill_between(g, y0, y1, alpha = 0.6, color = "aqua", edgecolor = "none")
            axs[i].fill_between(T, y0, y1, alpha = 0.6, color = "teal", edgecolor = "none")
            axs[i].fill_between(t, y0, y1, alpha = 0.6, color = "teal", edgecolor = "none")
            axs[i].fill_between(G, y0, y1, alpha = 0.6, color = "aqua", edgecolor = "none")
            axs[i].fill_between(A, y0, y1, alpha = 0.6, color = "darkturquoise", edgecolor = "none")
            if len(self.maxima["{}".format(f)]) == 1:
                axs[i].axvline(self.maxima["{}".format(f)][0], c = co, linestyle = "--")
                axs[i].text(self.minima["{}".format(f)][0],30, self.label["{}".format(f)][0], fontsize = 12, c = co, weight='bold')
            elif len(self.maxima["{}".format(f)]) == 2:
                axs[i].axvline(self.maxima["{}".format(f)][0], c = co, linestyle = "--")
                axs[i].axvline(self.maxima["{}".format(f)][1], c = co, linestyle = "--")
                axs[i].text(self.minima["{}".format(f)][0],30, self.label["{}".format(f)][0], fontsize = 12, c = co, weight='bold')
                axs[i].text(self.minima["{}".format(f)][1],30, self.label["{}".format(f)][1], fontsize = 12, c = co, weight='bold')
                axs[i].text(self.minima["{}".format(f)][2],30, self.label["{}".format(f)][2], fontsize = 12, c = co, weight='bold')

            elif len(self.maxima["{}".format(f)]) == 3:
                axs[i].axvline(self.maxima["{}".format(f)][0], c = co, linestyle = "--")
                axs[i].axvline(self.maxima["{}".format(f)][1], c = co, linestyle = "--")
                axs[i].axvline(self.maxima["{}".format(f)][2], c = co, linestyle = "--")
                axs[i].text(self.minima["{}".format(f)][0],30, self.label["{}".format(f)][0], fontsize = 12, c = co, weight='bold')
                axs[i].text(self.minima["{}".format(f)][1],30, self.label["{}".format(f)][1], fontsize = 12, c = co, weight='bold')
                axs[i].text(self.minima["{}".format(f)][2],30, self.label["{}".format(f)][2], fontsize = 12, c = co, weight='bold')
                axs[i].text(self.minima["{}".format(f)][3],30, self.label["{}".format(f)][3], fontsize = 12, c = co, weight='bold')

        plt.show()
    
    
    def run(self):
        
        binary, population, self.angles_separators, self.branches = self.create_binary(self.maxima, self.label, self.angles, self.separator_index, self.separator)
        self.perform_block_averages(binary, population, self.outputdir)
        
        return binary, population
        
    def _combine_dict(self, d1, d2, d3, d4):
        """
        Function that merges values from different 
        dictionaries having the same key into one.

        Parameters
        ----------
        d1,d2,d3,d4 : dict (hash table type)
            Up to 4 dictionaries 

        Returns
        -------
        dict
            1 merged dictionary
        """

        return {
            k: tuple(d[k] for d in (d1, d2, d3, d4) if k in d)
            for k in set(d1.keys()) | set(d2.keys()) | set(d3.keys()) | set(d4.keys())
        }

    def create_binary(self, maxima, label, angles, separator_index, separator):
        """
        Converts torsion angle values read from COVLAR file into IUPAC letters.

        Function that reads in a COLVAR file with torsion angle values stored 
        in each column and replaces each value by the corresponding IUPAC label, 
        which was previously defined in the label_dict.

        Parameters
        ----------
        maxima : dict
            Dictionary that contains the maxima for each feature
        label : dict
            Dictionary with features as keys and labels for each minima as values 
        angles : str
            List of torsion angle names (= column names) in the correct order
        separator: list ....
        loc3,loc4,loc5,loc6 : int 
            Location where to put the branch in the confromer string. Indicate
            the positon by an integer number
         bra3,bra4,bra5,bra6 : str
            Name of the separator to be used for the branch, e.g. "2──"

        Returns
        -------
        c : str
            Original COLVAR input as a pandas dataframe in which the torsion 
            angle values are replaced by labels
        ccf : mixed
            Pandas dataframe with a column for the conformer string and another 
            column for the occurance of the conformer (counts)
        angles : str
            Updated feature list with branch separators
        branches : str
            List of branches present in the processed glycan structure
        """

        colvar = self._readinputfile()
        colvar = colvar.iloc[::round(colvar.shape[0]/self.length), :]

        c = self._readinputfile()
        c = colvar.iloc[::round(colvar.shape[0]/self.length), :]

        for f in angles:
            if len(maxima["{}".format(f)]) == 1:

                c["{}".format(f)] = np.where((-3.5 <= colvar["{}".format(f)]) & (colvar["{}".format(f)] <= 3.5) , label["{}".format(f)][0], c["{}".format(f)])

            elif len(maxima["{}".format(f)]) == 2:

                c["{}".format(f)] = np.where((-3.5 <= colvar["{}".format(f)]) & (colvar["{}".format(f)] < maxima["{}".format(f)][0]) , label["{}".format(f)][0], c["{}".format(f)])
                c["{}".format(f)] = np.where((maxima["{}".format(f)][0] <= colvar["{}".format(f)]) & (colvar["{}".format(f)] < maxima["{}".format(f)][1]) , label["{}".format(f)][1], c["{}".format(f)])
                c["{}".format(f)] = np.where((maxima["{}".format(f)][1] <= colvar["{}".format(f)]) & (colvar["{}".format(f)] <= 3.5) , label["{}".format(f)][2], c["{}".format(f)])      

            elif len(maxima["{}".format(f)]) == 3:

                c["{}".format(f)] = np.where((-3.5 <= colvar["{}".format(f)]) & (colvar["{}".format(f)] < maxima["{}".format(f)][0]) , label["{}".format(f)][0], c["{}".format(f)])
                c["{}".format(f)] = np.where((maxima["{}".format(f)][0] <= colvar["{}".format(f)]) & (colvar["{}".format(f)] < maxima["{}".format(f)][1]) , label["{}".format(f)][1], c["{}".format(f)])
                c["{}".format(f)] = np.where((maxima["{}".format(f)][1] <= colvar["{}".format(f)]) & (colvar["{}".format(f)] < maxima["{}".format(f)][2]) , label["{}".format(f)][2], c["{}".format(f)])
                c["{}".format(f)] = np.where((maxima["{}".format(f)][2] <= colvar["{}".format(f)]) & (colvar["{}".format(f)] <= 3.5) , label["{}".format(f)][3], c["{}".format(f)])


        c = pd.DataFrame(c, columns = c.columns)
        c = c.replace(' t ', ' T ')
        c = c.replace(' TG', ' tg')

        branches = []

        if len(separator) == 1:
            c.insert(loc=separator_index[0], column='sep1', value=separator[0])
            branches.append(separator[0])
        elif len(separator) == 2:
            c.insert(loc=separator_index[0], column='sep1', value=separator[0])
            branches.append(separator[0])
            c.insert(loc=separator_index[1]+1, column='sep2', value=separator[1])
            branches.append(separator[1])
        elif len(separator) == 3:
            c.insert(loc=separator_index[0], column='sep1', value=separator[0])
            branches.append(separator[0])
            c.insert(loc=separator_index[1]+1, column='sep2', value=separator[1])
            branches.append(separator[1])
            c.insert(loc=separator_index[2]+2, column='sep3', value=separator[2])
            branches.append(separator[2])
        elif len(separator) == 4:
            c.insert(loc=separator_index[0], column='sep1', value=separator[0])
            branches.append(separator[0])
            c.insert(loc=separator_index[1]+1, column='sep2', value=separator[1])
            branches.append(separator[1])
            c.insert(loc=separator_index[2]+2, column='sep3', value=separator[2])
            branches.append(separator[2])
            c.insert(loc=separator_index[3]+3, column='sep4', value=separator[3])
            branches.append(separator[3])
        elif len(separator) == 5:
            c.insert(loc=separator_index[0], column='sep1', value=separator[0])
            branches.append(separator[0])
            c.insert(loc=separator_index[1]+1, column='sep2', value=separator[1])
            branches.append(separator[1])
            c.insert(loc=separator_index[2]+2, column='sep3', value=separator[2])
            branches.append(separator[2])
            c.insert(loc=separator_index[3]+3, column='sep4', value=separator[3])
            branches.append(separator[3])
            c.insert(loc=separator_index[4]+4, column='sep5', value=separator[4])
            branches.append(separator[4])
        elif len(separator) == 6:
            c.insert(loc=separator_index[0], column='sep1', value=separator[0])
            branches.append(separator[0])
            c.insert(loc=separator_index[1]+1, column='sep2', value=separator[1])
            branches.append(separator[1])
            c.insert(loc=separator_index[2]+2, column='sep3', value=separator[2])
            branches.append(separator[2])
            c.insert(loc=separator_index[3]+3, column='sep4', value=separator[3])
            branches.append(separator[3])
            c.insert(loc=separator_index[4]+4, column='sep5', value=separator[4])
            branches.append(separator[4])
            c.insert(loc=separator_index[5]+5, column='sep6', value=separator[5])
            branches.append(separator[5])
        elif len(separator) == 7:
            c.insert(loc=separator_index[0], column='sep1', value=separator[0])
            branches.append(separator[0])
            c.insert(loc=separator_index[1]+1, column='sep2', value=separator[1])
            branches.append(separator[1])
            c.insert(loc=separator_index[2]+2, column='sep3', value=separator[2])
            branches.append(separator[2])
            c.insert(loc=separator_index[3]+3, column='sep4', value=separator[3])
            branches.append(separator[3])
            c.insert(loc=separator_index[4]+4, column='sep5', value=separator[4])
            branches.append(separator[4])
            c.insert(loc=separator_index[5]+5, column='sep6', value=separator[5])
            branches.append(separator[5])
            c.insert(loc=separator_index[6]+6, column='sep7', value=separator[6])
            branches.append(separator[6])
        elif len(separator) == 8:
            c.insert(loc=separator_index[0], column='sep1', value=separator[0])
            branches.append(separator[0])
            c.insert(loc=separator_index[1]+1, column='sep2', value=separator[1])
            branches.append(separator[1])
            c.insert(loc=separator_index[2]+2, column='sep3', value=separator[2])
            branches.append(separator[2])
            c.insert(loc=separator_index[3]+3, column='sep4', value=separator[3])
            branches.append(separator[3])
            c.insert(loc=separator_index[4]+4, column='sep5', value=separator[4])
            branches.append(separator[4])
            c.insert(loc=separator_index[5]+5, column='sep6', value=separator[5])
            branches.append(separator[5])
            c.insert(loc=separator_index[6]+6, column='sep7', value=separator[6])
            branches.append(separator[6])
            c.insert(loc=separator_index[7]+7, column='sep8', value=separator[7])
            branches.append(separator[7])

        angles = list(c.columns)
        ccf = c.groupby(angles).size().to_frame(name = 'Count').reset_index()
        ccf['Conformer'] = ccf[ccf.columns[0:len(angles)]].apply(lambda x: ''.join(x.dropna().astype(str)),axis=1)
        ccf = ccf.drop(c.columns[0:len(angles)], axis=1)

        return c, ccf, angles, branches
    
    def perform_block_averages(self, c, ccf, outputdir):
        """
        Calculates probability of conformers for 10 separate blocks.

        Function that splits the torsion angle dataset in 10 
        blocks to do the conformer labelling and probability calculation 
        for each block.

        Parameters
        ---------- 
        c : str
            Binary COLVAR input file as a pandas dataframe in which the 
            torsion angle values are replaced by labels
        ccf : mixed
            Pandas dataframe with the first column for the conformer string 
            and the second column for the occurance of the conformer (counts)
        store_dir : str
            Directory name to store the conformer files in

        Returns
        -------
        """
        features = list(c.columns)
        ccf = ccf.set_index('Conformer')
        step = len(c)/10

        clear1 = np.arange(0,len(c), step)
        clear2 = np.arange(step, len(c) + step, step)
        clear1 = list(map(int, clear1))
        clear2 = list(map(int, clear2))

        for i,j,k in zip(range(1,11), clear1, clear2):
            c_short = c.iloc[j:k,:]
            cc_short = c_short.groupby(features).size().to_frame(name = 'Count').reset_index()
            norm = len(c_short)
            p = (cc_short.iloc[:,len(features)] / norm)
            cc_short["Probability"] = p
            cc_short['Conformer'] = cc_short[cc_short.columns[0:len(features)]].apply(lambda x: ''.join(x.dropna().astype(str)),axis=1)
            cc_short = cc_short.drop(cc_short.columns[0:len(features)], axis=1)
            # Join both dataframes taking 'Conformer' as key column for comparison
            cc_short = cc_short.set_index('Conformer')

            ccc_short = cc_short.join(ccf,how='outer',lsuffix='_block', rsuffix='_total')
            ccc_short = ccc_short.round(6)
            ccc_short = ccc_short.fillna(0.0)
            ccc_short = ccc_short.reset_index()

            # Print partially data sorted to file
            ccc_short.to_csv("{}/Cluster_conformer{}.dat".format(outputdir,i), sep= " ", header=True)
    
    def plot(self,
             threshold=2,
             ymax=100, 
             size=15, 
             colors=["#173c4d","#146b65","#4e9973","#a7c09f","#dfa790","#c76156","#9a2b4c","#600b4a"],
             dpi=300,
             file=None):
        """
        Plots a bar graph of the conformer density function.

        Function that plot the probability distribution of the 
        most occuring conformers. The limit defines the cutoff
        probability, where conformers are excluded from the 
        representation if below the cutoff.

        Parameters
        ----------
        input_dir: str
            Path to directory in which the "Cluster_conformer" 
            files from the block averaging step are stores
        limit: float
            Propability limit up to which conformers are included in the graph 
        branches: str
            List of branches present in the glycan
        features: str
            Updated features list
        ymax: float
            Maximum height of the y-axis (up to 100, as we plot probabilites), default is 100%

        Returns
        -------
        """

        hist1 = pd.read_csv("{}/Cluster_conformer1.dat".format(self.outputdir), delim_whitespace=True, header=0, names = ["Index", "Conformer", "Count_partial", "Prob", "Count_full"] )
        N, average, average2 = 1, hist1.iloc[:,3], hist1.iloc[:,3]*hist1.iloc[:,3]
        for i in range(2,11): 
            histn = pd.read_csv("{}/Cluster_conformer{}.dat".format(self.outputdir, i), delim_whitespace=True, header=0, names = ["Index", "Conformer", "Count_partial", "Prob", "Count_full"])
            N, average, average2 = N + 1, average + histn.iloc[:,3], average2 + histn.iloc[:,3]*histn.iloc[:,3]


        average = pd.DataFrame(average/N, columns=['Prob'])
        average = average.drop(average[average['Prob'] < threshold/100 ].index)
        average = average.sort_values(by='Prob',ascending=False)
        indexlist = average.index.tolist()

        hist = pd.read_csv("{}/Cluster_conformer1.dat".format(self.outputdir), names = ["Index", "Conformer", "Count_partial", "Prob", "Count_full"], sep = " ", dtype = str)
        hist = hist.filter(indexlist, axis = 0)
        namelist = hist['Conformer'].tolist()
        name_list = self._include_branch_conformer_string(namelist, self.branches, self.angles_separators)
        name_list = self._vertical_conformer_string(namelist)
        average2 = average2.to_frame()
        average2.columns = ["Error"]
        average2 = average2.filter(indexlist, axis = 0)

        # Final variances
        var = (N/(N-1))*( average2.Error / N - average.Prob*average.Prob ) 
        # Errors
        error = np.sqrt( var / N )

        pos_list = np.arange(len(name_list))

        colors = colors
        plt.rcParams['figure.dpi'] = dpi
        ax = plt.axes()
        ax.xaxis.set_major_locator(ticker.FixedLocator((pos_list)))
        ax.xaxis.set_major_formatter(ticker.FixedFormatter((name_list)))   
        bar = plt.bar(pos_list, average.Prob * 100, yerr = error * 100)
        for i in range(len(pos_list)):
            bar[i].set_color(colors[i])
        plt.ylim(0,ymax)
        plt.xticks(fontsize = size)
        plt.yticks(fontsize = size)
        plt.xlabel("Conformer", fontsize = size)
        plt.ylabel('Probability [%]', fontsize = size)
        
        if file is None:
            pass
        else:
            plt.savefig(file, bbox_inches='tight')
        plt.show()
    
    def _vertical_conformer_string(self, namelist):
        """
        Function that converts the conform string from a 
        horizontal format to a vertical one.

        Parameters
        ----------
        namelist : str
            List of conformer strings 

        Returns
        -------
        namelist_v : str
            Conformer strings as list in a vertical format
        """

        namelist_v = []
        for i in range(len(namelist)):
            res = ''
            for j,ele in enumerate(namelist[i]):
                if j in (range(2,1000,3)):
                    res += ele + "\n"
                else:
                    res += ele
            namelist_v += [res]
        return namelist_v

    def _include_branch_conformer_string(self, namelist, branches, features):
        """
        Function that inserts branch separators and replaces 
        labels by dots in the conform strings.

        Parameters
        ----------
        namelist : str
            List of conformer strings 
        branches : str
            List of branch names present in the glycan
        features : str
            List of features

        Returns
        -------
        namelist : str
            List of conformer strings with branch separators 
            and repetative labels exchange by dots
        """

        for i in range(1,len(namelist)):

            for j,k in zip(range(0,len(features)*3,3),range(3,len(features)*3+1,3)):
                if namelist[i][j:k] == namelist[0][j:k]:
                    newname = namelist[i][:j] + " • " + namelist[i][k:]
                    namelist[i] = newname
                else: 
                    pass 
            # Limitation: Can only replace branch separators which do occur max. twice
            for b in branches:
                index = namelist[0].find(b)
                newbranch = namelist[i][:index] + "───" + namelist[i][index+3:]
                namelist[i] = newbranch

            for b in branches:
                rindex = namelist[0].rfind(b)
                if rindex == "-1":
                    pass
                else:
                    rnewbranch = namelist[i][:rindex] + "───" + namelist[i][rindex+3:]
                    namelist[i] = rnewbranch

        return namelist




