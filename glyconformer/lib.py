__author__ = "Isabell Louise Grothaus"
__license__ = "GNU GENERAL PUBLIC LICENSE, Version 3"
__version__ = "1.0.1"
__email__ = "grothaus@uni-bremen.de"
__status__ = "Development"

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import plumed
from scipy.signal import argrelextrema
import json
import matplotlib.ticker as ticker
import importlib.resources
import csv
from sklearn.decomposition import PCA
from wpca import WPCA, EMPCA
import warnings
import os 

# Suppress FutureWarning messages
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

def _readinputfile(inputfile, angles):
    """
    Function reading a dataframe from file, using the self.inputfile.
    """
    try:
        colvar = plumed.read_as_pandas(inputfile)
    except:
        colvar = pd.read_csv(inputfile, delim_whitespace=True)
    finally:
        colvar = colvar[angles]
        
    return colvar

def _readfeature(path, file):
        """
        Function reading a list from file.

        Parameters
        ----------
        path, file : str
            Name of the path + file to read in,
            including file extension 

        Returns
        -------
        list
            Variable holding the list of strings
        """

        feature = importlib.resources.read_text(path, file)
        feature = feature.split()

        return feature

def _readdict(path, file):
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
    
        with importlib.resources.open_text(path, file) as f:
            data = f.read()
        dict = json.loads(data)
    
        return dict 
def _readseparator(path, file):
        """
        Function reading the separator information from file as lists.

        Parameters
        ----------
        path, file : str
            Name of the path + file to read in,
            including file extension 

        Returns
        -------
        index
        sep
            Variables holding the lists with the separator index and label
        """

        index = []
        sep = []
        with importlib.resources.open_text(path, file) as f:
            reader = csv.reader(f, delimiter=' ')
            for row in reader:
                index.append(int(row[0]))
                sep.append(str(row[1]))
        return index, sep

def _vertical_conformer_string(namelist):
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
def _include_branch_conformer_string(namelist, branches, features):
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

class Glyconformer():

    def __init__(self, inputfile=None, length=None, glycantype=None,
                 angles=None, omega_angles=None, separator_index=None, 
                 separator=None, fepdir=None, order_max=None, order_min=None, weights=None):
        # Instance variables
        """
        Initialize object's attributes, aka setting variables.

        Parameters
        ----------
        inputfile: str
            Name of the file to read in, including file extension 
        length: int
            Length to which the input file should be reduced
            Default: lenght is read from input file
        glycantype: str
            Select glycan type from "LIBRARY_GLYCANS"
        angles: list of str
            List of torsion angle names
        omega_angles: list of str
            List of omega torsion angle names
        separator_index: list of int
        separator: list of str
        fepdir: str
            Path to free energy profile files
        order_max: int
            How many points on each side of a maxima to use for 
            the comparison to consider comparator(n, n+x) to be True.
        order_min: int
            How many points on each side of a minima to use for 
            the comparison to consider comparator(n, n+x) to be True.

        Returns
        -------
        """
        # attribute reading mainly depends on if a glycantype is specified and 
        # whether information are read from the LIBRARY_GLYCANS or from 
        # user input variables. 
       
        if glycantype is None:
            self.inputfile = inputfile
            self.angles = angles
            self.omega_angles = omega_angles
            self.separator_index = separator_index
            self.separator = separator
            self.fepdir = fepdir
            self.order_max = order_max
            self.order_min = order_min
            self.maxima, self.minima = self._find_min_max() 
            self.weights = weights
        else:
            self.glycantype = glycantype
            self.inputfile = "TUTORIAL/{}_example/{}_angles.dat".format(self.glycantype,self.glycantype)

            self.angles = _readfeature("LIBRARY_GLYCANS.{}".format(self.glycantype), "angles.dat")
            self.omega_angles = _readfeature("LIBRARY_GLYCANS.{}".format(self.glycantype), "omega_angles.dat")
            self.separator_index, self.separator = _readseparator("LIBRARY_GLYCANS.{}".format(self.glycantype), "separator.dat")
            self.fepdir = "LIBRARY_GLYCANS/{}".format(self.glycantype)
            self.minima = _readdict("LIBRARY_GLYCANS.{}".format(self.glycantype), "minima.dat")
            self.maxima = _readdict("LIBRARY_GLYCANS.{}".format(self.glycantype), "maxima.dat")
            self.weights = None

        self.label = self._label_min()
        # read inputfile
        
        if length is None:
            self.colvar = _readinputfile(self.inputfile, self.angles)
            self.length = len(self.colvar) 
        else:
            self.length = length 
            colvar = _readinputfile(self.inputfile, self.angles)
            self.colvar = colvar.iloc[::round(colvar.shape[0]/self.length), :]

        self.binary, self.binary_compressed, self.count, self.angles_separator = self._create_binary()
    
    def _find_min_max(self):
        """ 
        Finds minima and maxima of a 1D array read from file.

        Function reading in the free energy profiles (as arrays 
        from fes_*.dat files) of each feature (torsion angle) and 
        determining the maxima and minima of the array. The minima 
        array is altered in a way that an array with 2 maxima 
        corresponds to 3 minima, considering periodicity.

        Parameters
        ----------
        
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

        for f in self.angles:
            profile = pd.read_csv("{}/fes_{}.dat".format(self.fepdir, f), 
                                  delim_whitespace=True, names=["x", "y"])
            profmin = profile.index[profile['y'] == 0.0]
            profmin_value = profile.iloc[profmin[0], 0]
            p1 = profile.iloc[:profmin[0], :]
            p2 = profile.iloc[profmin[0]:, :]
            p = pd.concat([p2, p1])
            p = p.to_numpy()

            maxima = argrelextrema(p[:, 1], np.greater, order=self.order_max)
            minima = argrelextrema(p[:, 1], np.less, order=self.order_min)

            if len(maxima[0]) == 0:
                maxima_dict["{}".format(f)] = [0]
            elif len(maxima[0]) == 1:
                maxima_dict["{}".format(f)] = [p[maxima[0][0], [0]].item()]
            elif len(maxima[0]) == 2:
                x = [p[maxima[0][0], [0]].item(), p[maxima[0][1], [0]].item()]
                x.sort()
                maxima_dict["{}".format(f)] = x
            elif len(maxima[0]) == 3:
                x = [p[maxima[0][0], [0]].item(), p[maxima[0][1], [0]].item(), p[maxima[0][2], [0]].item()]
                x.sort()
                maxima_dict["{}".format(f)] = x

            if len(minima[0]) == 0:
                minima_dict["{}".format(f)] = [profmin_value]

            elif len(minima[0]) == 1:   
                x = [profmin_value, p[minima[0][0], [0]].item()]
                x.sort()

                test = -3.5 <= x[0] <= maxima_dict["{}".format(f)][0]

                if test == True:                                  
                    x = [x[0], x[1], x[0]] 
                    minima_dict["{}".format(f)] = x

                else:
                    x = [x[1], x[0], x[1]] 
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

    def _label_min(self):
        """
        Function that labels minima by their corresponding IUPAC name.

        Reads in the created minima_dict, containing the position of 
        minimas along the torsion angles. These values are replaced by 
        capital letters originating from the IUPAC nomenclature for 
        torsion angles. 

        Parameters
        ----------

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
        name_dicts = [{} for _ in range(4)]  # Create a list of dictionaries for each iteration
        for f in self.angles:
            for i in range(min(4, len(self.minima[f]))):  # Iterate up to 4 times or the length of minima[f]
                angle_dict = oangle_dict if f in self.omega_angles else tangle_dict
                for key, value in angle_dict.items():
                    if value[0] <= self.minima[f][i] < value[1]:
                        name_dicts[i][f] = key
                    if self.minima[f][i] == 4:
                        name_dicts[i][f] = 'none'
        label_dict = {k: tuple(d.get(k, 'none') for d in name_dicts) for k in set().union(*name_dicts)}
        return label_dict

    def _create_binary(self):
        """
        Converts torsion angle values read from the input file into IUPAC letters.
    
        Returns
        -------
        binary: pd.DataFrame
            Original COLVAR input as a pandas dataframe in which the torsion 
            angle values are replaced by labels
        count : pd.DataFrame
            Pandas dataframe with a column for the conformer string and another 
            column for the occurrence of the conformer (counts)
        angles : list
            Updated feature list with branch separators
        branches : list
            List of branches present in the processed glycan structure
        """
    
        colvar = _readinputfile(self.inputfile, self.angles)
        colvar = colvar.iloc[::round(colvar.shape[0]/self.length), :]

        binary = _readinputfile(self.inputfile, self.angles)
        binary = binary.iloc[::round(binary.shape[0]/self.length), :]
    
        for f in self.angles:
            if len(self.maxima["{}".format(f)]) == 1:

                binary["{}".format(f)] = np.where((-3.5 <= colvar["{}".format(f)]) & (colvar["{}".format(f)] <= 3.5) , self.label["{}".format(f)][0], binary["{}".format(f)])

            elif len(self.maxima["{}".format(f)]) == 2:

                binary["{}".format(f)] = np.where((-3.5 <= colvar["{}".format(f)]) & (colvar["{}".format(f)] < self.maxima["{}".format(f)][0]) , self.label["{}".format(f)][0], binary["{}".format(f)])
                binary["{}".format(f)] = np.where((self.maxima["{}".format(f)][0] <= colvar["{}".format(f)]) & (colvar["{}".format(f)] < self.maxima["{}".format(f)][1]) , self.label["{}".format(f)][1], binary["{}".format(f)])
                binary["{}".format(f)] = np.where((self.maxima["{}".format(f)][1] <= colvar["{}".format(f)]) & (colvar["{}".format(f)] <= 3.5) , self.label["{}".format(f)][2], binary["{}".format(f)])      

            elif len(self.maxima["{}".format(f)]) == 3:

                binary["{}".format(f)] = np.where((-3.5 <= colvar["{}".format(f)]) & (colvar["{}".format(f)] < self.maxima["{}".format(f)][0]) , self.label["{}".format(f)][0], binary["{}".format(f)])
                binary["{}".format(f)] = np.where((self.maxima["{}".format(f)][0] <= colvar["{}".format(f)]) & (colvar["{}".format(f)] < self.maxima["{}".format(f)][1]) , self.label["{}".format(f)][1], binary["{}".format(f)])
                binary["{}".format(f)] = np.where((self.maxima["{}".format(f)][1] <= colvar["{}".format(f)]) & (colvar["{}".format(f)] < self.maxima["{}".format(f)][2]) , self.label["{}".format(f)][2], binary["{}".format(f)])
                binary["{}".format(f)] = np.where((self.maxima["{}".format(f)][2] <= colvar["{}".format(f)]) & (colvar["{}".format(f)] <= 3.5) , self.label["{}".format(f)][3], binary["{}".format(f)])


        binary = pd.DataFrame(binary, columns = binary.columns)
        binary = binary.replace({' t ': ' T ', ' TG': ' tg'}, regex=True)
    
        self.branches = []
    
        for i, sep in enumerate(self.separator):
            col_name = f'sep{i + 1}'
            binary.insert(loc=self.separator_index[i], column=col_name, value=sep)
            self.branches.append(sep)
    
        angles_separators = list(binary.columns)
        count = binary.groupby(angles_separators).size().to_frame(name='Count').reset_index()
        count['Conformer'] = count[count.columns[0:len(angles_separators)]].apply(lambda x: ''.join(x.dropna().astype(str)), axis=1)
        count = count.drop(binary.columns[0:len(angles_separators)], axis=1)
        count = count.sort_values("Count", ascending=False, ignore_index=True)
        
        binary_compressed = binary.copy()
        binary_compressed['Conformer'] = binary_compressed[binary_compressed.columns[0:len(binary_compressed.columns)]].apply(lambda x: ''.join(x.dropna().astype(str)),axis=1)
        binary_compressed = binary_compressed.drop(binary_compressed.columns[0:len(binary_compressed.columns)-1], axis=1)

        if self.weights == None:
            pass
        else:
            w = pd.read_csv(self.weights, names = ["weights"], header = 0)
            binary_w = pd.concat([binary_compressed, w], axis=1)

            indexlist = []
            count_weight = []
            for i in range(count.shape[0]):
                index = binary_w.index[binary_w['Conformer'] == count.loc[i,"Conformer"]].tolist()
                indexlist.append(index)
                count_weight.append(len(index))
                
            sum_weight = []
            for i in range(len(indexlist)):
                fraction = binary_w.loc[indexlist[i]]
                sum_fraction = fraction.loc[:, 'weights'].sum()
                sum_weight.append(sum_fraction)
            sum_weight = pd.DataFrame(sum_weight , columns = ['Count'])  
            weight_count = count.drop(columns=["Count"])
            count = pd.concat([weight_count, sum_weight], axis=1)
            count = count.sort_values("Count", ascending=False, ignore_index=True)
        
        return binary, binary_compressed, count, angles_separators
    
    def _perform_block_averages(self):
        """
        Calculates probability of conformers for 10 separate blocks.

        Function that splits the torsion angle dataset in 10 
        blocks to do the conformer labelling and probability calculation 
        for each block.

        Parameters
        ---------- 
        c: str
            Binary COLVAR input file as a pandas dataframe in which the 
            torsion angle values are replaced by labels
        count: mixed
            Pandas dataframe with the first column for the conformer string 
            and the second column for the occurance of the conformer (counts)
        Returns
        -------
        """
        features = list(self.binary.columns)
        count = self.count.set_index('Conformer')
        step = len(self.binary)/10

        self.blockdata = []

        clear1 = np.arange(0, len(self.binary), step)
        clear2 = np.arange(step, len(self.binary) + step, step)
        clear1 = clear1.astype(int)
        clear2 = clear2.astype(int)

        for i, j, k in zip(range(1, 11), clear1, clear2):
            binary_short = self.binary.iloc[j:k]
            bbinary_short = binary_short.groupby(features).size().reset_index(name='Count')
            bbinary_short['Probability'] = bbinary_short['Count'] / len(binary_short)
            bbinary_short['Conformer'] = bbinary_short[features].apply(lambda x: ''.join(x.dropna().astype(str)), axis=1)
            
            # Keep only the relevant columns
            bbinary_short = bbinary_short[['Conformer', 'Count', 'Probability']]
            
            # Join both dataframes taking 'Conformer' as key column for comparison
            bbinary_short = bbinary_short.set_index('Conformer').join(count, how='outer', lsuffix='_block', rsuffix='_total')
            bbinary_short = bbinary_short.round(6).fillna(0.0).reset_index()
            self.blockdata.append(bbinary_short)
            
        N = 1
        average_block = self.blockdata[0]['Probability'].copy()  # Assuming 'Probability' is the column name based on the earlier discussion
        average2_block = self.blockdata[0]['Probability'] ** 2
        for histn in self.blockdata[1:]:
            N += 1
            average_block += histn['Probability']
            average2_block += histn['Probability'] ** 2
        average_block /= N
        average = pd.DataFrame({'Conformer': self.blockdata[0]["Conformer"], 'Prob': average_block})
        error = pd.DataFrame({'Conformer': self.blockdata[0]["Conformer"], 'Error': average2_block})

        return average, error, N

    def _bootstrap(self, n_iterations):
        # Extract the single column of data
        data = self.binary_compressed.loc[:, "Conformer"]
        
        # Dictionary to store the results for each unique string
        unique_strings = data.unique()
        bootstrap_results = {string: [] for string in unique_strings}
        
        for _ in range(n_iterations):
            # Sample with replacement from the original data
            sample = data.sample(n=len(data), replace=True)
            
            # Count occurrences and calculate probability for each unique string
            counts = sample.value_counts(normalize=True)
            
            # Store probabilities for each string
            for string in unique_strings:
                if string in counts:
                    bootstrap_results[string].append(counts[string])
                else:
                    bootstrap_results[string].append(0)  # Add zero if the string didn't appear in the sample

        # Calculate average probabilities and standard deviations
        averages = {string: np.mean(probabilities) for string, probabilities in bootstrap_results.items()}
        errors = {string: np.std(probabilities) for string, probabilities in bootstrap_results.items()}
        
        # Creating DataFrames for averages and errors
        average = pd.DataFrame(list(averages.items()), columns=['Conformer', 'Prob'])
        error = pd.DataFrame(list(errors.items()), columns=['Conformer', 'Error'])

        return average, error
    
    def _order_conformer(self, average, error, threshold):

        average = average.drop(average[average['Prob'] < threshold/100 ].index)
        average = average.sort_values(by="Prob", ascending=False)
        indexlist = average.index.tolist()
        filtered_conformer = average["Conformer"].filter(indexlist, axis = 0)
        namelist = filtered_conformer.tolist()
        name_list = _include_branch_conformer_string(namelist, self.branches, self.angles_separator)
        name_list = _vertical_conformer_string(name_list)
        error = error.filter(indexlist, axis = 0)

        return name_list, average, error

    def _plot_distribution(self, name_list, colors, dpi, ymax, fontsize, file, average, error):
        
        pos_list = np.arange(len(name_list))

        plt.rcParams['figure.dpi'] = dpi
        ax = plt.axes()
        ax.xaxis.set_major_locator(ticker.FixedLocator((pos_list)))
        ax.xaxis.set_major_formatter(ticker.FixedFormatter((name_list)))   
        bar = plt.bar(pos_list, average.Prob * 100, yerr = error * 100)
        for i in range(len(pos_list)):
            bar[i].set_color(colors[i])
        plt.ylim(0,ymax)
        plt.xticks(fontsize = fontsize)
        plt.yticks(fontsize = fontsize)
        plt.xlabel("Conformer", fontsize = fontsize)
        plt.ylabel('Probability [%]', fontsize = fontsize)
        
        if file is None:
            pass
        else:
            plt.savefig(file, bbox_inches='tight')
        plt.show()

    def distribution(self,
                     threshold=2,
                     ymax=100, 
                     fontsize=15, 
                     colors=["#173c4d","#146b65","#4e9973","#a7c09f","#dfa790","#c76156","#9a2b4c","#600b4a"],
                     dpi=300,
                     file=None,
                     n_iterations=1000):
        """
        Plots a bar graph of the conformer density function.

        Function that plot the probability distribution of the 
        most occuring conformers. The limit defines the cutoff
        probability, where conformers are excluded from the 
        representation if below the cutoff.

        Parameters
        ----------
        threshold: float
            Propability limit up to which conformers are included in the graph 
        ymax: float
            Maximum height of the y-axis (up to 100, as we plot probabilites), default is 100%
        fontsize: int
            Text size
        colors: list of strings
            List of matplotlib color strings used to color the individual populations
        dpi: int
            Resolution of plot
        file: str
            Path+file, where to store the figure of the plot

        Returns
        -------
        """
        if self.weights is None:
            average, error, N = self._perform_block_averages()
            name_list, average, error = self._order_conformer(average, error, threshold) 
            var = (N/(N-1))*( (error.Error / N) - average.Prob*average.Prob ) 
            error = np.sqrt( var / N )
            self._plot_distribution(name_list, colors, dpi, ymax, fontsize, file, average, error)
       
        else:
            average, error = self._bootstrap(n_iterations)
            name_list, average, error = self._order_conformer(average, error, threshold) 
            error = error.Error
            self._plot_distribution(name_list, colors, dpi, ymax, fontsize, file, average, error)

    def validate_fep(self): 

        """
        Function that plots free energy profiles with annotated minima and maxima.

        Evaluation of how good the identification of minima and maxima worked.
        -------
        """
        # implement plotting also of reweighted torsion FEPs

        y0 = 0
        y1 = 40
        binary= np.arange(-0.5236, 0 + 0.001, 0.001)
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
            axs[i].fill_between(binary, y0, y1, alpha = 0.6, color = "paleturquoise", edgecolor = "none")
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

    def cumulative_average(self,
                           simulation_length,
                           ranks = 3,
                           label = ["1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th"],
                           color = ["#173c4d","#146b65","#4e9973","#a7c09f","#dfa790","#c76156","#9a2b4c","#600b4a"],
                           linestyle = ["-", "-", "-", "-", "-", "-", "-", "-"],
                           fontsize = 15,
                           dpi = 300,
                           ymax = 100,
                           file=None):
        if self.weights == None:
            #List top conformers
            top = []
            for r in range(0,ranks):
                top.append(self.count.loc[r, "Conformer"])
            #Count occurrance of top conformers
            occurrences = [[] for _ in range(len(top))]
            indices = [[] for _ in range(len(top))]
            for j, conf in enumerate(top):
                c = 0
                for i in self.binary_compressed.index:
                    if self.binary_compressed.loc[i, "Conformer"] == conf: 
                        c = c + 1
                        occurrences[j].append((c) / (i + 1))
                        indices[j].append(i)
                    else:
                        occurrences[j].append((c) / (i + 1))
                        indices[j].append(i)
                        
                occurrences[j][:] = [x * 100 for x in occurrences[j]]
                indices[j][:] = [x * (simulation_length/self.length) for x in indices[j]]
            #Plot
            plt.figure()
            plt.rcParams['figure.dpi'] = dpi
            plt.xlabel("Time [ns]", fontsize=fontsize)
            plt.ylabel("Probability [%]", fontsize=fontsize)
            plt.tick_params(axis='both', which='major', labelsize=fontsize)
            
            for j,l,c,line in zip(range(0,len(top)), label, color, linestyle):
                plt.plot(indices[j], occurrences[j], label = l, color = c, linestyle = line)
            
            plt.legend(fontsize = fontsize)
            plt.ylim(0,ymax)
            plt.plot()

            if file is None:
                pass
            else:
                plt.savefig(file, bbox_inches='tight')
            plt.show()
        else:
            print("Cumulative average can not be computed for weighted data")

    def moving_average(self,
                       simulation_length,
                       window = 5000,
                       ranks = 3,
                       label = ["1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th"],
                       color = ["#173c4d","#146b65","#4e9973","#a7c09f","#dfa790","#c76156","#9a2b4c","#600b4a"],
                       linestyle = ["-", "-", "-", "-", "-", "-", "-", "-"],
                       fontsize = 15,
                       dpi = 600,
                       ymax = 100,
                       file = None):
        if self.weights == None:
            #List top conformers
            top = []
            for r in range(0,ranks):
                top.append(self.count.loc[r, "Conformer"])
            
            #Count rolling occurrance of top conformers
            occurrences = [[] for _ in range(len(top))]
            
            for j, conf in enumerate(top):
                occurrences[j] = list(self.binary_compressed.loc[:,"Conformer"])
            
                for i in range(0, len(occurrences[j])):
                    if occurrences[j][i] == conf:
                        pass
                    else:
                        occurrences[j][i] = np.nan
                        
                occurrences[j] = pd.DataFrame(occurrences[j])   
                occurrences[j]["rollcount"] = occurrences[j].rolling(window, min_periods=1).count()
                occurrences[j].loc[:, "rollcount"] = [x / window * 100 for x in occurrences[j].loc[:, "rollcount"]]
            
            #Plot
            plt.figure()
            plt.rcParams['figure.dpi'] = dpi
            plt.xlabel("Time [ns]", fontsize=fontsize)
            plt.ylabel("Probability [%]", fontsize=fontsize)
            plt.tick_params(axis='both', which='major', labelsize=fontsize)
            
            for j,l,c,line in zip(range(0,len(top)), label, color, linestyle):
                plt.plot(occurrences[j].index * (simulation_length/self.length), occurrences[j].loc[:, "rollcount"], label = l, color = c, linestyle = line)
            
            plt.legend(fontsize = fontsize)
            plt.ylim(0,ymax)
            plt.show()

            if file is None:
                pass
            else:
                plt.savefig(file, bbox_inches='tight')
            plt.show()
        else:
            print("Moving average can not be computed for weighted data")

    def _plot_pca(self, top, label, color, components_plot, dpi, figsize, fontsize, all, all_color, all_label, marker, conformer, legend, file, biplot, biplot_fontsize,pick, datatopick, colorpick, coefficients, ticks, finalDf=None, pca_df_scaled=None, loadings_r=None):
            #Plot
            plt.rcParams['figure.dpi'] = dpi
            fig = plt.figure(figsize = (figsize))
            ax = fig.add_subplot(1,1,1) 
            ax.set_xlabel('Principal Component {}'.format(components_plot[0]), fontsize = fontsize)
            ax.set_ylabel('Principal Component {}'.format(components_plot[1]), fontsize = fontsize)
            #Ticks
            if ticks == False:
                ax.get_xaxis().set_ticks([])
                ax.get_yaxis().set_ticks([])
            else:
                plt.xticks(fontsize=fontsize)
                plt.yticks(fontsize=fontsize)
                pass
            #plot all
            if all == True:
                if biplot == True:
                    ax.scatter(pca_df_scaled.iloc[:,0], pca_df_scaled.iloc[:,1], c = all_color, alpha = 0.3, marker = marker, s = 5, label = all_label)
                else:
                    ax.scatter(finalDf.iloc[:,components_plot[0]-1], finalDf.iloc[:,components_plot[1]-1], c = all_color, alpha = 0.3, marker = marker, s = 5, label = all_label)
            #plot conformers
            else:
                pass
            if conformer == True:
                if biplot == True:
                    for i, t in enumerate(top):
                        indicesToKeep = pca_df_scaled.loc[:,"Conformer"] == t
                        ax.scatter(pca_df_scaled.loc[indicesToKeep,0], pca_df_scaled.loc[indicesToKeep,1]
                                   , c = color[i]
                                   , s = 30
                                   , label = label[i]
                                   , edgecolors = "black"
                                   , marker = marker
                                   , linewidth = 0.2
                                   , alpha = 1
                                       )
                else:
                    for i, t in enumerate(top):
                        indicesToKeep = finalDf.loc[:,"Conformer"] == t
                        ax.scatter(finalDf.loc[indicesToKeep,components_plot[0]-1], finalDf.loc[indicesToKeep,components_plot[1]-1]
                                   , c = color[i]
                                   , s = 30
                                   , label = label[i]
                                   , edgecolors = "black"
                                   , marker = marker
                                   , linewidth = 0.2
                                   , alpha = 1
                                       )
            else:
                pass
            if file is None:
                pass
            else:
                plt.savefig(file, bbox_inches='tight')
            if legend == True:
                ax.legend(fontsize = fontsize)
            else:
                pass
            if pick == True:
                if biplot == True:
                    for d in datatopick:
                        ax.scatter(pca_df_scaled.iloc[d,0], pca_df_scaled.iloc[d,0], c = colorpick, marker = "x", s = 5)
                else:
                    for d in datatopick:
                        ax.scatter(finalDf.iloc[d,components_plot[0]-1], finalDf.iloc[d,components_plot[1]-1], c = colorpick, marker = "x", s = 5)
            else:
                pass
            if biplot == True:
                #plot coefficient
                distances = []
                for row in range(0,len(loadings_r)):
                    distances.append(((0-loadings_r.iloc[row,components_plot[0]])**2 + (0-loadings_r.iloc[row,components_plot[1]])**2)**0.5)
                loadings_r['vector_norm'] = distances
                loadings_v = loadings_r.sort_values(by=['vector_norm'], ascending = False)
                #
                for i in range(0,coefficients):
                    plt.text(loadings_v.iloc[i,components_plot[0]]+0.025, loadings_v.iloc[i,components_plot[1]]+0.025, loadings_v.iloc[i,0], fontsize=biplot_fontsize, bbox=dict(boxstyle='round', alpha = 0.7, facecolor = "white", edgecolor="gray"))
                    ax.scatter(loadings_v.iloc[i,components_plot[0]], loadings_v.iloc[i,components_plot[1]], s=100, c = "darkgray", alpha = 0.7)
                    plt.arrow(
                        0, 0, # coordinates of arrow base
                        loadings_v.iloc[i,components_plot[0]], # length of the arrow along x
                        loadings_v.iloc[i,components_plot[1]], # length of the arrow along y
                        color='black', 
                        head_width=0.01,
                        linewidth=0.2
                        )
            else:
                pass
            plt.show()

    def pca(self,
            components = 2,
            ranks = 3,
            components_plot = [1,2], #only 2D supported #
            all = True,
            all_color = "gray",
            all_label = "all",
            marker = ".",
            conformer = True,
            label = ["1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th"],
            color = ["#173c4d","#146b65","#4e9973","#a7c09f","#dfa790","#c76156","#9a2b4c","#600b4a"],
            fontsize = 7,
            dpi = 600,
            figsize = [5,5],
            ticks = True,
            legend = True,
            pick = False,
            datatopick = [0,0],
            colorpick = "darkred",
            biplot = False,
            biplot_fontsize = 7,
            coefficients = 3,
            file = None):
        
        colvar_pca = pd.DataFrame({"index": self.colvar.index})
        for a in self.angles:
            colvar_pca.loc[:,"sin_{}".format(a)] = np.sin(self.colvar.loc[:,a])
            colvar_pca.loc[:,"cos_{}".format(a)] = np.cos(self.colvar.loc[:,a])
        colvar_pca = colvar_pca.drop(columns = ["index"])
        
        #List top conformers
        top = []
        for r in range(0,ranks):
            top.append(self.count.loc[r, "Conformer"])
        if self.weights == None: 
            # Choose number of principle components
            pca = PCA(n_components=components)
            # Compute PCA and transform into dataframe with target addition
            principalComponents = pca.fit_transform(colvar_pca)
        else:
            w = pd.read_csv(self.weights, names = ["weights"], header = 0)
            n = len(colvar_pca.columns)
            weights_pca = pd.concat([w] * (n+1), axis=1, ignore_index=True)
            weights_pca = pd.DataFrame(data = weights_pca, columns = np.linspace(0,n-1,n))
    
            # Choose number of principle components
            pca = EMPCA(n_components=components)
            # Compute PCA and transform into dataframe with target addition
            principalComponents = pca.fit_transform(colvar_pca, weights = weights_pca)

        principalDf = pd.DataFrame(data = principalComponents)
        finalDf = pd.concat([principalDf, self.binary_compressed[['Conformer']]], axis = 1)

        if biplot is True:
            # Principal components correlation coefficients, eigenvectors of CoVar Matrix
            loadings = pca.components_
            # Feature names before PCA
            feature_names = colvar_pca.columns
            # PC names
            pc_list = [f'PC{i}' for i in list(range(1, components + 1))]
            # Match PC names to loadings
            pc_loadings = dict(zip(pc_list, loadings))
            # Matrix of corr coefs between feature names and PCs
            loadings_df = pd.DataFrame.from_dict(pc_loadings)
            loadings_df['feature_names'] = feature_names
            loadings_df = loadings_df.set_index('feature_names')

            #Average vectors of sin and cos to make the interpretation more comprehensive  
            comp_average = [[] for _ in range(0,components)]
            loadings_r = pd.DataFrame({"feature_names": self.colvar.columns})
            for j, pc in enumerate(pc_list):
                for i in range(0,len(loadings_df),2):
                    comp_average[j].append((loadings_df.iloc[i,j] + loadings_df.iloc[i+1,j])/2)
                loadings_r[pc] = comp_average[j]

            # Scale PCS into a DataFrame
            pca_df_scaled = finalDf.copy()
            scaler_df = finalDf.iloc[:,components_plot[0]-1]
            scaler_df = pd.DataFrame(data = scaler_df)
            scaler_df = scaler_df.T.reset_index(drop=True).T
            scaler_df[len(scaler_df.columns)] = finalDf.iloc[:,components_plot[1]-1]
            pca_df_scaled = scaler_df.copy()
            pca_df_scaled["Conformer"] = finalDf.loc[:,"Conformer"]
            scaler = 1 / (scaler_df.max() - scaler_df.min())
            for index in scaler.index:
                pca_df_scaled[index] *= scaler[index]
       
            self._plot_pca(top, label, color, components_plot, dpi, figsize, fontsize, all, all_color, all_label, marker, conformer, legend, file, biplot, biplot_fontsize, pick, datatopick, colorpick, coefficients, ticks, pca_df_scaled=pca_df_scaled, loadings_r=loadings_r)
        else:
            self._plot_pca(top, label, color, components_plot, dpi, figsize, fontsize, all, all_color, all_label, marker, conformer, legend, file, biplot, biplot_fontsize, pick, datatopick, colorpick, coefficients, ticks, finalDf=finalDf)

    def pca_fep(self,
                components = 2,
                fontsize = 10,
                dpi = 600,
                components_plot = [1,2],
                figsize = [5,5],
                legend = True,
                ticks = True,
                bins = 50,
                cmap = "jet",
                colorrange = 50,
                alpha = 0.8,
                file = None):

        colvar_pca = pd.DataFrame({"index": self.colvar.index})
        for a in self.angles:
            colvar_pca.loc[:,"sin_{}".format(a)] = np.sin(self.colvar.loc[:,a])
            colvar_pca.loc[:,"cos_{}".format(a)] = np.cos(self.colvar.loc[:,a])
        colvar_pca = colvar_pca.drop(columns = ["index"])

        pca = PCA(n_components=components)
        # Compute PCA and transform into dataframe with target addition
        principalComponents = pca.fit_transform(colvar_pca)
        principalDf = pd.DataFrame(data = principalComponents)
        finalDf = pd.concat([principalDf, self.binary_compressed[['Conformer']]], axis = 1)
        
        H, xedges, yedges = np.histogram2d(finalDf.iloc[:,components_plot[0]-1], finalDf.iloc[:,components_plot[1]-1], bins=bins)

        np.seterr(divide='ignore')
        RT = 2.479
        F = - RT*np.log(H)
        
        F_min = F.min()
        if F_min > 0:
            F = F - F_min
        elif F_min == 0:
            print("Data is already shifted to 0 for its minimum value")
        else:
            F = F + abs(F_min)
        
        plt.rcParams['figure.dpi'] = dpi
        plt.figure(figsize = (figsize))
        
        ax = plt.gca()
        if ticks == False:
            ax.axes.xaxis.set_visible(False)
            ax.axes.yaxis.set_visible(False)
        else:
            pass
            
        ax.set_xlabel('Principal Component {}'.format(components_plot[0]), fontsize = fontsize)
        ax.set_ylabel('Principal Component {}'.format(components_plot[1]), fontsize = fontsize)
        
        cp = plt.contourf(yedges[1:], xedges[1:], F.T, colorrange, cmap=cmap, alpha = alpha)
        
        if legend == True:
            cbar = plt.colorbar(cp)
            cbar.set_label("kJ/mol", fontsize = fontsize)
        else:
            pass
            
        if file is None:
            pass
        else:
            plt.savefig(file, bbox_inches='tight')
        plt.show()

class Glycompare:
    # Class variables
    ## maybe add directory names ?
    def __init__(self, inputfile, inputdir, common_angles, common_branches, common_angles_separators, count, binary_compressed, outputdir = "./", n_blocks = 10, weights = None): # check later which parameters are really necessary?!
    
        # Instance variables
        """
        Initialize object's attributes, aka setting variables.

        Parameters
        ----------
        inputdir: str
            List of directory names to read in, storing the Cluster_conformer.dat files
            Files for each simulation should have a separate directory
        outputdir: str
            Path to output directory, where created files are stored
        common_angles: list of str
            List of common torsion angle names 

        # complete this!!!
        
        Returns
        -------
        """
        self.inputfile = inputfile
        self.inputdir = inputdir
        self.outputdir = outputdir 
        self.common_angles = common_angles
        self.common_branches = common_branches
        self.common_angles_separators = common_angles_separators
        self.n_blocks = n_blocks
        self.count = count
        self.binary_compressed = binary_compressed
        self.weights = weights

    def _join_data(self): #limited to 4 datasets, but comparing more does anyway makes not much sense
        """
        Function that reads in the Cluster_conformer data for n simulations and joins dataframes.
        
        Output:
            csv-files for each of the blocks (from block averaging) cointaining the joint dataframes
        """
        
        header = ["Index", "Conformer", "Count_partial", "Probability", "Count_full"]
        blocks = [[] for _ in range(len(self.inputdir))]
        
        # Load partially data sorted for n blocks
        for b in range(1, self.n_blocks + 1):
            for i, j in enumerate(self.inputdir):
                blocks[i] = pd.read_csv("{}/Cluster_conformer{}.dat".format(j,b), delim_whitespace=True, names = header, header = 0, dtype = str)
                blocks[i] = blocks[i].drop(columns = ['Count_full', 'Count_partial', 'Index'])
                blocks[i] = blocks[i].set_index('Conformer') 
            if len(blocks) == 2:
                    jointblocks = pd.concat([blocks[0], blocks[1]], axis=1)
            if len(blocks) == 3:
                    jointblocks = pd.concat([blocks[0], blocks[1], blocks[2]], axis=1)
            if len(blocks) == 4:
                    jointblocks = pd.concat([blocks[0], blocks[1], blocks[2], blocks[3]], axis=1)        
            jointblocks = jointblocks.fillna(0.0).reset_index()
            jointblocks.to_csv("{}/Joint_cluster_conformer{}.dat".format(self.outputdir,b), sep= " ", header = False)

    def plot(self, # so far only able to plot two comparative plots
             threshold = 2,
             ymax = 100, 
             fontsize = 15, 
             colors = ["#173c4d","#146b65","#4e9973","#a7c09f","#dfa790","#c76156","#9a2b4c","#600b4a","#6C7D85","#719F91","#96BBAB","#C9D8C6","#EAC9BC","#DBA099","#C38196","#996C92"],
             dpi = 300,
             label = ["Glycan 1", "Glycan 2", "Glycan 3", "Glycan 4"], 
             file = None):
        """ Function that reads in the Joint cluster conformer data, calculates the block averages, standard deviations and plots everything in a bar graph.
        """
        # Join dataframes
        self._join_data()
        
        # Initialize variables
        prob = [[] for _ in range(len(self.inputdir))]
        squared_prob = [[] for _ in range(len(self.inputdir))]
        average_prob = [[] for _ in range(len(self.inputdir))]
        errors = [[] for _ in range(len(self.inputdir))]
        
        # Load partially sorted data and calculate averages
        for d in range(len(self.inputdir)):    
            for i in range(1, self.n_blocks):
                blocks = pd.read_csv(f"{self.outputdir}/Joint_cluster_conformer{i}.dat", delim_whitespace=True, usecols=[d+2])
                if len(prob[d]) == 0:
                    prob[d] = blocks.iloc[:,0]
                    squared_prob[d] = blocks.iloc[:,0]**2
                else:
                    prob[d] += blocks.iloc[:,0]
                    squared_prob[d] += blocks.iloc[:,0]**2
            average_prob[d] = (prob[d] / self.n_blocks).to_frame()

        if len(prob) == 2:
            total_average = pd.concat([average_prob[0], average_prob[1]], axis=1)        
        elif len(prob) == 3:
            total_average = pd.concat([average_prob[0], average_prob[1], average_prob[2]], axis=1) 
        elif len(prob) == 4:
            total_average = pd.concat([average_prob[0], average_prob[1], average_prob[2], average_prob[3]], axis=1)
        else:
            print("GlyCOMPARE can only handle up to four datasets!")
        
        total_average = total_average.drop(total_average[total_average.iloc[:,0] < threshold/100].index, inplace=False)
        total_average = total_average.sort_values(by = "0.0",ascending=False)
        indexlist = total_average.index.tolist()

        for d in range(len(self.inputdir)):
            average_prob[d] = average_prob[d].filter(indexlist, axis = 0)
            squared_prob[d] = squared_prob[d].loc[indexlist].to_frame()
            errors[d] = abs((self.n_blocks / (self.n_blocks - 1)) * (squared_prob[d] / self.n_blocks - average_prob[d]**2))
            errors[d] = np.sqrt(errors[d] / self.n_blocks)

        # Create namelist
        dummy = pd.read_csv(f"{self.outputdir}/Joint_cluster_conformer1.dat", delim_whitespace=True, index_col = [0], header = None)
        dummy = dummy.filter(indexlist, axis = 0)
        namelist = dummy.iloc[:,0].tolist()

        # Make conformer string look nice
        pos_list = np.arange(len(namelist))
        namelist = _include_branch_conformer_string(namelist, self.common_branches, self.common_angles_separators)
        namelist_v = _vertical_conformer_string(namelist)
        
        plt.rcParams['figure.dpi'] = dpi
        ax = plt.axes()
        ax.xaxis.set_major_locator(ticker.FixedLocator((pos_list)))
        ax.xaxis.set_major_formatter(ticker.FixedFormatter((namelist_v)))   
        plt.tick_params(axis='both', which='major', labelsize=fontsize)
        plt.xlabel("Conformer", fontsize=fontsize)
        plt.ylabel('Probability [%]', fontsize=fontsize)
        if len(prob) == 2:
            width = 0.42
            for i,c in zip(range(0,len(pos_list)),colors):
                plt.bar(pos_list[i] - width/2, average_prob[0].iloc[i,0] * 100, width = width, color = c, yerr = errors[0].iloc[i,0] * 100, edgecolor = "black")
                plt.bar(pos_list[i] + width/2, average_prob[1].iloc[i,0] * 100, width = width, color = c, yerr = errors[1].iloc[i,0] * 100, edgecolor = "black", hatch = "//")
            l1_patch = mpatches.Patch(edgecolor= "black", facecolor="white", label=label[0])
            l2_patch = mpatches.Patch(edgecolor= "black", facecolor="white", label=label[1], hatch = "//")
            plt.legend(handles=[l1_patch, l2_patch], fontsize=fontsize)
        elif len(prob) == 3:
            width = 0.3
            for i,c in zip(range(0,len(pos_list)),colors):
                plt.bar(pos_list[i] - width, average_prob[0].iloc[i,0] * 100, width = width, color = c, yerr = errors[0].iloc[i,0] * 100, edgecolor = "black")
                plt.bar(pos_list[i], average_prob[1].iloc[i,0] * 100, width = width, color = c, yerr = errors[1].iloc[i,0] * 100, edgecolor = "black", hatch = "//")
                plt.bar(pos_list[i] + width, average_prob[2].iloc[i,0] * 100, width = width, color = c, yerr = errors[2].iloc[i,0] * 100, edgecolor = "black", hatch = "o")
            l1_patch = mpatches.Patch(edgecolor= "black", facecolor="white", label=label[0])
            l2_patch = mpatches.Patch(edgecolor= "black", facecolor="white", label=label[1], hatch = "//")
            l3_patch = mpatches.Patch(edgecolor= "black", facecolor="white", label=label[2], hatch = "o")
            plt.legend(handles=[l1_patch, l2_patch, l3_patch], fontsize=fontsize)
        elif len(prob) == 4:   
            width = 0.2
            for i,c in zip(range(0,len(pos_list)),colors):
                plt.bar(pos_list[i] - width*1.5, average_prob[0].iloc[i,0] * 100, width = width, color = c, yerr = errors[0].iloc[i,0] * 100, edgecolor = "black")
                plt.bar(pos_list[i] - width/2, average_prob[1].iloc[i,0] * 100, width = width, color = c, yerr = errors[1].iloc[i,0] * 100, edgecolor = "black", hatch = "//")
                plt.bar(pos_list[i] + width/2, average_prob[2].iloc[i,0] * 100, width = width, color = c, yerr = errors[2].iloc[i,0] * 100, edgecolor = "black", hatch = "o")
                plt.bar(pos_list[i] + width*1.5, average_prob[3].iloc[i,0] * 100, width = width, color = c, yerr = errors[3].iloc[i,0] * 100, edgecolor = "black", hatch = "x")
            l1_patch = mpatches.Patch(edgecolor= "black", facecolor="white", label=label[0])
            l2_patch = mpatches.Patch(edgecolor= "black", facecolor="white", label=label[1], hatch = "//")
            l3_patch = mpatches.Patch(edgecolor= "black", facecolor="white", label=label[2], hatch = "o")
            l4_patch = mpatches.Patch(edgecolor= "black", facecolor="white", label=label[3], hatch = "x")
            plt.legend(handles=[l1_patch, l2_patch, l3_patch, l4_patch], fontsize=fontsize)
        plt.ylim(0,ymax) 
        if file is None:
            pass
        else:
            plt.savefig(file, bbox_inches='tight')
        plt.show()
    
    def pca(self, #do we want to return data that is useful for PCA analysis? Also in Glyconformer PCA
            length = None,
            components = 2,
            ranks = 3,
            all = True,
            conformer = True,
            components_plot = [1,2], #only 2D supported 
            #
            dataset_label = ["Glycan 1", "Glycan 2", "Glycan 3", "Glycan 4"], 
            label = ["1st", "2nd", "3rd", "4th", "5th", "6th", "7th", "8th"],
            color = ["#173c4d","#146b65","#4e9973","#a7c09f","#dfa790","#c76156","#9a2b4c","#600b4a"],
            color_basic = ["gray","gray","gray","gray"],
            marker = [".","x","v","o"],
            fontsize = 7,
            dpi = 600,
            figsize = [5,5],
            ticks = True,
            legend = True,
            pick = False,
            datatopick = [0,0],
            colorpick = "darkred",
            biplot = False,
            coefficients = 3,
            file = None):

        colvars = [[] for _ in range(len(self.inputdir))]
        colvars_pca = [[] for _ in range(len(self.inputdir))]
        lengths = [[] for _ in range(len(self.inputdir))] 
        top = [[] for _ in range(len(self.inputdir))]
        finalDf = [[] for _ in range(len(self.inputdir))]
        w = [[] for _ in range(len(self.inputdir))]
        weights_pca = [[] for _ in range (len(self.inputdir))]
        
        for i,f in enumerate(self.inputfile):
            
            if length is None:
                colvars[i] = _readinputfile(f,self.common_angles)
                lengths[i] = len(colvars[i]) 
            else:
                lengths[i] = length 
                colvars[i] = _readinputfile(f,self.common_angles)
                colvars[i] = colvars[i].iloc[::round(colvars[i].shape[0]/lengths[i]), :]
            
            colvars_pca[i] = pd.DataFrame({"index": colvars[i].index})
            for a in self.common_angles:
                colvars_pca[i].loc[:,"sin_{}".format(a)] = np.sin(colvars[i].loc[:,a])
                colvars_pca[i].loc[:,"cos_{}".format(a)] = np.cos(colvars[i].loc[:,a])
            colvars_pca[i] = colvars_pca[i].drop(columns = ["index"])

        if len(self.inputdir) == 2:
            colvar_pca = pd.concat([colvars_pca[0], colvars_pca[1]], axis=0)
            conformer_pca = pd.concat([self.binary_compressed[0], self.binary_compressed[1]], axis=0)
            l_start = [0, lengths[0]+1]
            l_end = [lengths[0]+1, lengths[0]+lengths[1]+2]
        elif len(self.inputdir) == 3:
            colvar_pca = pd.concat([colvars_pca[0], colvars_pca[1], colvars_pca[2]], axis=0)
            conformer_pca = pd.concat([self.binary_compressed[0], self.binary_compressed[1], self.binary_compressed[2]], axis=0)
            l_start = [0, lengths[0]+1, lengths[0]+lengths[1]+2]
            l_end = [lengths[0]+1, lengths[0]+lengths[1]+2, lengths[0]+lengths[1]+lengths[2]+3]
        elif len(self.inputdir) == 4:
            colvar_pca = pd.concat([colvars_pca[0], colvars_pca[1], colvars_pca[2], colvars_pca[3]], axis=0)
            conformer_pca = pd.concat([self.binary_compressed[0], self.binary_compressed[1], self.binary_compressed[2], self.binary_compressed[3]], axis=0)
            l_start = [0, lengths[0]+1, lengths[0]+lengths[1]+2, lengths[0]+lengths[1]+lengths[2]+3]
            l_end = [lengths[0]+1, lengths[0]+lengths[1]+2, lengths[0]+lengths[1]+lengths[2]+3, lengths[0]+lengths[1]+lengths[2]+lengths[3]+4]

        if self.weights == None: 
            # Choose number of principle components
            pca = PCA(n_components=components)
            # Compute PCA and transform into dataframe with target addition
            principalComponents = pca.fit_transform(colvar_pca)
        else:
            for d in range(len(self.inputdir)):
                w[d] = pd.read_csv(self.weights[d], names = ["weights"], header = 0)
                n = len(colvar_pca.columns)
                weights_pca[d] = pd.concat([w[d]] * (n+1), axis=1, ignore_index=True)
                weights_pca[d] = pd.DataFrame(data = weights_pca[d], columns = np.linspace(0,n-1,n))
            if len(self.inputdir) == 2:
                weights_pca = pd.concat([weights_pca[0], weights_pca[1]], axis=0)
            elif len(self.inputdir) == 3:
                weights_pca = pd.concat([weights_pca[0], weights_pca[1], weights_pca[2]], axis=0)
            elif len(self.inputdir) == 4:
                weights_pca = pd.concat([weights_pca[0], weights_pca[1], weights_pca[2], weights_pca[3]], axis=0)
            # Choose number of principle components
            pca = WPCA(n_components=components)
            # Compute PCA and transform into dataframe with target addition
            principalComponents = pca.fit_transform(colvar_pca, weights = weights_pca)
        
        principalDf = pd.DataFrame(data = principalComponents)
        conformer_pca = conformer_pca.reset_index()
        principalDf = pd.concat([principalDf, conformer_pca[['Conformer']]], axis = 1)

        if biplot is True:
            # Principal components correlation coefficients, eigenvectors of CoVar Matrix
            loadings = pca.components_
            # Feature names before PCA
            feature_names = colvar_pca.columns
            # PC names
            pc_list = [f'PC{i}' for i in list(range(1, components + 1))]
            # Match PC names to loadings
            pc_loadings = dict(zip(pc_list, loadings))
            # Matrix of corr coefs between feature names and PCs
            loadings_df = pd.DataFrame.from_dict(pc_loadings)
            loadings_df['feature_names'] = feature_names
            loadings_df = loadings_df.set_index('feature_names')
            
            #Average vectors of sin and cos to make the interpretation more comprehensive  
            comp_average = [[] for _ in range(0,components)]
            loadings_r = pd.DataFrame({"feature_names": colvars[0].columns})
            for j, pc in enumerate(pc_list):
                for i in range(0,len(loadings_df),2):
                    comp_average[j].append((loadings_df.iloc[i,j] + loadings_df.iloc[i+1,j])/2)
                loadings_r[pc] = comp_average[j]
            
            
            #Plot
            plt.rcParams['figure.dpi'] = dpi
            fig = plt.figure(figsize = (figsize))
            ax = fig.add_subplot(1,1,1)
            ax.set_xlabel('Principal Component {}'.format(components_plot[0]), fontsize = fontsize)
            ax.set_ylabel('Principal Component {}'.format(components_plot[1]), fontsize = fontsize)

            if ticks == False:
                ax.get_xaxis().set_ticks([])
                ax.get_yaxis().set_ticks([])
            else:
                pass

            for d, l_s, l_e in zip(range(len(self.inputdir)),l_start,l_end):
                finalDf[d] = principalDf.iloc[l_s:l_e,:]
                
                for r in range(0,ranks):
                    top[d].append(self.count[d].loc[r, "Conformer"])
                
                # Scale PCS into a DataFrame
                pca_df_scaled = finalDf[d].copy()
                scaler_df = finalDf[d].iloc[:,components_plot[0]-1]
                scaler_df = pd.DataFrame(data = scaler_df)
                scaler_df = scaler_df.T.reset_index(drop=True).T
                scaler_df[len(scaler_df.columns)] = finalDf[d].iloc[:,components_plot[1]-1]
                pca_df_scaled = scaler_df.copy()
                pca_df_scaled["Conformer"] = finalDf[d].loc[:,"Conformer"]
                scaler = 1 / (scaler_df.max() - scaler_df.min())
                for index in scaler.index:
                    pca_df_scaled[index] *= scaler[index]

                #plot all
                if all == True: 
                    ax.scatter(pca_df_scaled.iloc[:,0], pca_df_scaled.iloc[:,1], c = color_basic[d], alpha = 0.5, marker = marker[d], s = 5, label = dataset_label[d])
                else:
                    pass
                if conformer == True:
                    for i, t in enumerate(top[d]):
                        indicesToKeep = pca_df_scaled.loc[:,"Conformer"] == t
                        ax.scatter(pca_df_scaled.loc[indicesToKeep,0], pca_df_scaled.loc[indicesToKeep,1]
                                   , c = color[i]
                                   , s = 30
                                   , label = label[i]
                                   , edgecolors = "black"
                                   , marker = marker[d]
                                   , linewidth = 0.2
                                   , alpha = 1
                                       )
                else:
                    pass
            #plot coefficient
            distances = []
            for row in range(0,len(loadings_r)):
                distances.append(((0-loadings_r.iloc[row,components_plot[0]])**2 + (0-loadings_r.iloc[row,components_plot[1]])**2)**0.5)
            loadings_r['vector_norm'] = distances
            loadings_v = loadings_r.sort_values(by=['vector_norm'], ascending = False)
            #
            for i in range(0,coefficients):
                plt.text(loadings_v.iloc[i,components_plot[0]]+0.025, loadings_v.iloc[i,components_plot[1]]+0.025, loadings_v.iloc[i,0], fontsize=5, bbox=dict(boxstyle='round', alpha = 0.7, facecolor = "white", edgecolor="gray"))
                ax.scatter(loadings_v.iloc[i,components_plot[0]], loadings_v.iloc[i,components_plot[1]], s=100, c = "darkgray", alpha = 0.7)
                plt.arrow(
                    0, 0, # coordinates of arrow base
                    loadings_v.iloc[i,components_plot[0]], # length of the arrow along x
                    loadings_v.iloc[i,components_plot[1]], # length of the arrow along y
                    color='black', 
                    head_width=0.01,
                    linewidth=0.2
                    )
            
            #Labels
            if legend == True:
                ax.legend(fontsize = fontsize)
            else:
                pass 
        
        else:
            #Plot
            plt.rcParams['figure.dpi'] = dpi
            fig = plt.figure(figsize = (figsize))
            ax = fig.add_subplot(1,1,1)
            ax.set_xlabel('Principal Component {}'.format(components_plot[0]), fontsize = fontsize)
            ax.set_ylabel('Principal Component {}'.format(components_plot[1]), fontsize = fontsize)
            if ticks == False:
                ax.get_xaxis().set_ticks([])
                ax.get_yaxis().set_ticks([])
            else:
                pass
            
            for d, l_s, l_e in zip(range(len(self.inputdir)),l_start,l_end):
                finalDf[d] = principalDf.iloc[l_s:l_e,:]
                for r in range(0,ranks):
                    top[d].append(self.count[d].loc[r, "Conformer"])
                
                if all == True:
                    ax.scatter(finalDf[d].iloc[:, components_plot[0]-1], finalDf[d].iloc[:,components_plot[1]-1], c = color_basic[d], alpha = 0.5, marker = marker[d], s = 5, label = dataset_label[d])
                else:
                    pass
                if conformer == True:
                    for i, t in enumerate(top[d]):
                        indicesToKeep = finalDf[d].loc[:,"Conformer"] == t
                        ax.scatter(finalDf[d].loc[indicesToKeep,components_plot[0]-1], finalDf[d].loc[indicesToKeep,components_plot[1]-1]
                                   , c = color[i]
                                   , s = 30
                                   , label = label[i]
                                   , edgecolors = "black"
                                   , marker = marker[d]
                                   , linewidth = 0.2
                                   , alpha = 1
                                       )
                else:
                    pass

            if legend == True:
                ax.legend(fontsize = fontsize)
            else:
                pass
        #
        if pick == True:
            for d in range(len(self.inputdir)):
                for p in datatopick:
                    ax.scatter(finalDf[d].iloc[p,components_plot[0]-1], finalDf[d].iloc[p,components_plot[1]-1], c = colorpick, marker = marker[d], s = 5)
        else:
            pass
            
        if file is None:
            pass
        else:
            plt.savefig(file, bbox_inches='tight')
        plt.show()

    def pca_fep(self,
                length = None,
                components = 2,
                fontsize = 10,
                dpi = 600,
                components_plot = [1,2],
                figsize = [5,5],
                legend = True,
                axislim = False,
                ticks = True,
                bins = 50,
                cmap = "jet",
                colorrange = 50,
                vmax = 25,
                alpha = 0.8,
                file = None):

        colvars = [[] for _ in range(len(self.inputdir))]
        colvars_pca = [[] for _ in range(len(self.inputdir))]
        lengths = [[] for _ in range(len(self.inputdir))]
        finalDf = [[] for _ in range(len(self.inputdir))]
        w = [[] for _ in range(len(self.inputdir))]
        weights_pca = [[] for _ in range (len(self.inputdir))]
        
        for i,f in enumerate(self.inputfile):
            
            if length is None:
                colvars[i] = _readinputfile(f,self.common_angles)
                lengths[i] = len(colvars[i]) 
            else:
                lengths[i] = length 
                colvars[i] = _readinputfile(f,self.common_angles)
                colvars[i] = colvars[i].iloc[::round(colvars[i].shape[0]/lengths[i]), :]
            
            colvars_pca[i] = pd.DataFrame({"index": colvars[i].index})
            for a in self.common_angles:
                colvars_pca[i].loc[:,"sin_{}".format(a)] = np.sin(colvars[i].loc[:,a])
                colvars_pca[i].loc[:,"cos_{}".format(a)] = np.cos(colvars[i].loc[:,a])
            colvars_pca[i] = colvars_pca[i].drop(columns = ["index"])

        if len(self.inputdir) == 2:
            colvar_pca = pd.concat([colvars_pca[0], colvars_pca[1]], axis=0)
            conformer_pca = pd.concat([self.binary_compressed[0], self.binary_compressed[1]], axis=0)
            l_start = [0, lengths[0]+1]
            l_end = [lengths[0]+1, lengths[0]+lengths[1]+2]
        elif len(self.inputdir) == 3:
            colvar_pca = pd.concat([colvars_pca[0], colvars_pca[1], colvars_pca[2]], axis=0)
            conformer_pca = pd.concat([self.binary_compressed[0], self.binary_compressed[1], self.binary_compressed[2]], axis=0)
            l_start = [0, lengths[0]+1, lengths[0]+lengths[1]+2]
            l_end = [lengths[0]+1, lengths[0]+lengths[1]+2, lengths[0]+lengths[1]+lengths[2]+3]
        elif len(self.inputdir) == 4:
            colvar_pca = pd.concat([colvars_pca[0], colvars_pca[1], colvars_pca[2], colvars_pca[3]], axis=0)
            conformer_pca = pd.concat([self.binary_compressed[0], self.binary_compressed[1], self.binary_compressed[2], self.binary_compressed[3]], axis=0)
            l_start = [0, lengths[0]+1, lengths[0]+lengths[1]+2, lengths[0]+lengths[1]+lengths[2]+3]
            l_end = [lengths[0]+1, lengths[0]+lengths[1]+2, lengths[0]+lengths[1]+lengths[2]+3, lengths[0]+lengths[1]+lengths[2]+lengths[3]+4]    

        if self.weights == None: 
            # Choose number of principle components
            pca = PCA(n_components=components)
            # Compute PCA and transform into dataframe with target addition
            principalComponents = pca.fit_transform(colvar_pca)
        else:
            for d in range(len(self.inputdir)):
                w[d] = pd.read_csv(self.weights[d], names = ["weights"], header = 0)
                n = len(colvar_pca.columns)
                weights_pca[d] = pd.concat([w[d]] * (n+1), axis=1, ignore_index=True)
                weights_pca[d] = pd.DataFrame(data = weights_pca[d], columns = np.linspace(0,n-1,n))
            if len(self.inputdir) == 2:
                weights_pca = pd.concat([weights_pca[0], weights_pca[1]], axis=0)
            elif len(self.inputdir) == 3:
                weights_pca = pd.concat([weights_pca[0], weights_pca[1], weights_pca[2]], axis=0)
            elif len(self.inputdir) == 4:
                weights_pca = pd.concat([weights_pca[0], weights_pca[1], weights_pca[2], weights_pca[3]], axis=0)
            # Choose number of principle components
            pca = WPCA(n_components=components)
            # Compute PCA and transform into dataframe with target addition
            principalComponents = pca.fit_transform(colvar_pca, weights = weights_pca)
        # Choose number of principle components

        principalDf = pd.DataFrame(data = principalComponents)
        conformer_pca = conformer_pca.reset_index()
        principalDf = pd.concat([principalDf, conformer_pca[['Conformer']]], axis = 1)

        for d, l_s, l_e in zip(range(len(self.inputdir)),l_start,l_end):
            finalDf[d] = principalDf.iloc[l_s:l_e,:]
            H, xedges, yedges = np.histogram2d(finalDf[d].iloc[:,components_plot[0]-1], finalDf[d].iloc[:,components_plot[1]-1], bins=bins)
    
            np.seterr(divide='ignore')
            RT = 2.479
            F = - RT*np.log(H)
            
            F_min = F.min()
            if F_min > 0:
                F = F - F_min
            elif F_min == 0:
                print("Data is already shifted to 0 for its minimum value")
            else:
                F = F + abs(F_min)
            
            plt.rcParams['figure.dpi'] = dpi
            plt.figure(figsize = (figsize))
            
            ax = plt.gca()
            if ticks == False:
                ax.axes.xaxis.set_visible(False)
                ax.axes.yaxis.set_visible(False)
            else:
                pass
                
            ax.set_xlabel('Principal Component {}'.format(components_plot[0]), fontsize = fontsize)
            ax.set_ylabel('Principal Component {}'.format(components_plot[1]), fontsize = fontsize)
            
            if axislim == False:
                pass
            else:
                ax.set_xlim(axislim[0], axislim[1])
                ax.set_ylim(axislim[2], axislim[3])

            cp = plt.contourf(yedges[1:], xedges[1:], F.T, colorrange, vmax = vmax, cmap=cmap, alpha = alpha)
            
            if legend == True:
                cbar = plt.colorbar(cp)
                cbar.set_label("kJ/mol", fontsize = fontsize)
            else:
                pass
                
            if file is None:
                pass
            else:
                plt.savefig(file, bbox_inches='tight')
            plt.show()