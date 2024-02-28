__author__ = "Isabell Louise Grothaus"
__license__ = "GNU GENERAL PUBLIC LICENSE, Version 3"
__version__ = "1.0.0"
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
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


class glyconformer:
    # Class variables
    ## maybe add directory names ?
    def __init__(self, inputfile, outputdir = "./", length = None, glycantype = None , angles = None, omega_angles = None, separator_index = None, separator = None, fepdir = None, order_max = None, order_min = None):
        # Instance variables
        """
        Initialize object's attributes, aka setting variables.

        Parameters
        ----------
        inputfile: str
            Name of the file to read in, including file extension 
        outputdir: str
            Path to output directory, where created files are stored
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
        
        self.inputfile = inputfile
        self.outputdir = outputdir 
        
        #attribute reading mainly depends on if a glycantype is specified and whether
        #information are read from the LIBRARY_GLYCANS or from user input variables. 
       
        if glycantype is None:
            self.angles = angles
            self.omega_angles = omega_angles
            self.separator_index = separator_index
            self.separator = separator
            self.fepdir = fepdir
            self.order_max = order_max
            self.order_min = order_min
            self.maxima,self.minima = self._find_min_max()    
        else:
            self.glycantype = glycantype
            self.angles = self._readfeature("LIBRARY_GLYCANS.{}".format(self.glycantype),"angles.dat")
            self.omega_angles = self._readfeature("LIBRARY_GLYCANS.{}".format(self.glycantype),"omega_angles.dat")
            self.separator_index, self.separator = self._readseparator("LIBRARY_GLYCANS.{}".format(self.glycantype),"separator.dat")
            self.fepdir = "LIBRARY_GLYCANS/{}".format(self.glycantype)
            self.minima = self._readdict("LIBRARY_GLYCANS.{}".format(self.glycantype),"minima.dat")
            self.maxima = self._readdict("LIBRARY_GLYCANS.{}".format(self.glycantype),"maxima.dat")
            
        self.label = self._label_min()
        
        #read inputfile
        
        if length is None:
            self.colvar = self._readinputfile()
            self.length = len(self.colvar) 
        else:
            self.length = length 
            colvar = self._readinputfile()
            self.colvar = colvar.iloc[::round(colvar.shape[0]/self.length), :]
            
    def _readinputfile(self):
        """
        Function reading a dataframe from file, using the self.inputfile.
        """

        try:
            colvar = plumed.read_as_pandas(self.inputfile)
        except:
            colvar = pd.read_csv(self.inputfile, delim_whitespace=True)
        finally:
            colvar = colvar[self.angles]
            
        return colvar
    
    def _readfeature(self, path, file):
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
        with importlib.resources.open_text(path,file) as f:
            reader = csv.reader(f, delimiter=' ')
            for row in reader:
                index.append(int(row[0]))
                sep.append(str(row[1]))
        return index, sep
    
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
            profile = pd.read_csv("{}/fes_{}.dat".format(self.fepdir, f), delim_whitespace=True, names = ["x","y"])
            profmin = profile.index[profile['y'] == 0.0]
            profmin_value = profile.iloc[profmin[0],0]
            p1 = profile.iloc[:profmin[0],:]
            p2 = profile.iloc[profmin[0]:,:]
            p = p2.append([p1])
            p = p.to_numpy()

            maxima = argrelextrema(p[:,1], np.greater, order=self.order_max)
            minima = argrelextrema(p[:,1], np.less, order=self.order_min)

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
    
    def run(self):
        """
        Function that converts the torsion angle values into corresponding letters, counts the occurance of individual conformer strings and assesses statistics by performing block averages.
        
        Parameters
        ----------

        Returns
        ------
        binary: dataframe
            shape of inputfile, with letters replacing the torsion angle values and including separators
        count: dataframe
            counting how often a conformer string occured
        """

        self.binary, self.count = self._create_binary()
        self._perform_block_averages()
        
        return self.binary, self.count

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
    
        colvar = self._readinputfile()
        colvar = colvar.iloc[::round(colvar.shape[0]/self.length), :]

        binary = self._readinputfile()
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
    
        self.angles_separators = list(binary.columns)
        count = binary.groupby(self.angles_separators).size().to_frame(name='Count').reset_index()
        count['Conformer'] = count[count.columns[0:len(self.angles_separators)]].apply(lambda x: ''.join(x.dropna().astype(str)), axis=1)
        count = count.drop(binary.columns[0:len(self.angles_separators)], axis=1)
        count = count.sort_values("Count", ascending=False, ignore_index=True)
     
        self.binary_compressed = binary.copy()
        self.binary_compressed['Conformer'] = self.binary_compressed[self.binary_compressed.columns[0:len(self.binary_compressed.columns)]].apply(lambda x: ''.join(x.dropna().astype(str)),axis=1)
        self.binary_compressed = self.binary_compressed.drop(self.binary_compressed.columns[0:len(self.binary_compressed.columns)-1], axis=1)
        
        return binary, count
    
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
        outputdir: str
            Directory name to store the conformer files in

        Returns
        -------
        """
        features = list(self.binary.columns)
        count = self.count.set_index('Conformer')
        step = len(self.binary)/10

        clear1 = np.arange(0,len(self.binary), step)
        clear2 = np.arange(step, len(self.binary) + step, step)
        clear1 = list(map(int, clear1))
        clear2 = list(map(int, clear2))

        for i,j,k in zip(range(1,11), clear1, clear2):
            binary_short = self.binary.iloc[j:k,:]
            bbinary_short = binary_short.groupby(features).size().to_frame(name = 'Count').reset_index()
            norm = len(binary_short)
            p = (bbinary_short.iloc[:,len(features)] / norm)
            bbinary_short["Probability"] = p
            bbinary_short['Conformer'] = bbinary_short[bbinary_short.columns[0:len(features)]].apply(lambda x: ''.join(x.dropna().astype(str)),axis=1)
            bbinary_short = bbinary_short.drop(bbinary_short.columns[0:len(features)], axis=1)
            # Join both dataframes taking 'Conformer' as key column for comparison
            bbinary_short = bbinary_short.set_index('Conformer')

            bbbinary_short = bbinary_short.join(count,how='outer',lsuffix='_block', rsuffix='_total')
            bbbinary_short = bbbinary_short.round(6)
            bbbinary_short = bbbinary_short.fillna(0.0)
            bbbinary_short = bbbinary_short.reset_index()

            # Print partially data sorted to file
            bbbinary_short.to_csv("{}/Cluster_conformer{}.dat".format(self.outputdir,i), sep= " ", header=True)
    
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
        threshold: float
            Propability limit up to which conformers are included in the graph 
        ymax: float
            Maximum height of the y-axis (up to 100, as we plot probabilites), default is 100%
        size: int
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

    def validate_fep(self):

        """
        Function that plots free energy profiles with annotated minima and maxima.

        Evaluation of how good the identification of minima and maxima worked.
        -------
        """

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

    def pca(self,
            components = 2,
            ranks = 3,
            components_plot = [1,2], #only 2D supported 
#
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
        #
        #colvar_scaled = StandardScaler().fit_transform(colvar_pca)
        # Choose number of principle components
        pca = PCA(n_components=components)
        # Compute PCA and transform into dataframe with target addition
        principalComponents = pca.fit_transform(colvar_pca)
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
                pass
            #plot all
            ax.scatter(pca_df_scaled.iloc[:,0], pca_df_scaled.iloc[:,1], c = "gray", alpha = 0.3, marker = ".", s = 5, label = "all")
            #plot conformers
            for i, t in enumerate(top):
                indicesToKeep = pca_df_scaled.loc[:,"Conformer"] == t
                ax.scatter(pca_df_scaled.loc[indicesToKeep,0], pca_df_scaled.loc[indicesToKeep,1]
                           , c = color[i]
                           , s = 30
                           , label = label[i]
                           , edgecolors = "black"
                           , marker = "."
                           , linewidth = 0.2
                           , alpha = 1
                               )
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
            #plot all
            ax.scatter(finalDf.iloc[:,components_plot[0]-1], finalDf.iloc[:,components_plot[1]-1], c = "darkgray", alpha = 0.3, marker = ".", s = 5, label = "all")
            #plot conformers
            for i, t in enumerate(top):
                indicesToKeep = finalDf.loc[:,"Conformer"] == t
                ax.scatter(finalDf.loc[indicesToKeep,components_plot[0]-1], finalDf.loc[indicesToKeep,components_plot[1]-1]
                           , c = color[i]
                           , s = 30
                           , label = label[i]
                           , edgecolors = "black"
                           , marker = "."
                           , linewidth = 0.2
                           , alpha = 1
                               )
            if pick == True:
                for d in datatopick:
                    ax.scatter(finalDf.iloc[d,components_plot[0]-1], finalDf.iloc[d,components_plot[1]-1], c = colorpick, marker = "x", s = 5)
            else:
                pass
            
            if legend == True:
                ax.legend(fontsize = fontsize)
            else:
                pass
                
        if file is None:
            pass
        else:
            plt.savefig(file, bbox_inches='tight')
        plt.show()

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