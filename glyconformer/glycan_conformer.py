"""
This module contains functions to generate N-glycan
conformer strings from a COLVAR file containing
torsion angle values recorded over a MD simulation.

List of functions:
conformer_limit_string
create_binary
find_min_max
label_min
load_dict
perform_block_averages
plot_distribution
safe_dict
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plumed
from scipy.signal import argrelextrema
import json
import matplotlib.ticker as ticker

def _combine_dict(d1, d2, d3, d4):
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

def safe_dict(dict_name, filename):   
    """
    Function saving a dictionary to file.
    
    Parameters
    ----------
    dict_name : dict
        Name of the dictionary to safe
    filename : str
        Path + name of the file with file extension     
    
    Returns
    -------
    """
    
    with open('{}'.format(filename), 'w') as file:
        file.write(json.dumps(dict_name))
        
def load_dict(filename):
    """
    Function loading a dictionary from file.
    
    Parameters
    ----------
    filename : str
        Name of the path + file to read in,
        including file extension 
    
    Returns
    -------
    dict
        Variable holding the dictionary
    """
    
    with open('{}'.format(filename)) as f:
        data = f.read()
        print(data)
    dict = json.loads(data)
    
    return dict

def find_min_max(input_dir, features, order_max, order_min):
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
    input_dir : str
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
        profile = pd.read_csv("{}/fes_{}.dat".format(input_dir, f), delim_whitespace=True, names = ["x","y"])
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

def label_min(cd, features, f_omega):
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


    for f in features:
        if f in f_omega:
            for key, value in oangle_dict.items():
                if value[0] <= cd["{}".format(f)][0] < value[1]:
                    name_dict["{}".format(f)] = key
                if cd["{}".format(f)][0] == 4:
                    name_dict["{}".format(f)] = 'none'      
        else:
            for key, value in tangle_dict.items():
                if value[0] <= cd["{}".format(f)][0] < value[1]:
                    name_dict["{}".format(f)] = key
                if cd["{}".format(f)][0] == 4:
                    name_dict["{}".format(f)] = 'none'

    for f in features:
        if len(cd["{}".format(f)]) > 1:
            if f in f_omega:
                for key, value in oangle_dict.items():
                    if value[0] <= cd["{}".format(f)][1] < value[1]:                    
                        name_dict1["{}".format(f)] = key
                    if cd["{}".format(f)][1] == 4:
                        name_dict1["{}".format(f)] = 'none'        
            else:
                for key,value in tangle_dict.items():
                    if value[0] <= cd["{}".format(f)][1] < value[1]:                    
                        name_dict1["{}".format(f)] = key
                    if cd["{}".format(f)][1] == 4:
                        name_dict1["{}".format(f)] = 'none'

    for f in features:                
        if len(cd["{}".format(f)]) > 2:
            if f in f_omega:
                for key, value in oangle_dict.items():
                    if value[0] <= cd["{}".format(f)][2] < value[1]:                    
                        name_dict2["{}".format(f)] = key
                    if cd["{}".format(f)][2] == 4:
                        name_dict2["{}".format(f)] = 'none'                    
            else:
                for key,value in tangle_dict.items():
                    if value[0] <= cd["{}".format(f)][2] < value[1]:                    
                        name_dict2["{}".format(f)] = key
                    if cd["{}".format(f)][2] == 4:
                        name_dict2["{}".format(f)] = 'none'

    for f in features:                
        if len(cd["{}".format(f)]) > 3:
            if f in f_omega:
                for key, value in oangle_dict.items():
                    if value[0] <= cd["{}".format(f)][3] < value[1]:                    
                        name_dict3["{}".format(f)] = key
            else:
                for key,value in tangle_dict.items():
                    if value[0] <= cd["{}".format(f)][3] < value[1]:                    
                        name_dict3["{}".format(f)] = key

    label_dict = _combine_dict(name_dict,name_dict1,name_dict2,name_dict3)
    return label_dict

def create_binary(maxima_dict, label_dict, colvar_dir, input_dir, 
                  colvar_length, features,* , loc1, loc2, loc3 = 0, 
                  loc4 = 0, loc5 = 0, loc6 = 0, locF = 0, 
                  locG = 0, bra3 = "───", bra4 = "───", bra5 = "───", 
                  bra6 = "───"):
    """
    Converts torsion angle values read from COVLAR file into IUPAC letters.
    
    Function that reads in a COLVAR file with torsion angle values stored 
    in each column and replaces each value by the corresponding IUPAC label, 
    which was previously defined in the label_dict.
    
    Parameters
    ----------
    maxima_dict : dict
        Dictionary that contains the maxima for each feature
    label_dict : dict
        Dictionary with features as keys and labels for each minima as values 
    colvar_dir : str
        Path to the COLVAR file with stored features as columns
    colvar_length : int
        Desired length of the COLVAR file
    features : str
        List of torsion angle names (= column names) in the correct order
    glycantype : str
        Define type of glycan from "complex", "mannose", or "hybride"
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
    features : str
        Updated feature list with branch separators
    branches : str
        List of branches present in the processed glycan structure
    """
   
    colvar = plumed.read_as_pandas(colvar_dir)
    colvar = colvar[features]
    colvar = colvar.iloc[::round(colvar.shape[0]/colvar_length), :]

    c = plumed.read_as_pandas(colvar_dir)
    c = colvar[features]
    c = colvar.iloc[::round(colvar.shape[0]/colvar_length), :]

    for f in features:
        if len(maxima_dict["{}".format(f)]) == 1:

            c["{}".format(f)] = np.where((-3.5 <= colvar["{}".format(f)]) & (colvar["{}".format(f)] <= 3.5) , label_dict["{}".format(f)][0], c["{}".format(f)])
            
 
        elif len(maxima_dict["{}".format(f)]) == 2:

            c["{}".format(f)] = np.where((-3.5 <= colvar["{}".format(f)]) & (colvar["{}".format(f)] < maxima_dict["{}".format(f)][0]) , label_dict["{}".format(f)][0], c["{}".format(f)])
            c["{}".format(f)] = np.where((maxima_dict["{}".format(f)][0] <= colvar["{}".format(f)]) & (colvar["{}".format(f)] < maxima_dict["{}".format(f)][1]) , label_dict["{}".format(f)][1], c["{}".format(f)])
            c["{}".format(f)] = np.where((maxima_dict["{}".format(f)][1] <= colvar["{}".format(f)]) & (colvar["{}".format(f)] <= 3.5) , label_dict["{}".format(f)][2], c["{}".format(f)])      
            
        elif len(maxima_dict["{}".format(f)]) == 3:

            c["{}".format(f)] = np.where((-3.5 <= colvar["{}".format(f)]) & (colvar["{}".format(f)] < maxima_dict["{}".format(f)][0]) , label_dict["{}".format(f)][0], c["{}".format(f)])
            c["{}".format(f)] = np.where((maxima_dict["{}".format(f)][0] <= colvar["{}".format(f)]) & (colvar["{}".format(f)] < maxima_dict["{}".format(f)][1]) , label_dict["{}".format(f)][1], c["{}".format(f)])
            c["{}".format(f)] = np.where((maxima_dict["{}".format(f)][1] <= colvar["{}".format(f)]) & (colvar["{}".format(f)] < maxima_dict["{}".format(f)][2]) , label_dict["{}".format(f)][2], c["{}".format(f)])
            c["{}".format(f)] = np.where((maxima_dict["{}".format(f)][2] <= colvar["{}".format(f)]) & (colvar["{}".format(f)] <= 3.5) , label_dict["{}".format(f)][3], c["{}".format(f)])
   
            
    c = pd.DataFrame(c, columns = c.columns)
    c = c.replace(' t ', ' T ')
    c = c.replace(' TG', ' tg')
    
    branches = []
    c.insert(loc=loc1, column='sep1', value="6──")
    branches.append("6──")
    c.insert(loc=loc2+1, column='sep2', value="3──")
    branches.append("3──")
    
    if loc3 == 0:
        pass
    else:
        c.insert(loc=loc3, column='sep3', value=bra3)
        branches.append(bra3)
    
    if loc4 == 0:
        pass
    else:
        if loc3 != 0:
            c.insert(loc=loc4, column='sep4', value=bra4)
            branches.append(bra4)
        else:
            raise ValueError("loc3 need to be defined before loc4")
            
    if loc5 == 0:
        pass
    else:
        if loc4 != 0:
            c.insert(loc=loc5, column='sep5', value=bra5)
            branches.append(bra5)
        else:
            raise ValueError("loc4 need to be defined before loc5")
    
    if loc6 == 0:
        pass
    else:
        if loc6 != 0:
            c.insert(loc=loc6, column='sep6', value=bra6)
            branches.append(bra6)
        else:
            raise ValueError("loc5 need to be defined before loc6")
    
    if locF == 0:
        pass
    else:
        c.insert(loc=locF, column='sepF', value="f──")
        branches.append("f──")
        
    if locG == 0:
        pass
    else:
        c.insert(loc=locG, column='sepG', value="g──")
        branches.append("g──")
        
    features = list(c.columns)
    ccf = c.groupby(features).size().to_frame(name = 'count').reset_index()
    ccf['Conformer'] = ccf[ccf.columns[0:len(features)]].apply(lambda x: ''.join(x.dropna().astype(str)),axis=1)
    ccf = ccf.drop(c.columns[0:len(features)], axis=1)
    c.to_csv("{}/COLVAR_binary".format(input_dir) ,sep = " ", index = False)
    
    return c, ccf, features, branches
    
def perform_block_averages(c, ccf, store_dir):
    """
    Calculates probability of conformers for 10 separate blocks.
    
    Function that splits the torsion angle COLVAR dataset in 10 
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
        cc_short = c_short.groupby(features).size().to_frame(name = 'count').reset_index()
        norm = len(c_short)
        p = (cc_short.iloc[:,len(features)] / norm)
        cc_short["Probability"] = p
        cc_short['Conformer'] = cc_short[cc_short.columns[0:len(features)]].apply(lambda x: ''.join(x.dropna().astype(str)),axis=1)
        cc_short = cc_short.drop(cc_short.columns[0:len(features)], axis=1)
        # Join both dataframes taking 'Conformer' as key column for comparison
        cc_short = cc_short.set_index('Conformer')

        ccc_short = cc_short.join(ccf,how='outer',lsuffix='_partial', rsuffix='_full')
        ccc_short = ccc_short.round(6)
        ccc_short = ccc_short.fillna(0.0)
        ccc_short = ccc_short.reset_index()

        # Print partially data sorted to file
        ccc_short.to_csv("{}/Cluster_conformer{}.dat".format(store_dir,i), sep= " ", header=False)

def plot_distribution(input_dir, limit, branches, features):
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
        
    Returns
    -------
    """

    hist1 = pd.read_csv("{}/Cluster_conformer1.dat".format(input_dir), delim_whitespace=True, names = ["Index", "Conformer", "Count_partial", "Prob", "Count_full"] )
    N, average, average2 = 1, hist1.iloc[:,3], hist1.iloc[:,3]*hist1.iloc[:,3]
    for i in range(2,11): 
        histn = pd.read_csv("{}/Cluster_conformer{}.dat".format(input_dir, i), delim_whitespace=True, names = ["Index", "Conformer", "Count_partial", "Prob", "Count_full"])
        N, average, average2 = N + 1, average + histn.iloc[:,3], average2 + histn.iloc[:,3]*histn.iloc[:,3]
    

    average = pd.DataFrame(average/N, columns=['Prob'])
    average = average.drop(average[average['Prob'] < limit ].index)
    average = average.sort_values(by='Prob',ascending=False)
    indexlist = average.index.tolist()

    hist = pd.read_csv("{}/Cluster_conformer1.dat".format(input_dir), names = ["Index", "Conformer", "Count_partial", "Prob", "Count_full"], sep = " ", dtype = str)
    hist = hist.filter(indexlist, axis = 0)
    namelist = hist['Conformer'].tolist()
    name_list = _include_branch_conformer_string(namelist, branches, features)
    name_list = _vertical_conformer_string(namelist)
    average2 = average2.to_frame()
    average2.columns = ["Error"]
    average2 = average2.filter(indexlist, axis = 0)

    # Final variances
    var = (N/(N-1))*( average2.Error / N - average.Prob*average.Prob ) 
    # Errors
    error = np.sqrt( var / N )
 
    pos_list = np.arange(len(name_list))
    
    colors = ["#173c4d","#146b65","#4e9973","#a7c09f","#dfa790","#c76156","#9a2b4c","#600b4a"]
    plt.rcParams['figure.dpi'] = 300
    ax = plt.axes()
    ax.xaxis.set_major_locator(ticker.FixedLocator((pos_list)))
    ax.xaxis.set_major_formatter(ticker.FixedFormatter((name_list)))   
    bar = plt.bar(pos_list, average.Prob * 100, color = "blue", linewidth = "1", yerr = error * 100)
    for i in range(len(pos_list)):
        bar[i].set_color(colors[i])
    
    plt.xlabel("Conformer")
    plt.ylabel('Probability [%]')
    plt.savefig("{}/Conformer_population_full.png".format(input_dir))

    
def conformer_limit_string(limit, directory):
    
    """
    Prints out the conform strings which are occuring more or less than the limit.
    
    Reads in constructed conformer probability records (Cluster_conformer*d.dat files) 
    from 10 blocks and averages them to sort conformers by their occurance. The limit 
    sets which conformers are discarded (below) and which are included for analysis (above).
    
    Parameters
    ----------
    limit: float
        Propability limit up to which conformers are excluded from coloured representation in DimRed
    directory: str
        Name of directory path from which to read the Cluster_conformer.dat files
        
    Returns
    -------
    targets_below: str
        List of conformer strings below the limit
    targets_above: str
        List of conformer strings above the limit
    """

    # Load partially data sorted 
    hist1 = pd.read_csv("{}/Cluster_conformer1.dat".format(directory), delim_whitespace=True, names = ["Index","Conformer","Count", "Prob", "Count_all"])
    N, average = 1, hist1.iloc[:,3]
    for i in range(2,11): 
        histn = pd.read_csv("{}/Cluster_conformer{}.dat".format(directory, i), delim_whitespace=True, names = ["Index","Conformer","Count", "Prob", "Count_all"])
        N, average = N + 1, average + histn.iloc[:,3]

    # Sort conformers by probability
    average = average.to_frame()
    average = average.div(N)
     
    # below
    below = average.drop(average[average['Prob'] > limit ].index, inplace=False)
    below = below.sort_values(by='Prob',ascending=False)
    indexlist = below.index.tolist()
    hist = pd.read_csv("{}/Cluster_conformer1.dat".format(directory), sep = " ", names = ["Index","Conformer","Count", "Prob", "Count_all"], dtype = str)
    hist = hist.filter(indexlist, axis = 0)
    targets_below = hist['Conformer'].tolist()

    # above
    above = average.drop(average[average['Prob'] < limit ].index, inplace=False)
    above = above.sort_values(by='Prob',ascending=False)
    indexlist = above.index.tolist()
    # Safe to use no.1 or any other of these files. All have all conformers included
    hist = pd.read_csv("{}/Cluster_conformer1.dat".format(directory), sep = " ", names = ["Index","Conformer","Count", "Prob", "Count_all"], dtype = str)
    hist = hist.filter(indexlist, axis = 0)
    targets_above = hist['Conformer'].tolist()

    return targets_below,targets_above
