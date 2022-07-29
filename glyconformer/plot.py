import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from glyconformer import process_files

def roundtriptime(replica, path):
    """
    Function that plots the (average)round trip time for each replica along the progression of time.
    
    Parameters
    ----------
    folder : str
        path to the desired directory 
        
    Returns
    -------
    """
    for i in range(0,replica):
        rtt0 = pd.read_csv("{}/rtt_rep{}.dat".format(path,i), delim_whitespace = True, comment='#')
        rtt0 = rtt0.index.tolist()
        average = round(sum(rtt0)/len(rtt0),2)

        progress = np.arange(0,len(rtt0))
        plt.rcParams['figure.dpi'] = 300
        plt.plot(progress, rtt0, "x")
        plt.axhline(average, linestyle = "--", c = "red")
        plt.xticks([])
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.ylabel("Duration [ns]", fontsize = 20)
        plt.text(3, 160, "Average rtt: {}".format(average), fontsize = 20)
        plt.title("Replica {}".format(i), fontsize = 20)
        plt.ylim(-10,200)
        plt.show()
    
def free_energy(path):
    """
    Function that plots the one-dimensional free energy of a torsion angle.
    
    The progression of the surface over time is plotted in different colours. 
    
    Parameters
    ----------
    path : str
        directory path leading to the final directory containing the fes.dat files  
    
    Returns
    -------
    """
    
    for filen in glob.glob(path):

        temp_dirList = process_files._read_in_folder(filen,'')
        dirList = []
        for folder in temp_dirList:
            dirList.append(filen+folder)
        del temp_dirList

        fig = plt.figure( figsize=(16, 9), facecolor='w')
        plt.subplots_adjust(left=0.16, bottom=0.15, right=0.9, top=None, wspace=None, hspace=0.05)

        counter = 0
        for ii in dirList:

            plt.rcParams['figure.dpi'] = 300
            col = np.loadtxt(ii, unpack=True, usecols=[0,1])
            plot = plt.plot(col[0], col[1],'-', linewidth=2, zorder=2, alpha=1, color = plt.cm.jet(counter/(len(dirList)-1)))

            plt.xlabel(r'Torsion angle [rad]', fontsize=25)
            plt.ylabel(r'Free energy (kJ/mol)', fontsize=25)

            ax1 = plt.gca()
            for tick in ax1.xaxis.get_major_ticks():
                tick.label1.set_fontsize(20)
            for tick in ax1.yaxis.get_major_ticks():
                tick.label1.set_fontsize(20)

            counter += 1

        plt.show()

def twodim_pucker(path, file, shape):
    """
    Function that plots a two-dimensional pucker free energy profile.
    
    Plots are along variables phi and theta of the Cremer-Pople coordinates. 
    Phi is plotted along the x-axis and theta along the y-axis. 
    The function was only test on files produced by modules of PLUMED.
    
    Parameters
    ----------
    path : str
        directory path to the input file
    file : str
        colvar file that contains the data points x,y and z in the first three consecutive columns
    shape : int
        final shape in x and y direction of the dataset , e.g.: [200,201]
    
    Returns
    -------
    """
    
    x,y,z = np.loadtxt("{}/{}".format(path,file), usecols=(0,1,2), unpack=True)

    X = np.reshape(x, (shape[0], shape[1]))
    Y = np.reshape(y, (shape[0], shape[1]))
    Z_unshifted = np.reshape(z, (shape[0], shape[1]))
    Z_minimum = Z_unshifted.min()

    xnew = np.linspace(-np.pi/2., np.pi/2.,shape[1])
    for i in range(len(X)):
        X[i] = xnew

    #Shift minimum to = 0 in FEP
    if Z_minimum > 0:
        Z = Z_unshifted - Z_minimum
    elif Z_minimum == 0:
        print("Data is already shifted to 0 for its minimum value")
        Z = Z_unshifted
    else:
        Z = Z_unshifted + abs(Z_minimum)
    Z_min = Z.min()
    Z_max = Z.max()

    # Plot 2D free energy surface with heatmap and legend 
    plt.rcParams['figure.dpi'] = 600
    fig = plt.figure()

    ax = fig.add_subplot(111, projection='mollweide')
    im = ax.pcolormesh(Y, X, Z, cmap=plt.cm.jet, vmin = 0, vmax = 40)

    ax.plot([0, 0], [-1.04, 1.04], color='black', linewidth = 0.5)
    ax.plot([0.523599, 0.523599], [-1.04, 1.04], color='black', linewidth = 0.5)
    ax.plot([-0.523599, -0.523599], [-1.04, 1.04], color='black', linewidth = 0.5)
    ax.plot([1.04, 1.04], [-1.04, 1.04], color='black', linewidth = 0.5)
    ax.plot([-1.04, -1.04], [-1.04, 1.04], color='black', linewidth = 0.5)
    ax.plot([1.57, 1.57], [-1.04, 1.04], color='black', linewidth = 0.5)
    ax.plot([-1.57, -1.57], [-1.04, 1.04], color='black', linewidth = 0.5)
    ax.plot([2.09, 2.09], [-1.04, 1.04], color='black', linewidth = 0.5)
    ax.plot([-2.09, -2.09], [-1.04, 1.04], color='black', linewidth = 0.5)
    ax.plot([2.61, 2.61], [-1.04, 1.04], color='black', linewidth = 0.5)
    ax.plot([-2.61, -2.61], [-1.04, 1.04], color='black', linewidth = 0.5)

    ax.plot([-np.pi, np.pi], [1.57, 1.57], color='black', linewidth = 0.5)
    ax.plot([-np.pi, np.pi], [1.04, 1.04], color='black', linewidth = 0.5)
    ax.plot([-np.pi, np.pi], [0.523, 0.523], color='black', linewidth = 0.5)
    ax.plot([-np.pi, np.pi], [0, 0], color='black', linewidth = 0.5)
    ax.plot([-np.pi, np.pi], [-0.523, -0.523], color='black', linewidth = 0.5)    
    ax.plot([-np.pi, np.pi], [-1.04, -1.04], color='black', linewidth = 0.5)
    ax.plot([-np.pi, np.pi], [-1.57, -1.57], color='black', linewidth = 0.5)

    cbar = plt.colorbar(im)
    cbar.set_label("kJ/mol", fontsize = 15)

    plt.xticks([])
    plt.yticks([])
    plt.show()

    
def graph_eval(input_dir, features, maxima_dict, cluster_dict, label_dict):
    """
    Function that plots free energy profiles with annotated minima and maxima.
    
    Evaluation of how good the identification of minima and maxima worked.
    
    Parameters
    ----------
    input_dir : str
        directory from which to read the files for the different torsion angle arrays
    features : str
        list of torsion angle names in the correct order
    maxima_dict : dict
        dictionary that stores the maxima for each feature
    cluster_dict : dict
        dictionary that holds the minima for each feature
    label_dict : dict
        dictionary with features as keys and labels for each minima as values 
    
    Returns
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
    
    fig, axs = plt.subplots(len(features), figsize = (6,60), constrained_layout=True, )
    
    for i, f in zip(range(len(features)),features):
        profile = pd.read_csv("{}/fes_{}.dat".format(input_dir, f), delim_whitespace=True, names = ["x","y"])
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
        if len(maxima_dict["{}".format(f)]) == 1:
            axs[i].axvline(maxima_dict["{}".format(f)][0], c = co, linestyle = "--")
            axs[i].text(cluster_dict["{}".format(f)][0],30, label_dict["{}".format(f)][0], fontsize = 12, c = co, weight='bold')
        elif len(maxima_dict["{}".format(f)]) == 2:
            axs[i].axvline(maxima_dict["{}".format(f)][0], c = co, linestyle = "--")
            axs[i].axvline(maxima_dict["{}".format(f)][1], c = co, linestyle = "--")
            axs[i].text(cluster_dict["{}".format(f)][0],30, label_dict["{}".format(f)][0], fontsize = 12, c = co, weight='bold')
            axs[i].text(cluster_dict["{}".format(f)][1],30, label_dict["{}".format(f)][1], fontsize = 12, c = co, weight='bold')
            axs[i].text(cluster_dict["{}".format(f)][2],30, label_dict["{}".format(f)][2], fontsize = 12, c = co, weight='bold')
    
        elif len(maxima_dict["{}".format(f)]) == 3:
            axs[i].axvline(maxima_dict["{}".format(f)][0], c = co, linestyle = "--")
            axs[i].axvline(maxima_dict["{}".format(f)][1], c = co, linestyle = "--")
            axs[i].axvline(maxima_dict["{}".format(f)][2], c = co, linestyle = "--")
            axs[i].text(cluster_dict["{}".format(f)][0],30, label_dict["{}".format(f)][0], fontsize = 12, c = co, weight='bold')
            axs[i].text(cluster_dict["{}".format(f)][1],30, label_dict["{}".format(f)][1], fontsize = 12, c = co, weight='bold')
            axs[i].text(cluster_dict["{}".format(f)][2],30, label_dict["{}".format(f)][2], fontsize = 12, c = co, weight='bold')
            axs[i].text(cluster_dict["{}".format(f)][3],30, label_dict["{}".format(f)][3], fontsize = 12, c = co, weight='bold')

    plt.show()
    
def convergence(binary_path,window):
    """
    Function that computes the moving average and plots the probability of the three most populated conformers over time, accessing convergence.
    
    Parameters
    ----------
    binary_path : str
        complete path to binary COLVAR file, including the file
    window : int
        length of the window used for computing the average
    
    Returns
    -------
    """
    
    binary = pd.read_csv(binary_path, delim_whitespace = True)
    binary['Conformer'] = binary[binary.columns[0:len(binary.columns)]].apply(lambda x: ''.join(x.dropna().astype(str)),axis=1)
    binary = binary.drop(binary.columns[0:len(binary.columns)-1], axis=1)
    
    count = binary.groupby('Conformer').size().to_frame(name = 'count').reset_index()
    count = count.sort_values(by=['count'], ascending=False, ignore_index = True)
    top = [count.loc[0, "Conformer"], count.loc[1, "Conformer"], count.loc[2, "Conformer"]]
    
    label = ["1st","2nd","3rd"]
    color = ["#173c4d","#146b65", "#4e9973"]
    
    plt.rcParams['figure.dpi'] = 600
    plt.figure()
    plt.xlabel("Time [µs]", fontsize=15)
    plt.ylabel("Probability [%]", fontsize=15)
    plt.tick_params(axis='both', which='major', labelsize=15)
    
    plt.ylim(0,100)
    
    for n in (range(0,3)):
        
        binary = pd.read_csv(binary_path, delim_whitespace = True)
        binary['Conformer'] = binary[binary.columns[0:len(binary.columns)]].apply(lambda x: ''.join(x.dropna().astype(str)),axis=1)
        binary = binary.drop(binary.columns[0:len(binary.columns)-1], axis=1)
        
        for i in binary.index:
            if binary.loc[i,"Conformer"] == top[n]:
                pass
            else:
                binary.loc[i,"Conformer"] = np.nan

        binary["rollcount"] = binary.rolling(window, min_periods=1).count()
        binary.loc[:, "rollcount"] = [x / window for x in binary["rollcount"]]
        binary.loc[:, "rollcount"] = [x * 100 for x in binary["rollcount"]]
        plt.plot(binary.index/125, binary.loc[:, "rollcount"], label = label[n], c = color[n], linestyle = "-")
        
    plt.legend(fontsize = 15)
    plt.show()
