import csv
import importlib
import streamlit as st
from glyconformer.lib import Glyconformer, Glycompare
import numpy as np
import pandas as pd
import plumed
import os as os




st.title("Hi, I am a streamlit Web-App")

def loadLocalGlycans(): 

    local_glycans = []

    for i in os.listdir("../LIBRARY_GLYCANS/"):
        if("__" not in i):
            local_glycans.append(i)

    return local_glycans

def readAnglesData(path, file):
    feature = importlib.resources.read_text(path, file)
    feature = feature.split()

    return feature

def readSeparatorData(path, file):
    dictionary = {'separator_index': [], 'separator': []}

    csv_reader = csv.reader(importlib.resources.open_text(path, file), delimiter=' ')
    
    for row in csv_reader:
        dictionary['separator_index'].append(int(row[0]))
        dictionary['separator'].append(str(row[1]))

    return dictionary

def initialize_Glycan(glycantype: str):
    st.write(glycantype)

    Glycan = Glyconformer(
        inputfile =         "../TUTORIAL/{}_example/{}_angles.dat".format(glycantype,glycantype), 
        angles =            readAnglesData("LIBRARY_GLYCANS.{}".format(glycantype), "angles.dat"), 
        omega_angles =      readAnglesData("LIBRARY_GLYCANS.{}".format(glycantype), "omega_angles.dat"), 
        separator_index =   readSeparatorData("LIBRARY_GLYCANS.{}".format(glycantype), "separator.dat")['separator_index'], 
        separator =         readSeparatorData("LIBRARY_GLYCANS.{}".format(glycantype), "separator.dat")['separator'], 
        fepdir =            "../LIBRARY_GLYCANS/{}".format(glycantype), 
        order_min =         5, 
        order_max =         5,

        length = 20, 
        weights = None
    )



    st.write(Glycan.binary)
    st.write(Glycan.minima["phi1_2"][0])
    st.write(Glycan.separator_index)
    st.write(Glycan.separator)


if 'glycan_state' not in st.session_state:
    st.session_state['glycan_state'] = None

def on_change(selected_glycan):
    if st.session_state['glycan_state'] != selected_glycan:

        print(selected_glycan)
        # try:
        #    position = loadLocalGlycans().index(selected_glycan)
        #    selected_glycan = loadLocalGlycans()[position]
        #except ValueError:
        #    print("Entry can not be found")

        st.session_state['glycan_state'] = selected_glycan
        initialize_Glycan(selected_glycan)

selected_glycan = st.selectbox("Wählen Sie eine der folgenden Optionen", loadLocalGlycans(), on_change=lambda: on_change(selected_glycan))




# colvar = plumed.read_as_pandas("../TUTORIAL/M5_example/M5_angles.dat")
# colvar = colvar[readData("LIBRARY_GLYCANS.{}".format(glycantype), "angles.dat")]
# st.write(readData("LIBRARY_GLYCANS.{}".format(glycantype), "angles.dat"))
# colvar = readColvar("TUTORIAL/{}_example/{}_angles.dat".format(glycantype, glycantype), readData("LIBRARY_GLYCANS.{}".format(glycantype), "angles.dat"))

