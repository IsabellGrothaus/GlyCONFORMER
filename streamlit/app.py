import csv
import importlib
import streamlit as st
from glyconformer.lib import Glyconformer, Glycompare
import os as os
import pandas as pd


# -------- INITIALIZATION --------- #


def setSessionStates():

    if 'glycan_state' not in st.session_state:
        st.session_state['glycan_state'] = None

setSessionStates()


# -------- MAIN PAGE --------- #

tab_main1, tab_main2, tab_main3 = st.tabs(["test1", "test2", "test3"])

with tab_main1:
    st.header("Headline")


# -------- FUNCTIONS--------- #


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

@st.cache_data                              # wenn der Glykan-String bereits abgerufen wurde, bezieht er die Daten aus dem Cache
def initialize_Glycan(glycantype: str):
    st.write(glycantype)

    Glycan = Glyconformer(
        inputfile =         "../TUTORIAL/{}_example/{}_angles.dat".format(glycantype, glycantype), 
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

def on_change(selected_glycan: str):
    
    st.session_state['glycan_state'] = selected_glycan

    with tab_main1: 
        initialize_Glycan(st.session_state['glycan_state'])


# -------- SIDEBAR --------- #


with st.sidebar:

    tab_sidebar1, tab_sidebar2 = st.tabs(["User Input", "Preselection"])

    with tab_sidebar1:
        st.subheader("select your dataframe")
        uploaded_file = st.file_uploader(
                                    "...",
                                    label_visibility = "collapsed",
        )

        if uploaded_file is not None:
            st.write("filename:", uploaded_file.name)
            file_content = pd.read_csv(
                                    uploaded_file, 
                                    sep='\s+',                           # \s+' -> mit Whitespace gertrennt. + -> (ein oder mehrmals)
                                    header=None, 
                                    names=['col1', 'col2'], 
                                    dtype={'col1': int, 'col2': str}
            )
            st.write(file_content)                                       # testimplementation
            st.write(file_content['col1'].to_numpy(), file_content['col2'].to_numpy())

    with tab_sidebar2:
        st.subheader("select a dataset")

        selected_glycan = st.selectbox(
                                    "...",
                                    ("M5", "A2G2S2"), 
                                    index = None,                                           # initialisiert eine leere Auswahlbox
                                    label_visibility = "collapsed",
                                    placeholder = "choose one of the following glycans"
        )

        if selected_glycan is not None:
            st.write("You selected:", selected_glycan)
            on_change(selected_glycan)


# -------- ### --------- #


st.write("""
    <style>
        .stTabs > div > div > div > button {
            width: 50%;
        }
    </style>
""", unsafe_allow_html=True)



'''
if "attendance" not in st.session_state:
    st.session_state.attendance = set()

def take_attendance():
    if st.session_state.name in st.session_state.attendance:
        st.info(f"{st.session_state.name} has already been counted.")
    else:
        st.session_state.attendance.add(st.session_state.name)


with st.form(key="my_form"):
    st.text_input("Name", key="name")
    st.form_submit_button("I'm here!", on_click=take_attendance)
'''