import csv
import importlib
import tempfile
import streamlit as st
import plumed # type: ignore
# from glyconformer.lib import Glyconformer, Glycompare
import pandas as pd # type: ignore
import os as os
import sys
from st_aggrid import AgGrid # type: ignore

sys.path.append("/home/eberl/bachelor_project/GlyCONFORMER/glyconformer/")
from lib import Glyconformer

# -------- INITIALIZATION --------- #


def setSessionStates():

    if 'glycan' not in st.session_state:
        st.session_state['glycan'] = None
    if 'progress' not in st.session_state:
        st.session_state['progress'] = []

    if 'dataframe' not in st.session_state:
        st.session_state['dataframe'] = None
    if 'angles' not in st.session_state:
        st.session_state['angles'] = None
    if 'fep' not in st.session_state:
        st.session_state['fep'] = {}

setSessionStates()


# -------- MAIN PAGE --------- #

tab_main1, tab_main2, tab_main3 = st.tabs(["test1", "test2", "test3"])

with tab_main1:
    st.header("Headline")


# -------- local FUNCTIONS--------- #


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
def initialize_local_Glycan(glycantype: str):

    Glycan = Glyconformer(
        inputfile =         "../TUTORIAL/{}_example/{}_angles.dat".format(glycantype, glycantype),

        angles =            readAnglesData("LIBRARY_GLYCANS.{}".format(glycantype), "angles.dat"), 
        omega_angles =      readAnglesData("LIBRARY_GLYCANS.{}".format(glycantype), "omega_angles.dat"), 
        separator_index =   readSeparatorData("LIBRARY_GLYCANS.{}".format(glycantype), "separator.dat")['separator_index'], 
        separator =         readSeparatorData("LIBRARY_GLYCANS.{}".format(glycantype), "separator.dat")['separator'], 
        fepdir =            "../LIBRARY_GLYCANS/{}".format(glycantype), 
        order_max =         5,
        order_min =         5, 
        weights =           None,

    )
 
    st.write(Glycan.colvar)
    st.write(Glycan.minima["phi1_2"][0])
    st.write(Glycan.separator_index)
    st.write(Glycan.separator)

def initialize_custom_Glycan(glycantype: str):

    Glycan = Glyconformer(
        inputfile =         None,

        angles =            st.session_state['angles'], 
        omega_angles =      readAnglesData("LIBRARY_GLYCANS.{}".format(glycantype), "omega_angles.dat"), 
        separator_index =   readSeparatorData("LIBRARY_GLYCANS.{}".format(glycantype), "separator.dat")['separator_index'], 
        separator =         readSeparatorData("LIBRARY_GLYCANS.{}".format(glycantype), "separator.dat")['separator'], 
        fepdir =            None,
        fep_files =         st.session_state['fep'],
        colvar =            st.session_state['dataframe'], 
        order_max =         5,
        order_min =         5, 
        weights =           None,

    )

    st.write(Glycan.colvar)


# -------- local FUNCTIONS--------- #

def on_change(selected_glycan: str):
    
    st.session_state['glycan'] = selected_glycan

    with tab_main1: 
        initialize_local_Glycan(st.session_state['glycan'])

def change_opacity(element_class, index, opacity):

    # script = f"<script>document.getElementsByClassName('{element_class}')[{index}].childNodes[0].style.opacity = '{opacity}%'; console.log(document.getElementsByClassName('{element_class}')[{index}].childNodes[0])</script>"
    # st.markdown(script, unsafe_allow_html=True)

    # st.markdown(
    #     """
    #     <script>
    #        console.log('Hi')
    #     </script>
    #     """,
    #     unsafe_allow_html=True,
    # )

    # st.markdown("<script src='script.js'></script>", unsafe_allow_html=True)
    return None

def custom_subheader(text: str):

    '''
    with open('styles.css') as f:
        if(index in st.session_state['progress']):
            st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
            st.write("Test")

        else:
            st.write("Test2")

            css_content = f.read().replace('--opacity: 50%;', '--opacity: 100%;')
            st.markdown(f"<style>{css_content}</style>", unsafe_allow_html=True)
    '''


    if(int(text[0]) in st.session_state['progress']):
        st.markdown(f"<h3 class='custom-subheader-active'>{text}</h3>", unsafe_allow_html=True)
    else:
        st.markdown(f"<h3 class='custom-subheader-inactive'>{text}</h3>", unsafe_allow_html=True)

def createKeyName(fileName: str):

    fileName = fileName.removeprefix("fes_")
    fileName = fileName.removesuffix(".dat")

    return fileName

def create_fep_directorie():
    st.session_state['fep'].clear()
    for angle in st.session_state['angles']:
        st.session_state['fep'][angle] = []

def extractAngles(colvar):

    angles: list = []
    for i in colvar.columns:
        if i != "time" and "pseudo" not in i:
            angles.append(i)

    return angles

def readDataframe(dataframe):
    
    # In Bytestring konvertieren
    bytes_data: bytes = dataframe.getvalue()

    # In temporäre Datei schreiben
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file.write(bytes_data)
        temp_file_path = temp_file.name

    colvar = plumed.read_as_pandas(temp_file_path)
    st.session_state['angles'] = extractAngles(colvar)

    colvar = colvar[st.session_state['angles']]

    # Temporäre Datei löschen
    os.remove(temp_file_path)

    return colvar

def checkProgress():
    st.session_state['progress'].clear()  
    st.session_state['progress'].append(1)             # Progress-Step "1" ist bereits enthalten


    custom_subheader("1. Select your dataframe")

    if(1 in st.session_state['progress']):

        dataframe = st.file_uploader(
                                        "...",
                                        label_visibility = "collapsed",
        )      

        if dataframe is not None:

            st.session_state['progress'].append(2)

            st.session_state['dataframe'] = readDataframe(dataframe)
   
            with tab_main2:
                 # st.subheader(f"filename: {dataframe.name}")
                 ''' '''

        else:
            if len(st.session_state['progress']) > 1:
                st.session_state['progress'].pop(1)
    

    st.write(st.session_state['progress'])
    custom_subheader("2. Select your fepfiles")

    if(2 in st.session_state['progress']):
       
        fep_files = st.file_uploader(   "t..",
                                        label_visibility = "collapsed",
                                        key = "fep_files",
                                        accept_multiple_files = True,
        )

        create_fep_directorie()

        if fep_files is not None:
            for file in fep_files:
 
                if createKeyName(file.name) in st.session_state['fep'] and len(st.session_state['fep'][createKeyName(file.name)]) == 0:
                    profile = pd.read_csv(file, sep='\s+', names=["x", "y"])
                    st.session_state['fep'][createKeyName(file.name)].append(profile)
                elif len(st.session_state['fep'][createKeyName(file.name)]) != 0:
                    st.error(f'fep_file bereits vorhanden: {createKeyName(file.name)}')
                else:
                    st.error(f'fep_file nicht vorhanden: {createKeyName(file.name)}')

        st.progress(len(st.session_state['fep']) * (1.0/len(st.session_state['angles'])), "progress ...")

    st.write(st.session_state['fep'])



# -------- SIDEBAR --------- #


with st.sidebar:

    tab_sidebar1, tab_sidebar2 = st.tabs(["User Input", "Preselection"])

    with tab_sidebar1:
        checkProgress()


    with tab_sidebar2:
        st.subheader("select a dataset")

        selected_glycan = st.selectbox(
                                    "...",
                                    loadLocalGlycans(), 
                                    index = None,                                           # initialisiert eine leere Auswahlbox
                                    label_visibility = "collapsed",
                                    placeholder = "choose one of the following glycans"
        )

        if selected_glycan is not None:
            on_change(selected_glycan)


# -------- ### --------- #


with open('styles.css') as f:
    st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)

'''
st.write("""
    <style>
        .stTabs > div > div > div > button {
            width: 50%;
        }
    </style>
""", unsafe_allow_html=True)
'''

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

    -------

    file_content = pd.read_csv(
                    dataframe_file, 
                    sep='\s+',              # \s+' -> mit Whitespace gertrennt. + -> (ein oder mehrmals)
                    header=None, 
                    names=['col1', 'col2'], 
                    dtype={'col1': int, 'col2': str}
    )

    st.write(file_content)
    st.write(file_content['col1'].to_numpy(), file_content['col2'].to_numpy())

    -------

    for key, value in st.session_state['fep'].items():
        st.write(key)
        st.write(value)
'''