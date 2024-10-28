import csv
import importlib
import tempfile
import time
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


st.markdown("""
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@4.1.3/dist/css/bootstrap.min.css" integrity="sha384-MCw98/SFnGE8fJT3GXwEOngsV7Zt27NXFoaoApmYm81iuXoPkFOJwJ8ERdknLPMO" crossorigin="anonymous">
""", unsafe_allow_html=True)


def setSessionStates():

    if 'glycan' not in st.session_state:
        st.session_state['glycan'] = None
    if 'progress' not in st.session_state:
        st.session_state['progress'] = []

    if 'dataframe' not in st.session_state:
        st.session_state['dataframe'] = None
    if 'angles' not in st.session_state:
        st.session_state['angles'] = []
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
        fepdir =            "../LIBRARY_GLYCANS/{}".format(glycantype), 

        angles =            readAnglesData("LIBRARY_GLYCANS.{}".format(glycantype), "angles.dat"), 
        omega_angles =      readAnglesData("LIBRARY_GLYCANS.{}".format(glycantype), "omega_angles.dat"), 
        separator_index =   readSeparatorData("LIBRARY_GLYCANS.{}".format(glycantype), "separator.dat")['separator_index'], 
        separator =         readSeparatorData("LIBRARY_GLYCANS.{}".format(glycantype), "separator.dat")['separator'], 
        fep_files =         {},
        order_max =         5,
        order_min =         5, 
        weights =           None,

        colvar =            None,
        length =            None
    )
 
    st.write(Glycan.colvar)
    st.write(Glycan.minima["phi1_2"][0])
    st.write(Glycan.separator_index)
    st.write(Glycan.separator)
    st.pyplot(Glycan.validate_fep())

def initialize_custom_Glycan(glycantype: str):

    Glycan = Glyconformer(
        inputfile =         None,
        fepdir =            None,

        angles =            st.session_state['angles'], 
        omega_angles =      readAnglesData("LIBRARY_GLYCANS.{}".format(glycantype), "omega_angles.dat"), 
        separator_index =   readSeparatorData("LIBRARY_GLYCANS.{}".format(glycantype), "separator.dat")['separator_index'], 
        separator =         readSeparatorData("LIBRARY_GLYCANS.{}".format(glycantype), "separator.dat")['separator'], 
        fep_files =         st.session_state['fep'],
        order_max =         5,
        order_min =         5, 
        weights =           None,

        colvar =            st.session_state['dataframe'], 
        length =            None
    )

    st.write(Glycan.colvar)
    st.write(Glycan.minima["phi1_2"][0])
    st.write(Glycan.separator_index)
    st.write(Glycan.separator)
    st.pyplot(Glycan.validate_fep())

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

def createProgressBar(percentage: str):

    return f"""
        <div class="progress_container">
            <div class='progress'>
                <div class='progress-bar progress-bar-striped progress-bar-animated' role='progressbar' style='width: {percentage}%' aria-valuemin='0' aria-valuemax='100'>{percentage}%</div>
            </div>
        </div>
        """

def rewindProgress(number: int):
    st.session_state['progress'].clear()  

    for i in range(1, number + 1):
        st.session_state['progress'].append(i)

    print(f"----------------------------- {st.session_state['progress']}")

def checkAmountOfExistingFeps():

    count: int = 0

    for key, value in st.session_state['fep'].items():
        if value:        # ['element'] => True, [] => False
            count += 1

    return count

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
        st.markdown(f"<h3 class='custom-subheader-active'>{text} </h3>", unsafe_allow_html=True)
    else:
        st.markdown(f"<h3 class='custom-subheader-inactive'>{text}</h3>", unsafe_allow_html=True)

def createPercentage(startValue: int, endValue: int):

    return int(round(startValue / endValue, 2) * 100)

def createKeyName(fileName: str):

    fileName = fileName.removeprefix("fes_")
    fileName = fileName.removesuffix(".dat")

    return fileName

def create_fep_directorie():
    st.session_state['fep'].clear()
    for angle in st.session_state['angles']:
        st.session_state['fep'][angle] = []

def extractAngles(colvar):

    for i in colvar.columns:
        if i != "time" and "pseudo" not in i:
            st.session_state['angles'].append(i)

def readDataframe(dataframe):
    
    # In Bytestring konvertieren
    bytes_data: bytes = dataframe.getvalue()

    # In temporäre Datei schreiben
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file.write(bytes_data)
        temp_file_path = temp_file.name

    colvar = plumed.read_as_pandas(temp_file_path)
    extractAngles(colvar)                       # füllt st.session_state['angles']

    colvar = colvar[st.session_state['angles']]
    colvar2 = colvar[['phi1_2', 'psi1_2', 'phi2_3', 'psi2_3', 'phi3_5', 'psi3_5', 'omega3_5', 'phi5_7', 'psi5_7', 'omega5_7', 'phi5_6', 'psi5_6', 'phi3_4', 'psi3_4']]
    with tab_main2:
        st.write(colvar)
        st.write(colvar2)

    # Temporäre Datei löschen
    os.remove(temp_file_path)

    return colvar

def checkProgress():
    st.session_state['progress'].clear()  
    st.session_state['progress'].append(1)             # Progress-Step "1" ist bereits enthalten

    # ------- 1 ------- #

    custom_subheader("1. Select your dataframe")

    if(1 in st.session_state['progress']):
        
        my_bar = st.progress(0, text = "0/1")

        file_upload = st.file_uploader(
                                        "...",
                                        label_visibility = "collapsed",
                                        key = "dataframe_upload",                           
        )


        if file_upload is not None:
            if True:
                st.session_state['progress'].append(2)
                my_bar.progress(createPercentage(1, 1), text = "1/1")

                st.session_state['dataframe'] = readDataframe(file_upload)
            else:
                st.error("test")

        else:
            rewindProgress(1)
    

    # ------- 2 ------- #

    custom_subheader("2. Select your fepfiles")

    if(2 in st.session_state['progress']):
       
        create_fep_directorie()
        my_bar = st.progress(0, text = f'{checkAmountOfExistingFeps()}/{len(st.session_state['fep'])}')

        file_upload = st.file_uploader( "...",
                                        label_visibility = "collapsed",
                                        key = "fep_upload",
                                        accept_multiple_files = True,
        )

        if file_upload != []:
            for file in file_upload:

                if createKeyName(file.name) in st.session_state['fep'] and not st.session_state['fep'][createKeyName(file.name)]:       # wenn Winkel im 'dictionary' vorhanden && nicht belegt ist
                    profile = pd.read_csv(file, delim_whitespace=True, names=["x", "y"])
                    st.session_state['fep'][createKeyName(file.name)].append(profile)
                elif createKeyName(file.name) in st.session_state['fep'] and st.session_state['fep'][createKeyName(file.name)]:
                    st.error(f'fep_file bereits vorhanden: {createKeyName(file.name)}')     # wenn Winkel im 'dictionary' belegt ist
                else:
                    st.error(f'fep_file nicht vorhanden: {createKeyName(file.name)}')       # wenn Winkel im 'dictionary' nicht exisitert

            my_bar.progress(createPercentage(checkAmountOfExistingFeps(), len(st.session_state['fep'])), text = f'{checkAmountOfExistingFeps()}/{len(st.session_state['fep'])}')
            if checkAmountOfExistingFeps() == len(st.session_state['fep']):
                st.session_state['progress'].append(3)
        else:
            rewindProgress(2)
            

    # ------- 3 ------- #

    custom_subheader("3. Select your Selectors")
    if(3 in st.session_state['progress']):

        file_upload = st.file_uploader(
                                        "...",
                                        label_visibility = "collapsed",
                                        key = "dataframe_upload",                           
        )

        if file_upload != None:
            ''' '''
        else:
            rewindProgress(3)



    if(checkAmountOfExistingFeps() == len(st.session_state['fep']) and 3 in st.session_state['progress']):
        with tab_main3:
            initialize_custom_Glycan("M5")
            ''' '''


    # st.markdown(createProgressBar(), unsafe_allow_html=True)
    # st.markdown("<div id='progress_start'></div>", unsafe_allow_html=True)



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