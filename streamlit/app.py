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
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.2.0/css/all.min.css">
""", unsafe_allow_html=True)


def setSessionStates():

    if 'glycan' not in st.session_state:
        st.session_state['glycan'] = None
    if 'progress' not in st.session_state:
        st.session_state['progress'] = []

    if 'request_state' not in st.session_state:
        st.session_state['request_state'] = False
    if 'request_access' not in st.session_state:
        st.session_state['request_access'] = False


    if 'dataframe' not in st.session_state:
        st.session_state['dataframe'] = None
    if 'angles' not in st.session_state:
        st.session_state['angles'] = None
    if 'separators' not in st.session_state:
        st.session_state['separators'] = None
    if 'fep' not in st.session_state:
        st.session_state['fep'] = {}
    if 'length' not in st.session_state:
        st.session_state['length'] = None

setSessionStates()


# -------- MAIN PAGE --------- #

container = st.container()

def buildHUD():
    if st.session_state['request_state']:

        Glycan = st.session_state['glycan']
        
        with container:
            st.header(f"Data for: {Glycan.glycantype}")
            tab1, tab2 = st.tabs(["test1", "test2"])

            with tab1:
                st.write(Glycan.colvar)
                st.write(st.session_state['separators'])
                st.write(Glycan.minima["phi1_2"][0])
                st.write(Glycan.separator_index)
                st.write(Glycan.separator)
                    
                st.pyplot(Glycan.robust_validate_fep())
            with tab2:
                st.pyplot(Glycan.pca(biplot = True))
                st.write(Glycan.distribution())
                st.pyplot(Glycan.pca_fep())
                st.pyplot(Glycan.moving_average(simulation_length = 500, window = 12500))
                st.pyplot(Glycan.cumulative_average(simulation_length = 500))

    else:
        ''' '''
        with container:
            st.markdown("""
                <div class='headline-container'>
                    <h1>GlyCONFORMER</h1>
                </div>
            """, unsafe_allow_html=True)
           
        # with container:
        #     st.markdown("<i class='fa-solid fa-solid fa-flask-vial'></i>", unsafe_allow_html=True)


# -------- local FUNCTIONS--------- #

def loadLocalGlycans(): 

    local_glycans = []
    for i in os.listdir("../LIBRARY_GLYCANS/"):
        if("__" not in i):
            local_glycans.append(i)

    return local_glycans

def readAnglesData(path, file):
    with importlib.resources.files(path).joinpath(file).open() as f:
        feature = f.read()
    feature = feature.split()

    return feature

def readSeparatorData(path, file):

    dictionary = {'separator_index': [], 'separator': []}

    if path is not None:            # liest aus einer lokalen Datei 
        with importlib.resources.files(path).joinpath(file).open('r') as f:         # öffnet Datei im Lesemodus 'r'
            reader = csv.reader(f, delimiter=' ')      

            for row in reader:
                dictionary['separator_index'].append(int(row[0]))
                dictionary['separator'].append(str(row[1]))

    else:                           # liest aus einer hochgeladenen Datei 
        reader = pd.read_csv(file, sep=" ", header=None)           

        dictionary['separator_index'] = reader[0].tolist()
        dictionary['separator'] = reader[1].tolist()
        
    return dictionary

# @st.cache_data                              # wenn der Glykan-String bereits abgerufen wurde, bezieht er die Daten aus dem Cache
def initialize_local_Glycan(glycantype: str):
    st.session_state['request_state'] = True

    Glycan = Glyconformer(
        glycantype =        glycantype,
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

    return Glycan

def initialize_custom_Glycan(glycantype: str):
    st.session_state['request_state'] = True

    Glycan = Glyconformer(
        glycantype =        glycantype,
        inputfile =         None,
        fepdir =            None,

        angles =            st.session_state['angles']['angles'].copy(), 
        omega_angles =      st.session_state['angles']['omega_angles'].copy(), 
        separator_index =   st.session_state['separators']['separator_index'].copy(), 
        separator =         st.session_state['separators']['separator'].copy(), 
        fep_files =         st.session_state['fep'].copy(),
        order_max =         5,
        order_min =         5, 
        weights =           None,

        colvar =            st.session_state['dataframe'].copy(), 
        length =            st.session_state['length']
    )

    return Glycan



# -------- local FUNCTIONS--------- #

def on_change():

    st.session_state['request_access'] = True

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

def determineLength(value: int):

    return int(len(st.session_state['dataframe']) * (value/100))

def rewindProgress(number: int):
    st.session_state['progress'].clear()  

    for i in range(1, number + 1):
        st.session_state['progress'].append(i)

    # print(f"----------------------------- {st.session_state['progress']}")

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
    for angle in st.session_state['angles']['angles']:
        st.session_state['fep'][angle] = []

def extractAngles(colvar):

    angles = {'angles': [], 'omega_angles': []}

    for i in colvar.columns:
        if i != "time" and "pseudo" not in i:
            angles['angles'].append(i)

            if "omega" in i:
                angles['omega_angles'].append(i)

    return angles

def readDataframe(uploaded_file):
    
    # In Bytestring konvertieren
    bytes_data: bytes = uploaded_file.getvalue()

    # In temporäre Datei schreiben
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        temp_file.write(bytes_data)
        temp_file_path = temp_file.name

    colvar = plumed.read_as_pandas(temp_file_path)
    st.session_state['angles'] = extractAngles(colvar)                               # füllt 'angles' und 'omega_angles' in st.session_state['angles']

    colvar = colvar[st.session_state['angles']['angles']]

    # colvar2 = colvar[['phi1_2', 'psi1_2', 'phi2_3', 'psi2_3', 'phi3_5', 'psi3_5', 'omega3_5', 'phi5_7', 'psi5_7', 'omega5_7', 'phi5_6', 'psi5_6', 'phi3_4', 'psi3_4']]


    # Temporäre Datei löschen
    os.remove(temp_file_path)

    return colvar

def checkProgress():
    is_completed: bool = False

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

    custom_subheader("3. Select your Separators")

    if(3 in st.session_state['progress']):

        my_bar = st.progress(0, text = "0/1")

        file_upload = st.file_uploader(
                                        "...",
                                        label_visibility = "collapsed",
                                        key = "selector_upload",                           
        )

        if file_upload is not None:
            st.session_state['progress'].append(4)
            st.session_state['separators'] = readSeparatorData(None, file_upload)

            my_bar.progress(createPercentage(1, 1), text = "1/1")

            is_completed = True
        else:
            rewindProgress(3)


    # ------- 4 ------- #

    custom_subheader("4. Select your Specifics (optional)")

    if(4 in st.session_state['progress']):
        slider_value = st.select_slider('Size', options = list(range(10, 101, 10)), value = 100)
        st.session_state['length'] = determineLength(slider_value)

    else:
        rewindProgress(3)


    # ------- BUTTON ------- #

    if st.button("start", use_container_width = True, disabled = not is_completed):
        st.session_state['glycan'] = initialize_custom_Glycan(st.session_state['dataframe_upload'].name)

def checkSelect():
    st.subheader("1. Select a dataset")

    selected_glycan = st.selectbox(
                                    "...",
                                    loadLocalGlycans(), 
                                    index = None,                              # initialisiert eine leere Auswahlbox                                           
                                    label_visibility = "collapsed",
                                    key = "glycans_select",
                                    placeholder = "choose one of the following glycans",
                                    on_change = on_change
    )

    if selected_glycan is not None and st.session_state['request_access']:               # falls Eintrag Glycan ist UND gewählt wurde -> Initialisierung ist möglich
        st.session_state['glycan'] = initialize_local_Glycan(selected_glycan)
        st.session_state['request_access'] = False


# -------- SIDEBAR --------- #

with st.sidebar:

    tab_sidebar1, tab_sidebar2 = st.tabs(["User Input", "Preselection"])

    with tab_sidebar1:
        checkProgress()
            

    with tab_sidebar2:
        checkSelect()
        
        
    buildHUD()


# -------- ### --------- #

with open('script.js') as f:
    st.markdown(f"<script>{f.read()}</script>", unsafe_allow_html=True)

with open('styles.css') as f:
    st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)