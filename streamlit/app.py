import csv
import importlib
from pathlib import Path
import tempfile
import time
import streamlit as st
import plumed # type: ignore
# from glyconformer.lib import Glyconformer, Glycompare
import pandas as pd # type: ignore
import os as os
import sys
from st_clickable_images import clickable_images

sys.path.append("/home/eberl/bachelor_project/GlyCONFORMER/glyconformer/")
from lib import Glyconformer

st.set_page_config(
    page_title = "GlyCONFORMER",
    page_icon = "🧪",
    layout = "centered",
    initial_sidebar_state = "expanded",
    menu_items = {
        'Get Help': 'https://github.com/IsabellGrothaus/GlyCONFORMER',
        'Report a bug': 'https://github.com/IsabellGrothaus/GlyCONFORMER',
    }
)

with open('script.js') as f:
    st.markdown(f"<script>{f.read()}</script>", unsafe_allow_html=True)

with open('styles.css') as f:
    st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)


# -------- INITIALIZATION --------- #

def setSessionStates():

    if 'glycan' not in st.session_state:
        st.session_state['glycan'] = None
    if 'progress' not in st.session_state:
        st.session_state['progress'] = []

    if 'request_state' not in st.session_state:
        st.session_state['request_state'] = False

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

    if 'container_position' not in st.session_state:
        st.session_state['container_position'] = 0
    if 'font_size' not in st.session_state:
        st.session_state['font_size'] = 15
    if 'colors' not in st.session_state:
        st.session_state['colors'] = {'first_color': '#173c4d', 'second_color': '#146b65', 'third_color': '#4e9973'}


setSessionStates()


background_text: str = f"""
    The GlyCONFORMER utility was initially designed as a nomenclature to label glycan conformers distinctively, providing an open-access community resource based on JupyterNotebooks.
    \nDevelopment of this Web application and graphical user interface is maintained by Isabell L. Grothaus, Theo Eberl, Klaus Reuter and Matt Sikora. 
    \nThe GlyCONFORMER WebApp source code as well as JupyterNotebooks to use GlyCONFORMER remotely are available on [GitHub]("https://github.com/IsabellGrothaus/GlyCONFORMER")
    \nPlease support us by citing our paper:
    [Grothaus]("https://pubs.acs.org/doi/full/10.1021/acs.jcim.2c01049") et al. DOI: 10.1021/acs.jcim.2c01049
                    """

caption_text: str = f"""
                Example of conformer string generation. 
                Assignment of a conformer string to the complex N-glycan A4G4S4 (blue: GlcNAc, green: Man, yellow: Gal, purple: Neu5Ac, red:Fuc) 
                drawn with DrawnGlycan 10.1093/glycob/cww115, based on its torsion angle conformations. 
                The various branches of a glycan are ranked according to their type of linkage and calculated free energy profiles of the 
                torsion angles aligned in the corresponding order. Their free energy minimum basins are labeled with respect to the IUPAC 
                nomenclature for torsion angles. For each adopted conformation along a trajectory, 
                the value of each torsion angle in the string is evaluated, assigned to a minimum and correspondingly labeled. 
                The outcome is an ordered string of digits that represents the global conformation of a glycan. 
                """

how_it_works_text: str = f"""
            N-glycans are multi-branched structures, characterized by the specific
            linkages between saccharides monomers. Each glycosidic linkage gives rise
            to at least two torsion angles (${r"\phi"}$ and ${r"\psi"}$), while 1${r"\rarr"}$6 and 2${r"\rarr"}$6 linkages 
            harbor an additional torsion angle ${r"\omega"}$. 
            \nBased on these structural characteristics, we constructed an unambiguous 
            labeling scheme to distinguish different conformers of the same N-glycan. 
            The scheme is also applicable to other glycans,
            independently of their size, number or type of branches and amount of substituents 
            such as fucosylation. Each conformer is identified by a digit string
            of length ${"N_z"}$, equal to the number of torsion angles in the glycan. 
            \nFor N-glycans, the string begins at the free reducing end, consisting of a ${r"\beta"}$ 1 ${r"\rarr"}$ 4
            - linked GlcNAc dimer followed by a mannose residue. For each linkage,
            the linear string reports digits assigned to ${r"\phi"}$, ${r"\psi"}$ and ${r"\omega"}$ (if applicable), in
            this order. In correspondence of a junction (leading e.g. to an ${r"\alpha"}$ 1${r"\rarr"}$6 and
            a ${r"\alpha"}$ 1${r"\rarr"}$3 branch after the first mannose), a string separator is introduced,
            labeled according to the C atom at the branch origin (e.g. ${r"\bold{6-}"}$ for 1${r"\rarr"}$6 linkages). 
            The string continues first along the branch of the higher C atom (6 in
            our case) until reaching the terminal residue, prior to returning to the last
            junction and following the branch of the next-lower C atom (3 in our case).
            Additional modifications like the attachment of fucose residues or bisecting
            GlcNAc residues are included after all other branches are assigned. 
            The separators of primary branches are labelled with bold numbers (${r"\bold{6-}"}$ or ${r"\bold{3-}"}$), for
            clarity. The string digits indicate in which intervals of values the torsion angle lies, 
            following the ${r"\mathrm{IUPAC}"}$ nomenclature for dihedrals 10.1351/goldbook.
            Namely, the digits for ${r"\phi"}$ and ${r"\psi"}$ and the corresponding interval of radian values
            are:

            \n${r"\mathrm{C = [-0.52, +0.52]}"}$
            \n${r"\mathrm{{G_+} = [+0.52, +1.57]}"}$
            \n${r"\mathrm{{A_+} = [+1.57, +2.62]}"}$
            \n${r"\mathrm{T = [+2.62, {\pi}]}"}$
            \n${r"\mathrm{{A_-} = [-2.62, -1.57]}"}$
            \n${r"\mathrm{{G_-} = [-1.57, -0.52]}"}$

            \nThe digits for ${r"\omega"}$ are:
            \n${r"\mathrm{gg = [-2.62, 0]}"}$
            \n${r"\mathrm{gt = [0, 2.62]}"}$
            \n${r"\mathrm{tg = [2.62, {\pi}]}"}$ or ${r"\mathrm{[-{\pi}, 2.62]}"}$.

            \nThe assignment of each torsion angle to a given interval is performed in the
            following way. First, the free-energy profile associated with rotation along the
            torsion angle is calculated from an MD trajectory (most often enhanced sampling MD is 
            necessary in order to achieve converged profiles). 
            \nThe positions of the free-energy minima are then labeled according to the nomenclature
            above. All angles belonging to the same free-energy basin (around a minimum between 
            the two neighboring maxima) are finally labeled equally to the minimum of their basin (Figure ${r"\mathrm{??}"}$). 
            \nAs a last step, each recorded set of torsion angles from a frame of a trajectory is translated into a conformer
            string, built according to the rules above. The scanning of free energy profiles for minima and maxima 
            and the subsequent assignment according to the ${r"\mathrm{IUPAC}"}$ nomenclature are a tedious task if done by hand. 
            Therefore, the workflow was automatized, categorizing free energy profiles according to the scheme in Figure ${r"\mathrm{??}"}$.
            \nThe package can also read in a feature matrix ${r"\bold{X}"}$ of
            shape ${r"{n_{samples}} \times {n^{features}}"}$, where the features are equal to all torsion angles of
            the analyzed glycan, and convert the torsion angle values of each frame into
            the respective conformer strings. This data set can be subsequently used to
            construct a histogram with the individual conformers selected as bins, yielding a conformer distribution. 
            In order to assess statistical features of this distribution, block averaging can be performed, separating the data set into
            evenly distributed blocks.
            \nThe average of all blocks ${r"\bar{X} = \frac{1}{N}\sum_{j-1}^N{X_j}"}$ is calculated over ${"N = 10"}$, where ${"X_j"}$ is the average 
            calculated within each ${"j"}$th block. Error bars are calculated as standard deviations of the sampling distribution 
            (standard error of the mean): ${r"std(X) = \sqrt{\frac{var(X)}{N}}"}$ with the variance of the 
            sampling distribution ${r"var(X) = (\frac{N}{N-1})[\frac{1}{N}\sum_{j-1}^N {X_{j}^2} - (\frac{1}{N}\sum_{j-1}^N {X_j})]"}$.
            \nThe outlined analysis represents a fundamental basis for the assessment of glycan simulations, as the conformer string 
            incorporates all degrees of freedom that are necessary to describe the global glycan conformation.
"""



# -------- MAIN PAGE --------- #

container = st.container()


def buildHUD():
    if st.session_state['request_state']:

        Glycan = st.session_state['glycan']

        with container:
            st.header(f"Data for: {Glycan.glycantype}")

            with st.expander("Configuration", expanded = False):
                st.session_state['font_size'] = st.select_slider('Font-Size', options = list(range(10, 16)), value = 12)
                
                st.write("Pick a color:")
                subcol1, subcol2, subcol3 = st.columns(3, gap = "medium")
                with subcol1:
                    st.session_state['colors']['first_color'] = st.color_picker("Pick A Color", "#173c4d", label_visibility = "collapsed")
                with subcol2:
                    st.session_state['colors']['second_color'] = st.color_picker("Pick A Color", "#146b65", label_visibility = "collapsed")
                with subcol3:
                    st.session_state['colors']['third_color'] = st.color_picker("Pick A Color", "#4e9973", label_visibility = "collapsed")

            tab1, tab2, tab3, tab4 = st.tabs(["pca", "fep", "distribution", "average"])

            with tab1:   
                st.pyplot(Glycan.pca(fontsize = st.session_state['font_size'], color = [st.session_state['colors']['first_color'],st.session_state['colors']['second_color'],st.session_state['colors']['third_color'],"#a7c09f","#dfa790","#c76156","#9a2b4c","#600b4a"]))  
                st.pyplot(Glycan.pca_fep(fontsize = st.session_state['font_size']))

            with tab2:
                st.pyplot(Glycan.robust_validate_fep(fontsize = st.session_state['font_size']))

            with tab3:
                st.pyplot(Glycan.distribution(fontsize = st.session_state['font_size'], colors = [st.session_state['colors']['first_color'],st.session_state['colors']['second_color'],st.session_state['colors']['third_color'],"#a7c09f","#dfa790","#c76156","#9a2b4c","#600b4a"]))
                
            with tab4:
                st.pyplot(Glycan.moving_average(simulation_length = 500, window = 12500, fontsize = st.session_state['font_size'], color = [st.session_state['colors']['first_color'],st.session_state['colors']['second_color'],st.session_state['colors']['third_color'],"#a7c09f","#dfa790","#c76156","#9a2b4c","#600b4a"]))
                st.pyplot(Glycan.cumulative_average(simulation_length = 500, fontsize = st.session_state['font_size'], color = [st.session_state['colors']['first_color'],st.session_state['colors']['second_color'],st.session_state['colors']['third_color'],"#a7c09f","#dfa790","#c76156","#9a2b4c","#600b4a"]))


    else:
        
        with container:

            # -------- Header --------- #

            st.markdown("""
                <div class='headline-container'>
                    <h1>GlyCONFORMER</h1>
                </div>
                
            """, unsafe_allow_html=True)

            # -------- Body --------- #

            main_col1, main_col2, main_col3 = st.columns(3)

            with main_col1:
                if st.button("back >", use_container_width=True):
                    st.session_state['container_position'] = 0
                with st.container():
                    st.subheader("How it works")
                    st.image("images/workflow.png", caption = caption_text)
                    st.write(how_it_works_text)
                

            with main_col2:
                left, right= st.columns(2)

                with left:
                    if st.button("< how it works", use_container_width=True):
                        st.session_state['container_position'] = 110

                with right:
                    if st.button("behind the project >", use_container_width=True):
                        st.session_state['container_position'] = -110


            with main_col3:
                if st.button("< back", use_container_width=True):
                    st.session_state['container_position'] = 0
                st.subheader("Background")
                st.write(background_text)


            # -------- Footer --------- #

            with st.container():
                col1, col2, col3, col4 = st.columns(4)

                with col1:
                    st.image("images/logo_university_bremen.png", caption=None, use_column_width=True)
                with col2:
                    st.image("images/logo_max_planck_biophysics.png", caption=None, use_column_width=True)
                with col3:
                    st.image("images/logo_max_planck_computing.png", caption=None, use_column_width=True)
                with col4:
                    st.image("images/logo_university_krakow.png", caption=None, use_column_width=True)

            change_slide()


# -------- local FUNCTIONS--------- #

def loadLocalTypes(): 

    local_glycans = []
    for i in os.listdir("../LIBRARY_GLYCANS"):
        if("__" not in i):
            local_glycans.append(i)

    return local_glycans

def loadLocalGlycans(glycan_type: str): 

    glycans: dict = {}

    for i in os.listdir(f"../LIBRARY_GLYCANS/{glycan_type}"):
        if("__" not in i):
            glycans[i] = {}
            for e in os.listdir(f"../LIBRARY_GLYCANS/{glycan_type}/{i}"):
                glycans[i][e] = {}

                glycans[i][e]['name'] = e

                for file in os.listdir(f"../LIBRARY_GLYCANS/{glycan_type}/{i}/{e}"):
                    
                    if file.endswith(".png"):
                        glycans[i][e]['picture'] = f"../LIBRARY_GLYCANS/{glycan_type}/{i}/{e}/{file}"
                        break
                    else:
                        glycans[i][e]['picture'] = f"../glycan_test.png"

    return glycans

def readAnglesData(path, file):

    with open(os.path.join(Path(os.getcwd()).parent, path, file), 'r') as f:        # benötigt den absoluten Pfad
        feature = f.read()

    feature = feature.split()

    return feature

def readSeparatorData(path, file):

    dictionary = {'separator_index': [], 'separator': []}

    if path is not None:            # liest aus einer lokalen Datei 
        with open(os.path.join(Path(os.getcwd()).parent, path, file), 'r') as f:        # benötigt den absoluten Pfad
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
def initialize_local_Glycan(type: str, subject: str, glycan: str):
    st.session_state['request_state'] = True

    Glycan = Glyconformer(
        glycantype =        glycan,
        inputfile =         f"../LIBRARY_GLYCANS/{type}/{subject}/{glycan}/colvar.dat",
        fepdir =            f"../LIBRARY_GLYCANS/{type}/{subject}/{glycan}", 

        angles =            readAnglesData(f"LIBRARY_GLYCANS/{type}/{subject}/{glycan}", "angles.dat"), 
        omega_angles =      readAnglesData(f"LIBRARY_GLYCANS/{type}/{subject}/{glycan}", "omega_angles.dat"), 

        separator_index =   readSeparatorData(f"LIBRARY_GLYCANS/{type}/{subject}/{glycan}", "separator.dat")['separator_index'], 
        separator =         readSeparatorData(f"LIBRARY_GLYCANS/{type}/{subject}/{glycan}", "separator.dat")['separator'],

        fep_files =         {},
        order_max =         5,
        order_min =         5, 
        weights =           None,

        colvar =            None,
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


# -------- advanced FUNCTIONS--------- #

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

def change_slide():

    st.write(
        f"""
        <style>
            .element-container:has(div > div > div.headline-container) + div {{
                transform: translateX({st.session_state['container_position']}%);
            }}
        </style>
        """
    , unsafe_allow_html=True)

def stream_data(text: str, id: str):
        
    yield st.session_state['previous_text'][id]

    for i in range(st.session_state['index'][id], len(text)):
        yield text[i]
        st.session_state['index'][id] += 1
        st.session_state['previous_text'][id] = st.session_state['previous_text'][id] + text[i]
        time.sleep(0.004)
    


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

    help_text: dict = {
        "1": "Upload your dataframe with the raw data of the torsion angles. The order of the torsion angles is essential for later analyses.", 
        "2": "The required upload is determined based on the identified angles in the data frame.", 
        "3": "Upload your file with the separators. The file must contain the indexed position and the character of the separator",
        "4": "Determine additional properties. These are set by default and are optional for initialising the data" }

    if(int(text[0]) in st.session_state['progress']):
        st.markdown(f"<div class='custom-subheader-active'></div>", unsafe_allow_html=True)
        st.subheader(text, help=help_text[text[0]])
    else:
        st.markdown(f"<div class='custom-subheader-inactive'></div>", unsafe_allow_html=True)
        st.subheader(text, help=help_text[text[0]])

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
        slider_value = st.select_slider('Size (%)', options = list(range(10, 101, 10)), value = 100)
        st.session_state['length'] = determineLength(slider_value)

    else:
        rewindProgress(3)


    # ------- BUTTON ------- #

    if st.button("start", use_container_width = True, type = "primary" ,disabled = not is_completed):
        st.session_state['glycan'] = initialize_custom_Glycan(st.session_state['dataframe_upload'].name)


# @st.fragment
def buildFragment(picture :str, type: str, subject: str, name: str):

    st.image(picture, caption = None)

    if st.button(f"{name}", use_container_width = True, key = f"{name}"):
        st.session_state['glycan'] = initialize_local_Glycan(type, subject, name)

def checkSelect():
    st.subheader("1. Select a dataset")

    selected_Type = st.selectbox(
                                    "...",
                                    loadLocalTypes(), 
                                    index = None,                              # initialisiert eine leere Auswahlbox                                           
                                    label_visibility = "collapsed",
                                    key = "glycans_select",
                                    placeholder = "choose one of the following glycans",
    )

    if selected_Type is not None:                    

        local_glycans: dict = loadLocalGlycans(selected_Type)

        for key, value in local_glycans.items():
            st.write(key)
            data = []

            for glycan in value.items():
                data.append(glycan)

            data_groups = [data[i:i+3] for i in range(0, len(data), 3)]


            for group in data_groups:
                col1, col2, col3 = st.columns(3)
                for i in range(0, len(group)):
                    # print(f"index: {i}, glycan: {group[i][1]}")
                    if i == 0:
                        with col1:                       
                            buildFragment(group[i][1]['picture'], selected_Type, key, group[i][1]['name'])

                    elif i == 1:
                        with col2:                           
                            buildFragment(group[i][1]['picture'], selected_Type, key, group[i][1]['name'])

                    else:
                        with col3:                         
                            buildFragment(group[i][1]['picture'], selected_Type, key, group[i][1]['name'])



# -------- SIDEBAR --------- #

with st.sidebar:

    tab_sidebar1, tab_sidebar2 = st.tabs(["User Input", "Preselection"])

    with tab_sidebar1:
        checkProgress()
            

    with tab_sidebar2:
        checkSelect()
        
        
    buildHUD()


# -------- ### --------- #

