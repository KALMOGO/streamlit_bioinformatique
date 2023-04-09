
######################## Bloc importation des librairies ####################

import streamlit as st # FrameWork Stream lit
from PIL import Image # Lecture d'image 
import subprocess      # execution des executable
from contextlib import contextmanager, redirect_stdout # gestioonaire de context
from ete3 import Tree  # visualisation de l'abre phylogenique

    #----------------Algorithme Biopython--------------
    
from Bio import pairwise2
from Bio import SeqIO
from io import StringIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline


    #----------------Algorithme WPGMA--------------------
from wpgma import * 


############################ Bloc Menu ############################################

st.set_page_config(layout="wide")
st.markdown("""
<style>
.css-1544g2n ,.e1fqkh3o4{
    padding-top:0 !important;
}

.css-1kyxreq ,.etr89bj2{
    justify-content:center;
    align-items:center;
}

.css-z5fcl4 {
    padding:10px !important
}
</style>
""", unsafe_allow_html=True)



col1 = st.sidebar
col2, col3 = st.columns((2,1))

st.sidebar.title('Bioinformatic Web application')
image_logo = Image.open('./images/logo.jpg')

st.sidebar.image(image_logo, width = 200)

st.sidebar.header("Effectuer Votre choix")
option = st.sidebar.radio(
    'Choose an option',
    ('Pairs Alignments', 'Multiple Alignments')
)


if option == 'Pairs Alignments':
    multiple_align_option=None
    st.sidebar.subheader('Alignement de séquences par paires')
    pair_align_option = st.sidebar.selectbox(
        'Choose an option',
        ('Local Alignment', 'Global Alignment', 'Global Alignment with Gap Model')
    )
    
    col2.title('Alignement de séquences par paire')
    image_logo = Image.open('./images/alignlogo.png')
    col3.image(image_logo, width=400)
    

else:
    pair_align_option=None
    st.sidebar.subheader(' Alignement de séquences Multiple & phylogenitique')
    multiple_align_option = st.sidebar.selectbox(
        'Choose an option',
        ('Exact method', 'Clustal method', 'Phylogenetic')
    )
    
    col2.title('Alignement de séquences Multiple & phylogenitique')
    image_logo = Image.open('./images/alignlogo.png')
    col3.image(image_logo, width=400)



# ###########################Page Principale ###########################

    ############# section 1: Pairs Alignment ############
    
            ###### Local Alignment #######
            
if pair_align_option == 'Local Alignment':
    
        # possibilité de rensigner la valeur du gap, de la correspondance ou de la non corrrepondance
    col1.subheader("Paramètres de l'alignement: vide par defaut")
    col_side_1, col_side_2, col_side_3 = st.sidebar.columns((1,1,1))
    match = col_side_1.text_input('match')
    mismatch = col_side_2.text_input('mismatch')
    gap  = col_side_3.text_input('Gap')
    col2.subheader('Local Alignment')
    
    col2.markdown("""
            Cette section permet d'effectuer un alignement local de deux sequence d'ADN  avec la librairie **Biopython de python**! \n
            Deux possibilités sont offertes: \n
            1. Renseigner manullement les sequences
            2. Charger les fichiers des séquences
    """)
    
        # -----------------------------------------
                # Renseigner manullement les sequences
        
    ### fonction de validation de la sequences saisie
    def validate_sequence(sequence):
        if sequence=="" :
            return False
        else:
            isvalid = True
            for char in list(sequence):
                if char not in ["A", "T", "G", "C"]:
                    isvalid =  False
                    break
                
        return isvalid
    
    ## Creation des champs de saisie des sequences
    expander_bar = col2.expander("**Renseigner manullement les sequences**")
    expander_bar.subheader("Only aplhabet: **{A G C T }** is supported & No space")
    expander_bar.markdown("""***Clicker sur Clt+Enter pour valider les séquences***""")
    
    sequence1 = expander_bar.text_area("Sequence1", height=100)
    sequence1 = sequence1.upper() # majuscule
    is_valid1 = validate_sequence(sequence1)
    # Validation seq1
    if not(is_valid1) :
        expander_bar.markdown("""***Invalid Sequence !!!***""")
        
    sequence2 = expander_bar.text_area("Sequence2", height=100)
    sequence2 = sequence2.upper() # seq en majuscule
    is_valid2 = validate_sequence(sequence2)
    
    # Validation seq2
    if not(is_valid2) :
        expander_bar.markdown("""***Invalid Sequence !!!***""")

    ############## Alignement local ###########################
    st.header('RESULT')
    
    if match == "" and mismatch=="" and gap=="" and is_valid1 == True and is_valid2==True:
        st.markdown(""" ***Alignement de séquence*** """)
        st.markdown(""" ***Avec paramètres:*** """)
        params = {'match':'default', 'mismatch':'default', 'gap':'default'}
        st.write(params)
        st.write("SeqA : ",sequence1.upper())
        st.write("SeqB : ", sequence2.upper())
        
        # localxx - matches score 1, mismatches 0 and no gap penalty.
        alignments = pairwise2.align.localxx(sequence1, sequence2) 
        for alignment in alignments:
            st.write(alignment)

    if match != "" and mismatch !="" and gap=="" and is_valid1 == True and is_valid2==True:
        match= int(match)
        mismatch=int(mismatch)
        
        st.markdown(""" ***Alignement de séquence*** """)
        st.markdown(""" ***Avec paramètres:*** """)
        params = {'match':match, 'mismatch':mismatch, 'gap':'default'}
        
        st.write(params )
        st.write("SeqA : ",sequence1.upper())
        st.write("SeqB : ", sequence2.upper())
        
        alignments = pairwise2.align.localmx(sequence1, sequence2, match=match, mismatch=mismatch) 
        for alignment in alignments:
            st.write(alignment)

    if match != "" and mismatch !="" and gap!="" and is_valid1 == True and is_valid2==True:
        match= int(match)
        mismatch=int(mismatch)
        gap = int(gap)
        
        st.markdown(""" ***Alignement de séquence*** """)
        st.markdown(""" ***Avec paramètres:*** """)
        params = {'match':match, 'mismatch':mismatch, 'gap':gap}
        
        st.write(params )
        st.write("SeqA : ",sequence1.upper())
        st.write("SeqB : ", sequence2.upper())
        
        # nous devons effectuer un alignement local une valeur gap. Nous fixons donc tous les type de gap avec les meme valeur
        # il ne s'agit pas d'alignement avec un gag affine 
        alignments = pairwise2.align.localms(sequence1, sequence2, match=match, mismatch=mismatch, open=gap, extend=gap) 
        for alignment in alignments:
            st.write(alignment)
            
    if (match == "" and mismatch !="") or (match !="" and mismatch =="")and is_valid1 == True and is_valid2==True:
        st.markdown(""" ***Erreur: match et mismatch doivent renseigner ensemble ou pas du tout!!*** """)


        # -----------------------------------------
                # Renseigner des fichiers
        # ----------------------------------------
    if is_valid1 == False:
        expander_bar = col2.expander("**Charger les fichiers des séquences**")
        expander_bar.subheader("Only FASTA file is supported")
        
        st.markdown(""" ***Only the default match, mismatch and gap is used*** """)
        
        sequence_file1 = expander_bar.file_uploader("***Sequence 1***", type=["fasta"])
        if sequence_file1 is not None:
            st.write("...chargement sequence 1: OK")
        else:
            st.write('Awaiting fasta file to be uploaded for sequence 1. No file Upload')

        sequence_file2 = expander_bar.file_uploader("***Sequence 2***", type=["fasta"])
        if sequence_file2 is not None:
            st.write("...chargement sequence 2: OK")
        else:
            st.write('Awaiting fasta file to be uploaded for sequence 2. No file Upload')

        # Alignement local
        if sequence_file2 is not None and sequence_file1 is not None:
            stringio1 = StringIO(sequence_file1.getvalue().decode("utf-8"))
            stringio2 = StringIO(sequence_file2.getvalue().decode("utf-8"))
            
            for record in SeqIO.parse(stringio1, 'fasta'):
                sequence1 = str(record.seq)
                
            for record in SeqIO.parse(stringio2, 'fasta'):
                sequence2 = str(record.seq)  
                
                # localxx - matches score 1, mismatches 0 and no gap penalty.
            alignments = pairwise2.align.localxx(sequence1, sequence2) 
            for alignment in alignments:
                st.write(alignment)





        ###### Gloab Alignment #########""

elif pair_align_option == 'Global Alignment':
    
        # possibilité de rensigner la valeur du gap, de la correspondance ou de la non corrrepondance
    col1.subheader("Paramètres de l'alignement: vide par defaut")
    col_side_1, col_side_2, col_side_3 = st.sidebar.columns((1,1,1))
    match = col_side_1.text_input('match')
    mismatch = col_side_2.text_input('mismatch')
    gap  = col_side_3.text_input('Gap')
    col2.subheader('Global Alignment')
    
    col2.markdown("""
            Cette section permet d'effectuer un alignement Global de deux sequence d'ADN  avec la librairie **Biopython de python**! \n
            Deux possibilités sont offertes: \n
            1. Renseigner manullement les sequences
            2. Charger les fichiers des séquences
    """)
    
        # -----------------------------------------
                # Renseigner manullement les sequences
        # ----------------------------------------
        
    # fonction de validation de la sequences saisie
    def validate_sequence(sequence):
        if sequence=="" :
            return False
        else:
            isvalid = True
            for char in list(sequence):
                if char not in ["A", "T", "G", "C"]:
                    isvalid =  False
                    break
                
        return isvalid
    
    # Creation des champs de saisie des sequences
    expander_bar = col2.expander("**Renseigner manullement les sequences**")
    expander_bar.subheader("Only aplhabet: **{A G C T }** is supported & No space")
    expander_bar.markdown("""***Clicker sur Clt+Enter pour valider les séquences***""")
    
    sequence1 = expander_bar.text_area("Sequence1", height=100)
    sequence1 = sequence1.upper() # majuscule
    is_valid1 = validate_sequence(sequence1)
    # Validation seq1
    if not(is_valid1) :
        expander_bar.markdown("""***Invalid Sequence !!!***""")
        
    sequence2 = expander_bar.text_area("Sequence2", height=100)
    sequence2 = sequence2.upper() # seq en majuscule
    is_valid2 = validate_sequence(sequence2)
    
    # Validation seq2
    if not(is_valid2) :
        expander_bar.markdown("""***Invalid Sequence !!!***""")

    ############## Alignement Global ###########################
    st.header('RESULT')
    
    if match == "" and mismatch=="" and gap=="" and is_valid1 == True and is_valid2==True:
        st.markdown(""" ***Alignement de séquence*** """)
        st.markdown(""" ***Avec paramètres:*** """)
        params = {'match':'default', 'mismatch':'default', 'gap':'default'}
        st.write(params)
        st.write("SeqA : ",sequence1.upper())
        st.write("SeqB : ", sequence2.upper())
        
        # localxx - matches score 1, mismatches 0 and no gap penalty.
        alignments = pairwise2.align.globalxx(sequence1, sequence2) 
        for alignment in alignments:
            st.write(alignment)

    if match != "" and mismatch !="" and gap=="" and is_valid1 == True and is_valid2==True:
        match= int(match)
        mismatch=int(mismatch)
        
        st.markdown(""" ***Alignement de séquence*** """)
        st.markdown(""" ***Avec paramètres:*** """)
        params = {'match':match, 'mismatch':mismatch, 'gap':'default'}
        
        st.write(params )
        st.write("SeqA : ",sequence1.upper())
        st.write("SeqB : ", sequence2.upper())
        
        alignments = pairwise2.align.globalmx(sequence1, sequence2, match=match, mismatch=mismatch) 
        for alignment in alignments:
            st.write(alignment)

    if match != "" and mismatch !="" and gap!="" and is_valid1 == True and is_valid2==True:
        match= int(match)
        mismatch=int(mismatch)
        gap = int(gap)
        
        st.markdown(""" ***Alignement de séquence*** """)
        st.markdown(""" ***Avec paramètres:*** """)
        params = {'match':match, 'mismatch':mismatch, 'gap':gap}
        
        st.write(params )
        st.write("SeqA : ",sequence1.upper())
        st.write("SeqB : ", sequence2.upper())
        
        # nous devons effectuer un alignement global une valeur gap. Nous fixons donc tous les type de gap avec les meme valeur
        # il ne s'agit pas d'alignement avec un gag affine 
        alignments = pairwise2.align.globalms(sequence1, sequence2, match=match, mismatch=mismatch, open=gap, extend=gap) 
        for alignment in alignments:
            st.write(alignment)
            
    if (match == "" and mismatch !="") or (match !="" and mismatch =="")and is_valid1 == True and is_valid2==True:
        st.markdown(""" ***Erreur: match et mismatch doivent renseigner ensemble ou pas du tout!!*** """)


        # -----------------------------------------
                # Renseigner des fichiers
        # ----------------------------------------
    if is_valid1 == False:
        expander_bar = col2.expander("**Charger les fichiers des séquences**")
        expander_bar.subheader("Only FASTA file is supported")
        
        st.markdown(""" ***Only the default match, mismatch and gap is used*** """)
        
        sequence_file1 = expander_bar.file_uploader("***Sequence 1***", type=["fasta"])
        if sequence_file1 is not None:
            st.write("...chargement sequence 1: OK")
        else:
            st.write('Awaiting fasta file to be uploaded for sequence 1. No file Upload')

        sequence_file2 = expander_bar.file_uploader("***Sequence 2***", type=["fasta"])
        if sequence_file2 is not None:
            st.write("...chargement sequence 2: OK")
        else:
            st.write('Awaiting fasta file to be uploaded for sequence 2. No file Upload')

        # Alignement local
        if sequence_file2 is not None and sequence_file1 is not None:
            stringio1 = StringIO(sequence_file1.getvalue().decode("utf-8"))
            stringio2 = StringIO(sequence_file2.getvalue().decode("utf-8"))
            
            for record in SeqIO.parse(stringio1, 'fasta'):
                sequence1 = str(record.seq)
                
            for record in SeqIO.parse(stringio2, 'fasta'):
                sequence2 = str(record.seq)  
                
                # localxx - matches score 1, mismatches 0 and no gap penalty.
            alignments = pairwise2.align.globalxx(sequence1, sequence2) 
            for alignment in alignments:
                st.write(alignment)



        ###### Alignement Global avec un model de gap #########""

elif pair_align_option == 'Global Alignment with Gap Model':
        # possibilité de rensigner la valeur du gap, de la correspondance ou de la non corrrepondance
    col1.subheader("Paramètres de l'alignement avec gag affine")
    
    col_side_1, col_side_2 = st.sidebar.columns((1,1))
    open_gap = col_side_1.text_input('Open gap')
    extend_gap = col_side_2.text_input('Extend gap')
    
    col_side_3, col_side_4 = st.sidebar.columns((1,1))
    match = col_side_3.text_input('match')
    mismatch = col_side_4.text_input('mismatch')
    
    col2.subheader('Global Alignment with Gap Model')
    
    col2.markdown("""
            Cette section permet d'effectuer un alignement Global avec un modèle de gap lineaire ou affine
            de deux sequence d'ADN  avec la librairie **Biopython de python**! \n
            Deux possibilités sont offertes: \n
            1. Renseigner manullement les sequences
            2. Charger les fichiers des séquences
            \n Par defaut un modèle de gap lineaire est effectué
            \n le modèle de gag affine s'effectue en reseignant 
            les valeurs des gap d'ouverture et d'extension""")
    
        # -----------------------------------------
                # Renseigner manullement les sequences
        # ----------------------------------------
        
    # fonction de validation de la sequences saisie
    def validate_sequence(sequence):
        if sequence=="" :
            return False
        else:
            isvalid = True
            for char in list(sequence):
                if char not in ["A", "T", "G", "C"]:
                    isvalid =  False
                    break
                
        return isvalid
    
    # Creation des champs de saisie des sequences
    expander_bar = col2.expander("**Renseigner manullement les sequences**")
    expander_bar.subheader("Only aplhabet: **{A G C T }** is supported & No space")
    expander_bar.markdown("""***Clicker sur Clt+Enter pour valider les séquences***""")
    
    sequence1 = expander_bar.text_area("Sequence1", height=100)
    sequence1 = sequence1.upper() # majuscule
    is_valid1 = validate_sequence(sequence1)
    # Validation seq1
    if not(is_valid1) :
        expander_bar.markdown("""***Invalid Sequence !!!***""")
        
    sequence2 = expander_bar.text_area("Sequence2", height=100)
    sequence2 = sequence2.upper() # seq en majuscule
    is_valid2 = validate_sequence(sequence2)
    
    # Validation seq2
    if not(is_valid2) :
        expander_bar.markdown("""***Invalid Sequence !!!***""")

    ############## Alignement Global ###########################
    st.header('RESULT')

    if open_gap != "" and extend_gap !="" and is_valid1 == True and is_valid2==True and match!="" and mismatch!="":
        match= int(match)
        mismatch=int(mismatch)
        open_gap = int(open_gap)
        extend_gap=int(extend_gap)
        
        st.markdown(""" ***Alignement de séquence: gap affine*** """)
        st.markdown(""" ***Avec paramètres:*** """)
        params = { 'mismatch':mismatch,'match':match,'open_gap':open_gap, 'extend_gap':extend_gap}
        
        st.write(params )
        st.write("SeqA : ",sequence1.upper())
        st.write("SeqB : ", sequence2.upper())
        
        # nous devons effectuer un alignement global une valeur gap. Nous fixons donc tous les type de gap avec les meme valeur
        alignments = pairwise2.align.globalms(sequence1, sequence2, match=match, mismatch=mismatch, open=open_gap, extend=extend_gap) 
        for alignment in alignments:
            st.write(alignment)

    if open_gap == "" and extend_gap =="" and is_valid1 == True and is_valid2==True and match=="" and mismatch=="":
        
        st.markdown(""" ***Alignement de séquence: gap lineaire*** """)
        st.markdown(""" ***Avec paramètres:*** """)

        st.write("SeqA : ",sequence1.upper())
        st.write("SeqB : ", sequence2.upper())
        
        # Gap lineaire : valeur par defaut
        alignments = pairwise2.align.globalxx(sequence1, sequence2) 
        for alignment in alignments:
            st.write(alignment)

    if open_gap != "" and extend_gap !="" and is_valid1 == True and is_valid2==True and match=="" and mismatch=="":
        open_gap = int(open_gap)
        extend_gap=int(extend_gap)
        
        st.markdown(""" ***Alignement de séquence: gap affine*** """)
        st.markdown(""" ***Avec paramètres:*** """)
        params = { 'mismatch':'default','match':'default','open_gap':open_gap, 'extend_gap':extend_gap}
        
        st.write(params )
        st.write("SeqA : ",sequence1.upper())
        st.write("SeqB : ", sequence2.upper())
        
        # Gap lineaire : valeur par defaut
        alignments = pairwise2.align.globalxs(sequence1, sequence2, open=open_gap, extend=extend_gap) 
        for alignment in alignments:
            st.write(alignment)

    if open_gap == "" and extend_gap =="" and is_valid1 == True and is_valid2==True and match!="" and mismatch!="":
        
        st.markdown(""" ***Alignement de séquence: gap lineaire*** """)
        st.markdown(""" ***Avec paramètres:*** """)
        params = { 'mismatch':mismatch,'match':match,'open_gap':open_gap, 'extend_gap':extend_gap}
        
        st.write(params )
        st.write("SeqA : ",sequence1.upper())
        st.write("SeqB : ", sequence2.upper())
        
        # Gap lineaire : valeur par defaut
        alignments = pairwise2.align.globalmx(sequence1, sequence2, match=match,mismatch=mismatch) 
        for alignment in alignments:
            st.write(alignment)





    #Section 2
            #Multiple Alignment

if multiple_align_option == 'Exact method':

    col2.subheader('Multiple Alignment with ***An Exact method*** ')
    
    col2.markdown("""
            Cette section permet d'effectuer un alignement multiple  
            avec une methode Exacte \n
            ***Les mehtodes exactes sont des methode d'alignement permettant d'obtenir obligatoirement une solution optimale*** \n
            
            ***Comme Methode exacte nous avons implementer la **methde MUSCLE*****
            
            A l'instar de la methode clustal, cette methode neccessite deux fichier:
            1. L'excecutable Muscle.exe (Nous avons utilisé la version 3.8.31 )
            2. Un fichier .fasta contenant la liste séquences à aligner
            
            \n
            """)
    st.header('RESULT')

    # telechargement des fichiers neccessaire
    st.subheader("Selectionner le fichier contenant les sequences")
    
    sequence_file = st.file_uploader("***Sequence ***", type=["fasta"])
    if sequence_file is not None:
        st.write("...chargement sequence : OK")
    else:
        st.write('Awaiting fasta file to be uploaded for sequence 1. No file Upload')

    # Alignement local
    if sequence_file is not None:
        with open("file_muscle.fasta", "wb") as f:
            f.write(sequence_file.getbuffer())
        
        muscle_exe = "./muscle3.8.31_i86win32.exe"
        infile = "./file_muscle.fasta"
        muscle_cline = MuscleCommandline(muscle_exe, input=infile, out="file_muscle.fasta" )
        
        child = subprocess.Popen(str(muscle_cline),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        
        st.subheader("Télécharger le resultat de l'alignement ")
        
        with open('file_muscle.aln', 'r') as f:
            text = f.read()

        st.download_button(
            label="Download file",
            data=text,
            file_name='./file_muscle.aln',
            mime='text/plain'
        )
        
elif multiple_align_option == 'Clustal method':
    col2.subheader('Multiple Alignment with ***the Clustal method*** ')
    
    col2.markdown("""
            Cette section permet d'effectuer un alignement multiple  
            avec la methode Clustal \n
            
            Pour ce faire deux fichier sont necessaires:
            1. L'excecutable clustalw2.exe (Version clustalw2.exe)
            2. Un fichier .fasta contenant la liste séquences à aligner
            
            \n
            """)
    st.header('RESULT')

    # telechargement des fichiers neccessaire
    st.subheader("Selectionner le fichier contenant les sequences")

    sequence_file = st.file_uploader("***Sequence 1***", type=["fasta"])
    if sequence_file is not None:
        st.write("...chargement sequence 1: OK")
    else:
        st.write('Awaiting fasta file to be uploaded for sequence 1. No file Upload')

    # Alignement local
    if sequence_file is not None:
        with open("file.fasta", "wb") as f:
            f.write(sequence_file.getbuffer())
        
        clustalw_exe = "./clustalw2.exe"
        infile = "./file.fasta"
        clustalw_cline = ClustalwCommandline(clustalw_exe, infile=infile)
        
        child = subprocess.Popen(str(clustalw_cline),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        
        st.subheader("Télécharger le resultat de l'alignement ")
        
        with open('file.aln', 'r') as f:
            text = f.read()

        st.download_button(
            label="Download file",
            data=text,
            file_name='./file.aln',
            mime='text/plain'
        )
        
        st.subheader("Télécharger l'arbre correspondant à l'alignement ")
        
        with open('file.dnd', 'r') as f:
            text_dnd = f.read()

        st.download_button(
            label="Download file",
            data=text_dnd,
            file_name='./file.dnd',
            mime='text/plain'
        )
        
        with open("file.dnd", "r") as f:
            result = f.read()

elif multiple_align_option == 'Phylogenetic':
    
    col2.subheader('phylogenetic Tree with **WPGMA** ')
    
    col2.markdown("""
            Cette section permet de construire un arbre phylogenique en utilisant **l'agorithme WPGMA** en python \n
            
            Cette algorithme prend en entrer : la matrice de distance obtenue après avoir effectuer l'alignement multiple des séquences.
            
            ***(Les fichiers .aln que nous générons avec les algorithmes 
            Clustal et Muscle dans les autres sections permettent grâce la librairie Biopython de générer facilement cette matrice)***\n
            
            **Input**:
            
            1. Fichier contenant la matrice devra être renseigner avec un format txt
            2. La matrice dans le fichier doit être obligatoirement renseigner comme l'exemple de la figure ci-dessous
            
            **Exemple du format de matrice danns le fichier txt:**
            """)
    
    img_matrix = Image.open("./files/exampleMatrix.png")
    col2.image(img_matrix, width=200)
    
        # telechargement des fichiers neccessaire
    st.subheader("Selectionner le fichier contenant la matrice de distance")
    
    sequence_file = st.file_uploader("***matrice***", type=["txt"])
    if sequence_file is not None:
        st.write("...chargement sequence 1: OK")
    else:
        st.write('Awaiting fasta file to be uploaded for sequence 1. No file Upload')

    # Arbre phylogenique
    if sequence_file is not None:
        with open("WPGMA_Input.txt", "wb") as f:
            f.write(sequence_file.getbuffer())
            
        infile = "./WPGMA_Input.txt"
        matrix, length = readInput(infile)
        dictionary={}
        result_wpgma = wpgma(matrix, length,dictionary)
        tree_result=contruireArbre(dictionary,result_wpgma)
        tree_result = ''.join(tree_result)
        tree=tree_result+";"
        
        st.header('RESULT')
        st.subheader("Notation de Newick de l'Arbre phylogénetique obtenue: ")
        st.write(tree)
        st.subheader("Télecharger")
        st.download_button('Download', tree)
        
        st.header('Visualisation graphique')
        tree=Tree(tree)
                
        @contextmanager
        def st_capture(output_func):
            with StringIO() as stdout, redirect_stdout(stdout):
                old_write = stdout.write

                def new_write(string):
                    ret = old_write(string)
                    output_func(stdout.getvalue())
                    return ret
                
                stdout.write = new_write
                yield

    
        output = st.empty()
        with st_capture(output.code):
            print(f" ********************************************** {tree}")

