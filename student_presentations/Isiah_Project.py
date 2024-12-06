from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
import streamlit as st
import random
from PIL import Image
#Show the amino acids and their SMILES
data_folder = os.path.join('data')
print(data_folder)
file_path = os.path.join(data_folder, 'amino_acid_SMILES.txt')
df = pd.read_csv(file_path, skiprows=2)
df
AminoAcids = [Chem.MolFromSmiles(SMILES) for SMILES in df['SMILES']]
AminoAcids    
# Create a folder to save the images
output_folder = os.path.join(data_folder, 'amino_acid_images')
os.makedirs(output_folder, exist_ok=True)

# Iterate over the rows in the DataFrame
for index, row in df.iterrows():
    amino_acid = row['name']  # Assuming column name for amino acids
    smiles = row['SMILES']  # Assuming column name for SMILES strings
    
    # Convert the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    
    if mol:
        # Generate the 2D image for the molecule
        img = Draw.MolToImage(mol)
        
        # Define the output image path
        image_path = os.path.join(output_folder, f'{amino_acid}.png')
        
        # Save the image as PNG
        img.save(image_path)
        print(f'Saved image for {amino_acid} at {image_path}')
    else:
        print(f'Invalid SMILES for {amino_acid}: {smiles}')

# List of amino acids with their full names, three-letter codes, single-letter codes, and image file paths
amino_acids = [
    {"name": "Alanine", "three_letter": "Ala", "one_letter": "A", "image": "alanine.png"},
    {"name": "Arginine", "three_letter": "Arg", "one_letter": "R", "image": "arginine.png"},
    {"name": "Asparagine", "three_letter": "Asn", "one_letter": "N", "image": "asparagine.png"},
    {"name": "Aspartic acid", "three_letter": "Asp", "one_letter": "D", "image": "aspartate.png"},
    {"name": "Cysteine", "three_letter": "Cys", "one_letter": "C", "image": "cysteine.png"},
    {"name": "Glutamic acid", "three_letter": "Glu", "one_letter": "E", "image": "glutamate.png"},
    {"name": "Glutamine", "three_letter": "Gln", "one_letter": "Q", "image": "glutamine.png"},
    {"name": "Glycine", "three_letter": "Gly", "one_letter": "G", "image": "glycine.png"},
    {"name": "Histidine", "three_letter": "His", "one_letter": "H", "image": "histidine.png"},
    {"name": "Isoleucine", "three_letter": "Ile", "one_letter": "I", "image": "isoleucine.png"},
    {"name": "Leucine", "three_letter": "Leu", "one_letter": "L", "image": "leucine.png"},
    {"name": "Lysine", "three_letter": "Lys", "one_letter": "K", "image": "lysine.png"},
    {"name": "Methionine", "three_letter": "Met", "one_letter": "M", "image": "methionine.png"},
    {"name": "Phenylalanine", "three_letter": "Phe", "one_letter": "F", "image": "phenylalanine.png"},
    {"name": "Proline", "three_letter": "Pro", "one_letter": "P", "image": "proline.png"},
    {"name": "Serine", "three_letter": "Ser", "one_letter": "S", "image": "serine.png"},
    {"name": "Threonine", "three_letter": "Thr", "one_letter": "T", "image": "threonine.png"},
    {"name": "Tryptophan", "three_letter": "Trp", "one_letter": "W", "image": "tryptophan.png"},
    {"name": "Tyrosine", "three_letter": "Tyr", "one_letter": "Y", "image": "tyrosine.png"},
    {"name": "Valine", "three_letter": "Val", "one_letter": "V", "image": "valine.png"},
]

# Initialize session state
if 'score' not in st.session_state:
    st.session_state.score = 0
if 'question_index' not in st.session_state:
    st.session_state.question_index = 0
if 'shuffled_amino_acids' not in st.session_state:
    st.session_state.shuffled_amino_acids = random.sample(amino_acids, len(amino_acids))
if 'feedback' not in st.session_state:
    st.session_state.feedback = ""
if 'user_answer' not in st.session_state:
    st.session_state.user_answer = ""

# Helper function to check the user's answer
def check_answer(amino_acid, user_answer, quiz_type):
    user_answer = user_answer.strip().lower()  # Normalize the user input
    if quiz_type == "Name to One-letter Code":
        return user_answer.upper() == amino_acid["one_letter"]
    elif quiz_type == "Name to Three-letter Code":
        return user_answer.capitalize() == amino_acid["three_letter"]
    elif quiz_type == "Three-letter Code to Name":
        return user_answer == amino_acid["name"].lower()
    elif quiz_type == "One-letter Code to Name":
        return user_answer == amino_acid["name"].lower()
    elif quiz_type == "Identify by Structure":
        return user_answer == amino_acid["name"].lower()  # Match structure quiz answer
    return False

# Streamlit interface
st.title("Amino Acid Quiz")

# Quiz type selection
quiz_type = st.selectbox(
    "Select the quiz type:",
    ("Name to One-letter Code", "Name to Three-letter Code", "Three-letter Code to Name", "One-letter Code to Name", "Identify by Structure")
)

# Current question
if st.session_state.question_index < len(st.session_state.shuffled_amino_acids):
    current_amino_acid = st.session_state.shuffled_amino_acids[st.session_state.question_index]

    # Display the quiz question
    if quiz_type == "Name to One-letter Code":
        st.write(f"What is the one-letter code for {current_amino_acid['name']}?")
    elif quiz_type == "Name to Three-letter Code":
        st.write(f"What is the three-letter code for {current_amino_acid['name']}?")
    elif quiz_type == "Three-letter Code to Name":
        st.write(f"What is the name of the amino acid with the three-letter code {current_amino_acid['three_letter']}?")
    elif quiz_type == "One-letter Code to Name":
        st.write(f"What is the name of the amino acid with the one-letter code {current_amino_acid['one_letter']}?")
    elif quiz_type == "Identify by Structure":
        st.write("Identify the amino acid by its structure:")
        try:
            img = Image.open(current_amino_acid["image"])
            st.image(img, caption="Amino Acid Structure")
        except FileNotFoundError:
            st.write("Image not found. Please check the file path.")

    # User answer input
    user_answer = st.text_input("Your answer:", value=st.session_state.user_answer, key="answer_input")

    # Submit button logic
    if st.button("Submit"):
        if check_answer(current_amino_acid, user_answer, quiz_type):
            st.session_state.feedback = "Correct!"
            st.session_state.score += 1
        else:
            correct_answer = {
                "Name to One-letter Code": current_amino_acid["one_letter"],
                "Name to Three-letter Code": current_amino_acid["three_letter"],
                "Three-letter Code to Name": current_amino_acid["name"],
                "One-letter Code to Name": current_amino_acid["name"],
                "Identify by Structure": current_amino_acid["name"]
            }[quiz_type]
            st.session_state.feedback = f"Incorrect. The correct answer was {correct_answer}."
        st.session_state.user_answer = ""  # Clear input field after submission

    # Display feedback
    if st.session_state.feedback:
        if "Correct!" in st.session_state.feedback:
            st.markdown('<p style="color:green;">{}</p>'.format(st.session_state.feedback), unsafe_allow_html=True)
        else:
            st.markdown('<p style="color:red;">{}</p>'.format(st.session_state.feedback), unsafe_allow_html=True)

# Next button logic
if st.button("Next"):
    if st.session_state.question_index < len(st.session_state.shuffled_amino_acids) - 1:
        st.session_state.question_index += 1
        st.session_state.feedback = ""  # Clear feedback for the next question
        st.session_state.user_answer = ""  # Reset the input field
        st.experimental_rerun()  # Force rerun to display the next question
    else:
        st.write(f"Quiz completed! Your final score: {st.session_state.score}/{len(amino_acids)}")
        if st.button("Restart Quiz"):
            st.session_state.score = 0
            st.session_state.question_index = 0
            st.session_state.shuffled_amino_acids = random.sample(amino_acids, len(amino_acids))
            st.experimental_rerun()
