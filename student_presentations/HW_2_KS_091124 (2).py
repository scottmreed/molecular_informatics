from rdkit import Chem
from rdkit.Chem import AllChem
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# 1 (a)
morphine = Chem.MolFromSmiles('CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O')
morphine
# 1 (a)
Chem.Draw.MolToFile(morphine, 'morphine.svg', size=(500,500),imageType='svg')
# 1 (b)
Chiral_Center_Count = Chem.FindMolChiralCenters(morphine)
len(Chiral_Center_Count)
# 1 (c)
Chem.rdMolDescriptors.CalcNumHBA(morphine)
# 2
data_folder = os.path.join('data')
file_path = os.path.join(data_folder, 'amino_acid_SMILES.txt')
df = pd.read_csv(file_path, skiprows=2)
df
# 2 (a)
AminoAcids = [Chem.MolFromSmiles(SMILES) for SMILES in df['SMILES']]
AminoAcids
# 2 (a) - used ChatGPT; couldn't figure out how to get chirality from a list of SMILES, but could do it individually

from rdkit import Chem

def determine_absolute_configuration(smiles: str) -> str:
    """
    Determines the absolute configuration (R/S) of the alpha-carbon
    for the given SMILES string representing an amino acid.
    """
    # Parse the molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return "Invalid SMILES"
    
    # Ensure stereochemistry is assigned properly
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    
    # Identify the chiral centers (returns a list of (atom_index, R/S/?))
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False, useLegacyImplementation=False)
    
    # Return the configuration (R or S) of the first chiral center (alpha-carbon)
    if chiral_centers:
        # We're only interested in the alpha-carbon, which is typically the first chiral atom
        atom_idx, chirality = chiral_centers[0]
        return chirality  # 'R' or 'S'
    else:
        return "No chiral center found"

# List of SMILES strings for common amino acids
amino_acids_smiles = {
    'Alanine': 'CC[C@H](C(=O)O)N',
    'Cysteine': 'C([C@@H](C(=O)O)N)S',
    'Aspartic acid': 'C([C@H](C(=O)O)N)C(=O)O',
    'Glutamic acid': 'C(CC([C@H](C(=O)O)N)C(=O)O)',
    'Phenylalanine': 'C([C@H](C(=O)O)N)Cc1ccccc1',
    'Glycine': 'C(C(=O)O)N',  # Not chiral
    'Histidine': 'C([C@H](C(=O)O)N)Cc1c[nH]cn1',
    'Isoleucine': 'CC[C@@H](C)[C@H](C(=O)O)N',
    'Leucine': 'CC(C)CC([C@H](C(=O)O)N)',
    'Lysine': 'C(CCN)CC([C@H](C(=O)O)N)',
    'Methionine': 'C([C@H](C(=O)O)N)CSCC',
    'Asparagine': 'C([C@H](C(=O)O)N)C(=O)N',
    'Proline': 'C1CC(NC1)C(=O)O',  # Technically chiral, but often treated as achiral
    'Glutamine': 'C(CC([C@H](C(=O)O)N)C(=O)N)',
    'Arginine': 'C(CCCN=C(N)N)C([C@H](C(=O)O)N)',
    'Serine': 'C([C@H](C(=O)O)N)O',
    'Threonine': 'C([C@@H](C(=O)O)N)O',
    'Tryptophan': 'C([C@H](C(=O)O)N)Cc1c[nH]c2ccccc12',
    'Tyrosine': 'C([C@H](C(=O)O)N)Cc1ccc(O)cc1',
    'Valine': 'CC[C@H](C(=O)O)N',
}

# Iterate over each amino acid and determine its absolute configuration
for name, smiles in amino_acids_smiles.items():
    config = determine_absolute_configuration(smiles)
    print(f"{name}: {config}")
ChatGPT has opposite configuration outputs for chirality, suggesting the SMILEs for each amino acid were of the D-isomers instead of L-isomers
# 2 (a) individually done instead of through list
alanine	= Chem.MolFromSmiles('C[C@@H](C(=O)[O-])[NH3+]')
arginine = Chem.MolFromSmiles('[NH3+][C@@H](CCCNC(=[NH2+])N)C(=O)[O-]')
asparagine = Chem.MolFromSmiles('O=C(N)C[C@H]([NH3+])C(=O)[O-]')
aspartate = Chem.MolFromSmiles('C([C@@H](C(=O)[O-])[NH3+])C(=O)[O-]')
cysteine = Chem.MolFromSmiles('C([C@@H](C(=O)[O-])[NH3+])S')
glutamine = Chem.MolFromSmiles('[NH3+][C@@H](CCC(=O)N)C([O-])=O')
glutamate = Chem.MolFromSmiles('C(CC(=O)[O-])[C@@H](C(=O)[O-])[NH3+]')
glycine = Chem.MolFromSmiles('C(C(=O)[O-])[NH3+]')
histidine = Chem.MolFromSmiles('O=C([C@H](CC1=CNC=N1)[NH3+])[O-]')
isoleucine = Chem.MolFromSmiles('CC[C@H](C)[C@@H](C(=O)[O-])[NH3+]')
leucine = Chem.MolFromSmiles('CC(C)C[C@@H](C(=O)[O-])[NH3+]')
lysine = Chem.MolFromSmiles('C(CC[NH3+])C[C@@H](C(=O)[O-])[NH3+]')
methionine = Chem.MolFromSmiles('CSCC[C@H]([NH3+])C(=O)[O-]')
phenylalanine = Chem.MolFromSmiles('[NH3+][C@@H](CC1=CC=CC=C1)C([O-])=O')
proline	= Chem.MolFromSmiles('[O-]C(=O)[C@H](CCC2)[NH2+]2')
serine = Chem.MolFromSmiles('C([C@@H](C(=O)[O-])[NH3+])O')
threonine = Chem.MolFromSmiles('C[C@H]([C@@H](C(=O)[O-])[NH3+])O')
tryptophan = Chem.MolFromSmiles('c1[nH]c2ccccc2c1C[C@H]([NH3+])C(=O)[O-]')
tyrosine = Chem.MolFromSmiles('[NH3+][C@@H](Cc1ccc(O)cc1)C([O-])=O')
valine = Chem.MolFromSmiles('CC(C)[C@@H](C(=O)[O-])[NH3+]')

Chem.FindMolChiralCenters(alanine)
Chem.FindMolChiralCenters(asparagine)
Chem.FindMolChiralCenters(aspartate)
Chem.FindMolChiralCenters(cysteine)
Chem.FindMolChiralCenters(glutamate)
Chem.FindMolChiralCenters(glutamine)
Chem.FindMolChiralCenters(glycine)
Chem.FindMolChiralCenters(histidine)
Chem.FindMolChiralCenters(isoleucine)
Chem.FindMolChiralCenters(leucine)
Chem.FindMolChiralCenters(lysine)
Chem.FindMolChiralCenters(methionine)
Chem.FindMolChiralCenters(phenylalanine)
Chem.FindMolChiralCenters(proline)
Chem.FindMolChiralCenters(serine)
Chem.FindMolChiralCenters(threonine)
Chem.FindMolChiralCenters(tryptophan)
Chem.FindMolChiralCenters(tyrosine)
Chem.FindMolChiralCenters(valine)
2 (a)/(b)

All amino acids, except for glycine (no chiral center, and not a chiral molecule) and cysteine (R) are the S-conformation at their alpha-carbons. Cysteine is the R conformation because the sulfhydral group has higher priority. 

Isoleucine and threonine have two chiral carbons.
# 3
data_folder = os.path.join('data')
file_path = os.path.join(data_folder, 'organic_molecules.txt')
df = pd.read_csv(file_path, skiprows=2)
df
# 3 (a) Couldn't figure out how to determine this using a list, so used ChatGPT, but didn't work

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_primary_or_secondary_alcohol(smiles: str) -> bool:
    """
    Determines if a molecule is a primary or secondary aliphatic alcohol.
    """
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return False  # Invalid SMILES, ignore
    
    # Iterate over atoms to find the hydroxyl group (-OH)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen atom
            # Check if it is part of a hydroxyl group
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 1 and neighbors[0].GetAtomicNum() == 6:  # Attached to a single carbon
                # Get the carbon it's attached to
                carbon = neighbors[0]
                # Count the number of carbon neighbors (for primary or secondary alcohol)
                carbon_neighbors = [nbr for nbr in carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
                if len(carbon_neighbors) == 1:  # Primary alcohol (1 other carbon attached)
                    return True
                elif len(carbon_neighbors) == 2:  # Secondary alcohol (2 other carbons attached)
                    return True
    return False  # Not a primary or secondary alcohol

def filter_alcohols_from_smiles_file(filename: str) -> list:
    """
    Reads a file of SMILES strings and returns a list of SMILES that are
    primary or secondary aliphatic alcohols.
    """
    primary_secondary_alcohols = []
    
    # Read the SMILES strings from the file
    with open(filename, 'r') as f:
        smiles_list = f.readlines()
    
    # Iterate through the SMILES strings and filter primary/secondary alcohols
    for smiles in smiles_list:
        smiles = smiles.strip()  # Clean up any extra whitespace/newlines
        if is_primary_or_secondary_alcohol(smiles):
            primary_secondary_alcohols.append(smiles)
    
    return primary_secondary_alcohols

# Example usage
# filename = 'organic_molecules.txt'  # Replace with your actual filename
# result = filter_alcohols_from_smiles_file(filename)
# print(result)


# 3 (b)
Substructure = Chem.MolFromSmiles('CO')
Substructure
A = Chem.MolFromSmiles('CCC(C)CC')
B = Chem.MolFromSmiles('OCC')
C = Chem.MolFromSmiles('NC(CCO1)CC1')
D = Chem.MolFromSmiles('NCC(=CC=C1)C=C1')
E = Chem.MolFromSmiles('NC(=O)C(=CC=C1)C=C1')
F = Chem.MolFromSmiles('NC(C)C(C)O')
G = Chem.MolFromSmiles('NCC(CC)CC')
H = Chem.MolFromSmiles('NCC(C)CC')
I = Chem.MolFromSmiles('NC')
J = Chem.MolFromSmiles('SC(CCC1)CC1')
K = Chem.MolFromSmiles('SCC(C=CC1)=CC=1')
L = Chem.MolFromSmiles('OCCCCC')
M = Chem.MolFromSmiles('OC(CCO1)CC1')
N = Chem.MolFromSmiles('OCC(C=CC1)=CC=1')
O = Chem.MolFromSmiles('OC(=O)C(=CC=C1)C=C1')
P = Chem.MolFromSmiles('O[C@@H](CC[C@]1(C)[C@@H](CC[C@]2(C)[C@H]3[C@H](C)CCCC(C)C)[C@@H]4[C@@H]2CC3)CC1=CC4')
Q = Chem.MolFromSmiles('OC(=CC1)C(=CC=1)C(O)=O')
R = Chem.MolFromSmiles('O=C(O)/C=C/C(O)=O')
S = Chem.MolFromSmiles('O=C(C)OCC')
T = Chem.MolFromSmiles('O=C(C)C')
U = Chem.MolFromSmiles('O=C(C)OC(=O)C')
V = Chem.MolFromSmiles('OCC(=O)[C@@H](O)[C@H](O)[C@H](O)CO')
A.HasSubstructMatch(Substructure)
B.HasSubstructMatch(Substructure)
C.HasSubstructMatch(Substructure)
D.HasSubstructMatch(Substructure)
E.HasSubstructMatch(Substructure)
F.HasSubstructMatch(Substructure)
G.HasSubstructMatch(Substructure)
H.HasSubstructMatch(Substructure)
I.HasSubstructMatch(Substructure)
J.HasSubstructMatch(Substructure)
K.HasSubstructMatch(Substructure)
L.HasSubstructMatch(Substructure)
M.HasSubstructMatch(Substructure)
N.HasSubstructMatch(Substructure)
O.HasSubstructMatch(Substructure)
P.HasSubstructMatch(Substructure)
Q.HasSubstructMatch(Substructure)
R.HasSubstructMatch(Substructure)
S.HasSubstructMatch(Substructure)
T.HasSubstructMatch(Substructure)
U.HasSubstructMatch(Substructure)
V.HasSubstructMatch(Substructure)
# 4
DMC = Chem.MolFromSmiles('CC1CCCCC1C')
isomers = list(Chem.EnumerateStereoisomers.EnumerateStereoisomers(DMC))
Chem.Draw.MolsToGridImage(isomers, useSVG=True)