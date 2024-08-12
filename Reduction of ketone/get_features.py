import os
# Change working directory
os.chdir('/Users/suongsuong/Documents/GitHub/Reactivity-based-metric-of-complexity/Reduction of ketone/')

import pandas as pd
import numpy as np
from count_view_Reaction import *

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors,GraphDescriptors, AllChem

import dbstep.Dbstep as db
from calculate_published_features.bottchscore3 import *      #developed by Forli Lab https://github.com/forlilab/bottchscore
import calculate_published_features.spacial_score as sps     #from GitHub page of published paper https://github.com/frog2000/Spacial-Score
from calculate_published_features.pychem.src.pychem.topology import * # from https://github.com/cosylabiiit/chemopy?tab=readme-ov-file


### ------ GET BASIC FEATURES ----------

def count_HA(smiles: str) -> int:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("Invalid SMILES string")
    
    Num_heteroatoms = 0
    for atom in molecule.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6]: 
            Num_heteroatoms += 1
    return Num_heteroatoms

def count_sp3_C(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    Num_sp3_carbons = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
            Num_sp3_carbons += 1
    return Num_sp3_carbons

def count_C(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    Num_carbons = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            Num_carbons += 1
    return Num_carbons

def count_stereogenic_C(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    Num_chiral_centers = rdMolDescriptors.CalcNumAtomStereoCenters(mol)
    return Num_chiral_centers

def count_NumRing(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    NumRing = rdMolDescriptors.CalcNumRings(mol)
    return NumRing

def count_Num_aroRing(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    Num_aroRing = rdMolDescriptors.CalcNumAromaticRings(mol)
    return Num_aroRing

def count_Num_ketone(smiles: str) -> int:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    Num_ketone = 0
    for bond in mol.GetBonds():
        # Check if the bond is a double bond and if it connects a carbon to an oxygen
        if bond.GetBondTypeAsDouble() == 2:
            atoms = (bond.GetBeginAtom(), bond.GetEndAtom())
            if (atoms[0].GetAtomicNum() == 6 and atoms[1].GetAtomicNum() == 8) or \
               (atoms[1].GetAtomicNum() == 6 and atoms[0].GetAtomicNum() == 8):
                Num_ketone += 1
    return Num_ketone

def get_partialCharge_carbonylC(smiles, carbonylC_idx):
    molecule = Chem.MolFromSmiles(smiles)
    AllChem.ComputeGasteigerCharges(molecule)
    atom = molecule.GetAtomWithIdx(carbonylC_idx - 1) #idex in table count from 1
    charge = atom.GetProp('_GasteigerCharge')
    return charge

def get_basic_features(data, col_smiles, carbonyl_idx):
    data = data.copy()
    data['Ring'] = data.apply(lambda x: count_NumRing( x[col_smiles]), axis=1)
    data['aroma_Ring'] = data.apply(lambda x: count_Num_aroRing( x[col_smiles]), axis=1)

    data['Carbons'] = data.apply(lambda x: count_C( x[col_smiles]), axis=1)
    data['Chiral Carbons'] = data.apply(lambda x: count_stereogenic_C( x[col_smiles]), axis=1)

    data['sp3 Carbons'] = data.apply(lambda x: count_sp3_C( x[col_smiles]), axis=1)
    data['HA'] = data.apply(lambda x: count_HA( x[col_smiles]), axis=1)

    data['Number of ketone'] = data.apply(lambda x: count_Num_ketone( x[col_smiles]), axis=1)
    data['Partial Charge of C carbonyl'] = data.apply(lambda x: get_partialCharge_carbonylC( x[col_smiles], x[carbonyl_idx]), axis=1)
    return data

### ------ GET published index ----------

def calculate_bertz_index(smiles):
    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    bertz_index = GraphDescriptors.BertzCT(mol)
    return bertz_index

def get_published_index(data, col_smiles):
    data['Cm'] = data.apply(lambda x: calculate_bottchscore_from_smiles( x[col_smiles],verbose_response=False, debug_arg=False, disable_mesomer=False, automorp_memory_maxsize=3000000), axis=1)
    data['Cm/HA'] = data.apply(lambda x: x['Cm']/x['HA'], axis=1)

    data['SPS'] = data.apply(lambda x: sps.calculate_score_from_smiles(x[col_smiles], per_atom=False), axis=1)
    data['nSPS'] = data.apply(lambda x:sps.calculate_score_from_smiles(x[col_smiles], per_atom=True), axis=1)

    data['F_sp3'] = data.apply(lambda x: x['sp3 Carbons']/x['Carbons'], axis=1)
    data['F_Cstereo'] = data.apply(lambda x: x['Chiral Carbons']/x['Carbons'], axis=1)
    
    data['C_T'] = data.apply(lambda x: calculate_bertz_index(x[col_smiles]), axis=1)

    return data

### ------ GET topological index ----------

def get_topology_index(data, col_smiles):
    
    data['Harary'] = data.apply(lambda x: CalculateHarary(Chem.MolFromSmiles(x[col_smiles])), axis=1)
    data['GraphDistance'] = data.apply(lambda x: CalculateGraphDistance(Chem.MolFromSmiles(x[col_smiles])), axis=1)
    data['Diameter'] = data.apply(lambda x: CalculateDiameter(Chem.MolFromSmiles(x[col_smiles])), axis=1)
    data['Platt'] = data.apply(lambda x: CalculatePlatt(Chem.MolFromSmiles(x[col_smiles])), axis=1)
    data['SimpleTopo'] = data.apply(lambda x: CalculateSimpleTopoIndex(Chem.MolFromSmiles(x[col_smiles])), axis=1)
    data['GeometricTopo'] = data.apply(lambda x: CalculateGeometricTopoIndex(Chem.MolFromSmiles(x[col_smiles])), axis=1)
    data['ArithmeticTopo'] = data.apply(lambda x: CalculateArithmeticTopoIndex(Chem.MolFromSmiles(x[col_smiles])), axis=1)
    data['MolSizeTotalInf'] = data.apply(lambda x: CalculateMolSizeTotalInf(Chem.MolFromSmiles(x[col_smiles])), axis=1)

    return data

### ------ GET FEATURES at alpha positions----------

def find_alpha_carbon_index(smiles, carbonyl_C_index):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")

    atom = mol.GetAtomWithIdx(carbonyl_C_index - 1) #index in table counting from 1
    atom_neighbors = atom.GetNeighbors()

    alpha_atom_indices = []
    for neighbor in atom_neighbors:
        if neighbor.GetSymbol() != 'O':
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
            if bond and bond.GetBondType() != Chem.BondType.DOUBLE:  #not getting the oxygen of the carbonyl
                neighbor_symbol = neighbor.GetSymbol()
                neighbor_idx = neighbor.GetIdx() + 1
                neighbor_info = [neighbor_symbol,neighbor_idx]
                alpha_atom_indices.append(neighbor_info) 
    return alpha_atom_indices #return a list of [atom, idx]

def get_isAromatic_atAlpha(smiles, index_list):
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    is_aromatic = 0
    for index in index_list:
        index -= 1 #  Convert index to 0-based if necessary (RDKit uses 0-based indexing)
        atom = mol.GetAtomWithIdx(index) 

        if atom.GetIsAromatic():# Check for aromaticity
            is_aromatic += 1
        else:
            is_aromatic += 0
    return is_aromatic

def get_subs_group_atAlpha(smiles, index):
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    index -= 1#  Convert index to 0-based if necessary (RDKit uses 0-based indexing)
    atom = mol.GetAtomWithIdx(index)

    # Count the number of bond on alpha position
    num_bonds = atom.GetDegree()
    if num_bonds == 4:
        return 4 - 1 #minus 1 of the main branch
    elif num_bonds == 3:
        return 3 - 1
    elif num_bonds == 2:
        return 2 - 1 
    elif num_bonds == 1:
        return 0
    else:
        raise ValueError(smiles, "Number of substituent groups cannot be determined")

def get_alpha_features(data, col_smiles):
    data = data.copy()
    data['C_alpha_indexes'] = data.apply(lambda x: find_alpha_carbon_index(x[col_smiles], x['C_idx']), axis=1)
    data['C_alpha_index 1'] = data.apply(lambda x: x['C_alpha_indexes'][0][1], axis = 1)
    data['C_alpha_index 2'] = data.apply(lambda x: x['C_alpha_indexes'][1][1], axis = 1)

    data['Number aromatic ring at alpha'] = data.apply(lambda x: get_isAromatic_atAlpha(x[col_smiles], [x['C_alpha_index 1'],x['C_alpha_index 2']]), axis=1)
    data['Number sub group at index 1'] = data.apply(lambda x: get_subs_group_atAlpha(x[col_smiles], x['C_alpha_index 1']), axis=1)
    data['Number sub group at index 2'] = data.apply(lambda x: get_subs_group_atAlpha(x[col_smiles], x['C_alpha_index 2']), axis=1)
    data['Number of sub group both side'] = data['Number sub group at index 1'] + data['Number sub group at index 2']
    return data

### ------ GET Sterimol, PBV ----------

def calculate_Sterimol_fromXYZ(Reaction_ID, atom1_idx, atom2_idx, folder_path):
    file = folder_path + str(Reaction_ID) + '.xyz'
    if file:
        mol = db.dbstep(file, atom1=atom1_idx, atom2=atom2_idx, commandline=True, verbose=True, sterimol=True, volume=False, measure='grid')
        L = mol.L
        Bmin = mol.Bmin
        Bmax = mol.Bmax
        return L, Bmin, Bmax
    else:
        return None, None, None

def calculate_PBV(Reaction_ID, atom1_idx, atom2_idx, folder_path):
    file = folder_path + str(Reaction_ID) + '.xyz'
    if file:
        mol = db.dbstep(file,atom1=atom1_idx,atom2=atom2_idx,commandline=True,verbose=True,sterimol=False, volume = True, measure='grid') 
        PBV = mol.bur_vol
        return PBV
    else:
        return None

###  At C carbonyl 

def calculate_SterimolPBV_OtoC_data(data, carbonylC_idxColumn, carbonylO_idxColumn, folder_path):
    #set atom1 = O and atom 2 = C
    data[['L','Bmin', 'Bmax']] = data.apply(lambda x: pd.Series(calculate_Sterimol_fromXYZ(x['Reaction ID'], x[carbonylO_idxColumn], x[carbonylC_idxColumn], folder_path)), axis = 1)
    #set atom1 = C (doesn't matter atom 2)
    data['PBV'] = data.apply(lambda x: calculate_PBV(x['Reaction ID'], x[carbonylC_idxColumn], x[carbonylO_idxColumn], folder_path), axis = 1)
    return data

###  At alpha positions

def calculate_SterimolPBV_CtoAlpha_data(data, carbonylC_idxColumn, folder_path):
    #set atom1 = C-carbonyl and atom 2 = atom at alpha
    data[['L_alpha1','Bmin_alpha1', 'Bmax_alpha1']] = data.apply(lambda x: pd.Series(calculate_Sterimol_fromXYZ(x['Reaction ID'], x[carbonylC_idxColumn], x['C_alpha_index 1'], folder_path)), axis = 1)
    data[['L_alpha2','Bmin_alpha2', 'Bmax_alpha2']] = data.apply(lambda x: pd.Series(calculate_Sterimol_fromXYZ(x['Reaction ID'], x[carbonylC_idxColumn], x['C_alpha_index 2'], folder_path)), axis = 1)
    data['L_alpha Sum'] = data.apply(lambda x: (x['L_alpha1'] + x['L_alpha2']), axis = 1)
    data['Bmin_alpha Sum'] = data.apply(lambda x: (x['Bmin_alpha1'] + x['Bmin_alpha2']), axis = 1)
    data['Bmax_alpha Sum'] = data.apply(lambda x: (x['Bmax_alpha1'] + x['Bmax_alpha2']), axis = 1)

    #set atom1 = atom at alpha (doesn't matter atom 2)
    data['PBV_alpha1'] = data.apply(lambda x: calculate_PBV(x['Reaction ID'], x['C_alpha_index 1'], x[carbonylC_idxColumn], folder_path), axis = 1)
    data['PBV_alpha2'] = data.apply(lambda x: calculate_PBV(x['Reaction ID'], x['C_alpha_index 2'], x[carbonylC_idxColumn], folder_path), axis = 1)
    data['PBV_alpha Sum'] = data.apply(lambda x: (x['PBV_alpha1'] + x['PBV_alpha2']), axis = 1)
    data['PBV_alpha Diff.'] = data.apply(lambda x: np.abs(x['PBV_alpha1'] - x['PBV_alpha2']), axis = 1)
    return data

### ---- calculate feature all in 1 function ----

def calculate_features(data, col_smiles, carbonylC_idxColumn, carbonylO_idxColumn,  folder3D_path):
    data = get_basic_features(data, col_smiles, carbonylC_idxColumn)
    data = get_published_index(data, col_smiles)
    data = get_topology_index(data, col_smiles)
    data = get_alpha_features(data, col_smiles)
    data = calculate_SterimolPBV_OtoC_data(data, carbonylC_idxColumn, carbonylO_idxColumn, folder3D_path)
    data = calculate_SterimolPBV_CtoAlpha_data(data, carbonylC_idxColumn, folder3D_path)
    return data