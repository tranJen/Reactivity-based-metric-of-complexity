import pandas as pd
from rdkit import Chem

def count_C_O_bonds(molecule_SMILES, bond_type):
    mol = Chem.MolFromSmiles(molecule_SMILES)
    num_bonds = 0
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 8) or (a1.GetAtomicNum() == 8 and a2.GetAtomicNum() == 6):
            if bond.GetBondType() == bond_type:
                num_bonds += 1            
    return num_bonds

def change_C_O_bonds(reactant_SMILES, product_SMILES, bond_type):
    change = count_C_O_bonds(product_SMILES, bond_type) - count_C_O_bonds(reactant_SMILES, bond_type)
    return change

def change_single_and_double_C_O_bond(data):
    data['change in C-O single bond'] = data.apply(
        lambda x: change_C_O_bonds(x['Reactant SMILES'][0], x['Product SMILES'][0], Chem.rdchem.BondType.SINGLE), axis=1)
    data['change in C=O double bond'] = data.apply(
        lambda x: change_C_O_bonds(x['Reactant SMILES'][0], x['Product SMILES'][0], Chem.rdchem.BondType.DOUBLE), axis=1)
    return data

def verify_SingleReduc_CObond(data):
    data = change_single_and_double_C_O_bond(data)
    data = data[(data['change in C-O single bond'] == 1) & (data['change in C=O double bond'] == -1)]
    return data
