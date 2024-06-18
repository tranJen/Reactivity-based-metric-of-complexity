
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

def getReactandProduct(data):
    '''take a dataframe have reaction SMILES as column named 'Reaction' 
    extract reactants and products into a lists'''
    
    data = data.copy() #to avoid SettingWithCopyWarning
    
    # Split reactions into reactants and product 
    data[['Reactant SMILES', 'Product SMILES']] = data['Reaction'].str.split('>>', expand=True)

    # Split reactants and products into list
    data['Reactant SMILES'] = data['Reactant SMILES'].astype(str)
    data['Product SMILES'] = data['Product SMILES'].astype(str)
    data['Reactant SMILES'] = data['Reactant SMILES'].apply(lambda x: x.split('.'))
    data['Product SMILES'] = data['Product SMILES'].apply(lambda x: x.split('.'))
    return data

def calculate_change_MW(reactants_list,products_list):
    '''calculate change in MW after the reaction'''
    reactants_MW = 0
    products_MW = 0
    # Calculate mass of reactants
    for reactant in reactants_list:
        reactant_molecule = Chem.MolFromSmiles(reactant)
        if reactant_molecule is not None:
            reactants_MW += Descriptors.MolWt(reactant_molecule)
    # Calculate mass of products       
    for product in products_list:
        product_molecule = Chem.MolFromSmiles(product)
        if product_molecule is not None:
            products_MW += Descriptors.MolWt(product_molecule)
    #Calculate change
    change_MW = products_MW - reactants_MW
    return change_MW


def get_change_MW(data):
    '''get change in MW after reaction and append to data frame'''
    data = data.dropna(subset=['Reactant SMILES', 'Product SMILES'])
    data['Change_MW'] = data.apply(
        lambda x: round(calculate_change_MW(x['Reactant SMILES'], x['Product SMILES']),3), axis=1)
    return data

def verify_SingleReduc_MW(data):
    '''filter reaction have change MW by 2 and have 1 reactant, 1 product'''
    # For the single reduction, verified by MW:
    data = data[data['Change_MW'] == 2.016]
    #only take those have 1 reactant and 1 product
    data = data[(data['Reactant SMILES'].apply(len) == 1) & (data['Product SMILES'].apply(len) == 1)]
    return data

def verify_SingleReduc_Stereo_MW(data):
    '''filter reaction have change MW by 4 and have 2 stereoisomers products'''
    # For the double reduction, verified by MW:
    data = data[data['Change_MW'] == 4.032]

    # Only take those have 2 reactant and 2 product 
    # and products' SMILES are the same after removing stereochemistry
    # These are reactions making stereoisomers product

    data  = data[(data['Reactant SMILES'].apply(len) == 2) | 
        (data['Product SMILES'].apply(len) == 2)]

        # Remove stereochemistry for the product 
    data = data.copy() #to avoid SettingWithCopyWarning
    data['Product SMILES'] = data['Product SMILES'].apply(lambda x: [i.replace('@', '') for i in x])

        # Keep the row if the 2 product SMILES (after removing stereochemistry) are the same.
    data = data[data['Product SMILES'].apply(lambda x: x[0]== x[1])]
    
        # Only keep rows having 2 yield-numbers reported by counting 'percent'
    data = data[data['Yield'].str.count('percent') == 2]
    return data
        
