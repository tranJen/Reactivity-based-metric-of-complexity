from count_view_Reaction import *

import pandas as pd
import re
from rdkit import Chem
from rdkit.Chem import Descriptors



#---------STEP 1: General cleaning data-----------

def GeneralCleaning(data):
    '''Step 1 general cleaning'''
    # Only takes these information:
    columns_to_keep = [ 
        'Reaction', 
        'Reactant', 
        'Product', 
        'Reagent', 
        'Catalyst', 
        'Solvent (Reaction Details)',
        
        'Time (Reaction Details) [h]',
        'Temperature (Reaction Details) [C]',
        
        'Yield',
        
        'Reaction ID', 
        'Links to Reaxys',
        'Reaction: Links to Reaxys',
        'References'
    ]

    data = data[columns_to_keep]
    
    # Drop NA for Yield, Solvent, Reagent, Reaction, Time, Temperature
    data = data.dropna(subset=['Reagent','Solvent (Reaction Details)','Yield', 'Reaction', 'Temperature (Reaction Details) [C]','Time (Reaction Details) [h]'])

    # Remove rows that have catalyst
    data = data[data['Catalyst'].isna() == True]

    # Only keep row that has "Article" for "References""
    data = data[data['References'].str.contains('Article')]
    # Remove row having ';' or '-'
    data = data[ data['Time (Reaction Details) [h]'].str.contains(';|-') == False].copy()

    # Convert time to numeric
    data['Time (Reaction Details) [h]'] = pd.to_numeric(data['Time (Reaction Details) [h]'], errors='raise')

    # Remove yield reported with '>' (either >99 or >95 - remove)
    data = data[data['Yield'].str.contains('>') == False].copy()

    # Only get the row that has 'percent'
    data = data[data['Yield'].str.contains('percent', case = False)]

    # Drop duplicate rows
    data = data.drop_duplicates()
    return data

def extract_yield(string):
    ''' extract yield number from a string then add up all yield numbers'''
    # Extract yield
    numbers = re.findall(r'\d+\.\d+|\d+', string)
    # If there are 2 numbers, sum the yield (!!!: only for stereoisomers)
    if len(numbers) == 2: 
        return sum(float(num) for num in numbers)
    elif len(numbers) == 1:  # If there is only 1 number, return the number
        return float(numbers[0])
    else:
        return None

#---------STEP 2: Verify reaction by change in MW-----------

# Process the reaction SMILES
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

def calculate_mw(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Descriptors.MolWt(mol)
    else:
        return 0
    
def calculate_change_MW(reactants_list,products_list):
    '''calculate change in MW after the reaction'''
    reactants_MW = 0
    products_MW = 0
    valid_reaction = True
    # Calculate mass of reactants
    for reactant in reactants_list:
        if reactant is not None:
            mw = calculate_mw(reactant)
            if mw == 0:
                valid_reaction = False
                break
            reactants_MW += mw
            
    # Calculate mass of products       
    for product in products_list:
        if product is not None:
            mw = calculate_mw(product)
            if mw == 0:
                valid_reaction = False
                break
            products_MW += mw
    #Calculate change
    if valid_reaction:
        change_MW = products_MW - reactants_MW
    else:
        change_MW = -100 
    
    return change_MW

def get_largest_reactant_MW(reactants_list):
    '''calculate MW of largest reactant'''
    reactant_MW = 0
    # Calculate mass of reactants
    for reactant in reactants_list:
        if reactant_MW < calculate_mw(reactant):
            reactant_MW = calculate_mw(reactant)
    return reactant_MW


def get_change_MW(data):
    '''get change in MW after reaction and append to data frame'''
    data = data.dropna(subset=['Reactant SMILES', 'Product SMILES'])
    data['Change_MW'] = data.apply(
        lambda x: round(calculate_change_MW(x['Reactant SMILES'], x['Product SMILES']),3), axis=1)
        # remove invalid change_MW
    data = data[data['Change_MW'] != -100]
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

    data  = data[(data['Reactant SMILES'].apply(len) == 2) & (data['Product SMILES'].apply(len) == 2)]

    # Remove stereochemistry for the product 
    data = data.copy() #to avoid SettingWithCopyWarning
    data['Product SMILES'] = data['Product SMILES'].apply(lambda x: [i.replace('@', '') for i in x])

    # Keep the row if the 2 product SMILES (after removing stereochemistry) are the same.
    data = data[data['Product SMILES'].apply(lambda x: x[0]== x[1])]
    
     # Only keep rows having 2 yield-numbers reported by counting 'percent'
    data = data[data['Yield'].str.count('percent') == 2]

    #give back the strereochemistry for the product
    data = getReactandProduct(data)

    return data

#---------STEP 3: Verify reaction by change in C-O bonds-----------

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

#---------STEP 4: Ensure consistent reaction conditions-----------

def filter_same_condition(data,reagent, solvent, temp, less_than_time):

    data = data.copy()
    
    data['Reagent'] = data['Reagent'].astype(str)
    data = data[data['Reagent'].isin({reagent})]
    print(' - After filtering by reagent:')
    count_num_reaction(data)
    data= data[data['Solvent (Reaction Details)'] == solvent]
    print(' - After filtering by solvent:')
    count_num_reaction(data)
    data = data[data['Temperature (Reaction Details) [C]'] == temp]
    print(' - After filtering by temperature:')
    count_num_reaction(data)
    data = data[data['Time (Reaction Details) [h]'] < less_than_time ]
    print(' - After filtering by time:')
    count_num_reaction(data) 

    return data

#---------GENERALIZE 4 STEP ABOVE-----------

def step_1_to_4(data):

    # STEP 1: General cleaning data

    # clean
    data = GeneralCleaning(data)
    # extract yield (number)
    data = data.copy()
    data['Yield (number)'] = data['Yield'].apply(extract_yield)

    print('STEP 1 - general cleaning:')
    count_num_reaction(data)
    print('---------------------------------------')

    # STEP 2: Verify reaction by change in MW

    # get reactant and product SMILES from the reaction SMILES
    data = getReactandProduct(data)
    # Calculate change in MW after reaction
    data = get_change_MW(data)
    # For the single reduction: filter reaction have change MW by 2 and have 1 reactant, 1 product
    data_single = verify_SingleReduc_MW(data)
    # For the double reduction: filter reaction have change MW by 4 and have 2 stereoisomers products
    data_double = verify_SingleReduc_Stereo_MW(data)
    # Concatenate both data
    data = pd.concat([data_single, data_double], axis = 0)

    # Add a new column 'Largest Reactant MW' based on 'Reactant SMILES'
    data['Largest Reactant MW'] = data.apply(lambda x: round(get_largest_reactant_MW( x['Reactant SMILES']),3), axis=1)

    print('STEP 2 - Verify reaction by change in MW:')
    count_num_reaction(data)
    print('---------------------------------------')

    #STEP 3: Verify reaction by change in C-O bonds
   
    data = verify_SingleReduc_CObond(data)
    print('STEP 3 - Verify reaction by change in C-O bond:')
    count_num_reaction(data)
    print('---------------------------------------')

    # STEP 4: Ensure consistent reaction conditions
    print('STEP 4 - Ensure consistent reaction conditions:')
    data = filter_same_condition(data,'sodium tetrahydroborate', 'methanol', '0', 24)
    count_num_row(data)

    return data



