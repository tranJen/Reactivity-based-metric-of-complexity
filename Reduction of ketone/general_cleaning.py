import pandas as pd
import re
from rdkit.Chem import Draw
from rdkit.Chem import rdChemReactions

def count_num_reaction(data):
    '''Count how many reactions there are'''
    rxnCount_data = data['Reaction ID'].nunique()
    print('Number of Reactions:', rxnCount_data)
    print('Number of Rows:', data.shape[0])

def view_reactionScheme(data, NumReaction_to_view):
    ''' randomly pick and show reaction scheme'''
    # Get 1 sample for each reaction ID, Remove duplicated ID
    Reaction_data = data.drop_duplicates(subset=['Reaction ID'], keep='first')
    Reaction_data = Reaction_data.reset_index(drop = True)
    
    if NumReaction_to_view > Reaction_data.shape[0]:
        raise ValueError('Number of reactions to view is more than the total number of reactions in the dataset')
    else:
        # Draw
        for idx, row in Reaction_data.sample(n=NumReaction_to_view).iterrows():
            reaction_smarts = row['Reaction']
            rxn = rdChemReactions.ReactionFromSmarts(reaction_smarts, useSmiles=True)
            if rxn:
                print('Reaction ID:', Reaction_data.iloc[idx]['Reaction ID'] )
                Draw.ReactionToImage(rxn).show() # image pop up
                print(f'Link: {row["Links to Reaxys"]}\n')

def GeneralCleaning(data):
    '''Step 1 general cleaning'''
    # Drop NA for Yield, Solvent, Reagent
    data = data.dropna(subset=['Reagent','Solvent (Reaction Details)','Yield', 'Reaction'])

    # Remove rows that have catalyst
    data = data[data['Catalyst'].isna() == True]
    
    # Only get the row that has 'percent'
    data = data[data['Yield'].str.contains('percent', case = False)]

    # Only keep row that has "Article" for "References""
    data = data[data['References'].str.contains('Article')]
    
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