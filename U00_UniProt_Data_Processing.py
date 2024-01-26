#!/usr/bin/env python
# coding: utf-8
#====================================================================================================#
# The following code ensures the code work properly in 
# MS VS, MS VS CODE and jupyter notebook on both Linux and Windows.
import os 
import sys
from os import path
from sys import platform
from pathlib import Path

if __name__ == "__main__":
    print("\n\n")
    print("="*80)
    if os.name == 'nt' or platform == 'win32':
        print("Running on Windows")
        if 'ptvsd' in sys.modules:
            print("Running in Visual Studio")

    if os.name != 'nt' and platform != 'win32':
        print("Not Running on Windows")

    if "__file__" in globals().keys():
        print('CurrentDir: ', os.getcwd())
        try:
            os.chdir(os.path.dirname(__file__))
        except:
            print("Problems with navigating to the file dir.")
        print('CurrentDir: ', os.getcwd())
    else:
        print("Running in python jupyter notebook.")
        try:
            if not 'workbookDir' in globals():
                workbookDir = os.getcwd()
                print('workbookDir: ' + workbookDir)
                os.chdir(workbookDir)
        except:
            print("Problems with navigating to the workbook dir.")
#====================================================================================================#
# Imports

import re
import pandas as pd



#====================================================================================================#
# Global Functions
def beautiful_print(df):
    # Print the dataset in a well-organized format.
    with pd.option_context('display.max_rows'      , 20   , 
                           'display.min_rows'      , 20   , 
                           'display.max_columns'   , 5    , 
                           #"display.max_colwidth" , None , 
                           "display.width"         , None , 
                           "expand_frame_repr"     , True , 
                           "max_seq_items"         , None , ) :  # more options can be specified
        
        # Once the display.max_rows is exceeded, 
        # the display.min_rows options determines 
        # how many rows are shown in the truncated repr.

        print(df)
        
    return 


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

UniProt_Data_file    = "_DATABASE_uniprotkb_AND_reviewed_true_2023_10_23.tsv"
UniProt_Data_columns = ["Sequence", "Catalytic activity"]

output_folder        = Path("Processed_Dataset/")
output_folder.mkdir(parents = True, exist_ok = True)

output_file_dataset  = output_folder / "U00_Processed_Dataset.csv"
output_file_fasta    = output_folder / "U00_SwissProt_Sequences.fasta" 







#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

df_data = pd.read_csv(filepath_or_buffer = UniProt_Data_file     , 
                      sep                =  '\t'                 , 
                      usecols            =  UniProt_Data_columns ,
                      )


# Process one catalytic activity in one row.
def get_rxn_info(one_catalytic_activity):

    # Extracting the Reaction after `` CATALYTIC ACTIVITY: ... Reaction= ``
    reaction = re.search(r"CATALYTIC ACTIVITY:.*Reaction=(.*?);", one_catalytic_activity)
    reaction = reaction.group(1) if reaction else ""

    # Extracting the numbers after ChEBI:CHEBI:
    chebi_numbers = re.findall(r"ChEBI:CHEBI:(\d+)", one_catalytic_activity)
    chebi_numbers = chebi_numbers if chebi_numbers else [""]

    # Extracting the string starting with Rhea:
    rhea = re.findall(r"Rhea:RHEA(?:-COMP)?:(\d+)", one_catalytic_activity)
    rhea = rhea if rhea else [""]

    # Extracting the number after EC=
    ec = re.search(r"EC=(\d+\.\d+\.\d+\.\d+);", one_catalytic_activity)
    ec = ec.group(1) if ec else ""

    return [reaction, chebi_numbers, rhea, ec]



# Process one row of the dataframe.
def process_row(one_row):

    one_sequence    = one_row["Sequence"]
    activities_list = []

    if isinstance(one_row["Catalytic activity"], str):
        # Split by the pattern but without removing it
        activities = re.split(r"(?=CATALYTIC ACTIVITY)", one_row["Catalytic activity"])
        activities = [activity for activity in activities if activity]  # Removing any empty strings
        activities_data = [get_rxn_info(activity) for activity in activities]
        activities_list.extend(activities_data)

    return [one_sequence, activities_list]




# Explode the dataframe into a list of lists
all_data_rxn_info_list          = df_data.apply(process_row, axis = 1).tolist()
all_data_rxn_info_list_expanded = []

for one_pair_seqs_rxns in all_data_rxn_info_list:
    one_sequence, reactions_info_list = one_pair_seqs_rxns
    for one_reaction_info in reactions_info_list:
        all_data_rxn_info_list_expanded.append([one_sequence, ] + one_reaction_info)




# Print first 10 exploded data points to inspect
for i in range(min(6, len(all_data_rxn_info_list_expanded))):
    print(all_data_rxn_info_list_expanded[i])

print(len(all_data_rxn_info_list_expanded))

# Convert all_data_rxn_info_list_expanded to a dataframe
df_all_data_rxn_info_list_expanded = pd.DataFrame(all_data_rxn_info_list_expanded, 
                                                  columns = ["Sequence"    , 
                                                             "Reaction"    ,
                                                             "ChEBI_list"  ,
                                                             "Rhea_list"   ,
                                                             "EC"          ,
                                                             ])

beautiful_print(df_all_data_rxn_info_list_expanded)


# Save the dataframe to a csv file
df_all_data_rxn_info_list_expanded.to_csv(output_file_dataset, index = False)


# Print the number of unique sequences and reactions
seqs_list = [one_data[0]    for one_data in all_data_rxn_info_list_expanded]
rxn_list  = [one_data[1]    for one_data in all_data_rxn_info_list_expanded]

seqs_set_list = list(set(seqs_list))
rxn_set_list  = list(set(rxn_list ))

print("\n\nlen(seqs_set_list) : ", len(seqs_set_list))
print(    "len(rxn_set_list)  : ", len(rxn_set_list ))


# Write a fasta file with all the sequences
with open(output_file_fasta, "w") as f:
    for i in range(len(seqs_set_list)):
        f.write(">SwissProt_seq_" + str(i) + "\n")
        f.write(seqs_set_list[i] + "\n")



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#    .g8"""bgd `7MM       `7MM"""YMM  `7MM"""Yp, `7MMF'                    .M"""bgd `7MMM.     ,MMF'`7MMF'`7MMF'      `7MM"""YMM   .M"""bgd            #
#  .dP'     `M   MM         MM    `7    MM    Yb   MM         `MM.        ,MI    "Y   MMMb    dPMM    MM    MM          MM    `7  ,MI    "Y            #
#  dM'       `   MMpMMMb.   MM   d      MM    dP   MM           `Mb.      `MMb.       M YM   ,M MM    MM    MM          MM   d    `MMb.                #
#  MM            MM    MM   MMmmMM      MM"""bg.   MM    MMMMMMMMMMMMD      `YMMNq.   M  Mb  M' MM    MM    MM          MMmmMM      `YMMNq.            #
#  MM.           MM    MM   MM   Y  ,   MM    `Y   MM            ,M'      .     `MM   M  YM.P'  MM    MM    MM      ,   MM   Y  , .     `MM            #
#  `Mb.     ,'   MM    MM   MM     ,M   MM    ,9   MM          .M'        Mb     dM   M  `YM'   MM    MM    MM     ,M   MM     ,M Mb     dM            #
#    `"bmmmd'  .JMML  JMML.JMMmmmmMMM .JMMmmmd9  .JMML.                   P"Ybmmd"  .JML. `'  .JMML..JMML..JMMmmmmMMM .JMMmmmmMMM P"Ybmmd"             #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# 

import pickle
import requests
import xml.etree.ElementTree as ET


# Function to convert a list of ChEBI codes to a dictionary mapping ChEBI codes to SMILES strings
def chebi_list_to_smiles_dict(chebi_ids):
    """
    Convert a list of ChEBI codes to a dictionary mapping ChEBI codes to SMILES strings.
    
    :param chebi_list: List of ChEBI codes (as integers or strings without "CHEBI:" prefix).
    :return: Dictionary {ChEBI code: SMILES string}
    """
    base_url = "https://www.ebi.ac.uk/chebi/saveStructure.do?xml=true&chebiId="
    smiles_dict = {}
    
    for index, chebi_id in enumerate(chebi_ids, start=1):
        print()
        print("Fetching SMILES for CHEBI:", chebi_id, "|", index, "out of", len(chebi_ids))
        try:
            # Construct the URL for the API request
            url = f"{base_url}{chebi_id}"
            response = requests.get(url)
            
            if response.status_code == 200:
                # Parse the XML response
                root = ET.fromstring(response.content)
                smiles = None
                
                # Find the SMILES string in the XML
                for entry in root.findall('./ENTRY'):
                    smiles_element = entry.find('SMILES')
                    if smiles_element is not None:
                        smiles = smiles_element.text
                        break
                
                if smiles:
                    smiles_dict[chebi_id] = smiles
                    print(f"Found SMILES: {chebi_id} : {smiles}")
                else:
                    smiles_dict[chebi_id] = "No SMILES found"
                    print(f"No SMILES found for CHEBI:{chebi_id}")
            else:
                print(f"Failed to fetch data for CHEBI:{chebi_id}")
        except Exception as e:
            print(f"Error fetching SMILES for CHEBI:{chebi_id}: {e}")
    
    return smiles_dict



# Check if ChEBI2SMILES dictionary already exists

if (output_folder / "chebi_smiles_dict.pickle").exists():

    print("ChEBI2SMILES dictionary already exists. Loading from file... ")

    with open(output_folder / "chebi_smiles_dict.pickle", "rb") as f:
        chebi_smiles_dict = pickle.load(f)

    # Print the number of None values in the dictionary
    print("Number of None values in the dictionary: ", list(chebi_smiles_dict.values()).count(None))

else:

    print("ChEBI2SMILES dictionary does not exist. Fetching from ChEBI API... ")

    # Get a list of all ChEBI numbers from the ChEBI_list column in the dataframe.
    all_chebi_numbers = df_all_data_rxn_info_list_expanded["ChEBI_list"].tolist()

    # Flatten the list of lists into a list
    all_chebi_numbers = [chebi_number for chebi_numbers in all_chebi_numbers for chebi_number in chebi_numbers] 

    # Get a list of all unique ChEBI numbers
    all_chebi_numbers = list(set(all_chebi_numbers))

    # Print the number of unique ChEBI numbers
    print("len(all_chebi_numbers) : ", len(all_chebi_numbers))

    # Get a dictionary mapping ChEBI numbers to SMILES strings
    chebi_smiles_dict = chebi_list_to_smiles_dict(all_chebi_numbers)

    # Print the number of None values in the dictionary
    print("Number of None values in the dictionary: ", list(chebi_smiles_dict.values()).count(None))

    # Save the dictionary to a pickle file
    with open(output_folder / "chebi_smiles_dict.pickle", "wb") as f:
        pickle.dump(chebi_smiles_dict, f)




#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#    .g8"""bgd `7MM       `7MM"""YMM  `7MM"""Yp, `7MMF'                    .M"""bgd `7MMM.     ,MMF'`7MMF'`7MMF'      `7MM"""YMM   .M"""bgd            #
#  .dP'     `M   MM         MM    `7    MM    Yb   MM         `MM.        ,MI    "Y   MMMb    dPMM    MM    MM          MM    `7  ,MI    "Y            #
#  dM'       `   MMpMMMb.   MM   d      MM    dP   MM           `Mb.      `MMb.       M YM   ,M MM    MM    MM          MM   d    `MMb.                #
#  MM            MM    MM   MMmmMM      MM"""bg.   MM    MMMMMMMMMMMMD      `YMMNq.   M  Mb  M' MM    MM    MM          MMmmMM      `YMMNq.            #
#  MM.           MM    MM   MM   Y  ,   MM    `Y   MM            ,M'      .     `MM   M  YM.P'  MM    MM    MM      ,   MM   Y  , .     `MM            #
#  `Mb.     ,'   MM    MM   MM     ,M   MM    ,9   MM          .M'        Mb     dM   M  `YM'   MM    MM    MM     ,M   MM     ,M Mb     dM            #
#    `"bmmmd'  .JMML  JMML.JMMmmmmMMM .JMMmmmd9  .JMML.                   P"Ybmmd"  .JML. `'  .JMML..JMML..JMMmmmmMMM .JMMmmmmMMM P"Ybmmd"             #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# 









































































#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
###################################################################################################################
###################################################################################################################
#====================================================================================================#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#--------------------------------------------------#
#------------------------------

#                  __                  A             M       
#      `MM.        MM                 MMM            M       
#        `Mb.      MM   `MM.         MMMMM           M       
# MMMMMMMMMMMMD    MM     `Mb.     ,MA:M:AM.     `7M'M`MF'   
#         ,M'      MMMMMMMMMMMMD       M           VAMAV     
#       .M'                ,M'         M            VVV      
#                        .M'           M             V     




