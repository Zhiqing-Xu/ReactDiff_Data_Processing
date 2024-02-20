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
import ast
import csv
import copy
import pandas as pd


#====================================================================================================#
# Global Functions
def beautiful_print(df):
    # Print the dataset in a well-organized format.
    with pd.option_context('display.max_rows'      , 20   , 
                           'display.min_rows'      , 20   , 
                           'display.max_columns'   , 4    , 
                           #"display.max_colwidth" , None , 
                           "display.width"         , None , 
                           "expand_frame_repr"     , True , 
                           "max_seq_items"         , None , ) :  # more options can be specified
        
        # Once the display.max_rows is exceeded, 
        # the display.min_rows options determines 
        # how many rows are shown in the truncated repr. 

        print(df)
        
    return 



# 
Step_code = "U00B_"

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#   `7MMF'      ,M'   .g8""8q.        `7MM"""Mq.    db   MMP""MM""YMM`7MMF'  `7MMF'     ,gM""bg         `7MM"""Yb. `7MMF'`7MM"""Mq.                    #
#     MM        MV  .dP'    `YM.        MM   `MM.  ;MM:  P'   MM   `7  MM      MM       8MI  ,8           MM    `Yb. MM    MM   `MM.                   #
#     MM       AW   dM'      `MM        MM   ,M9  ,V^MM.      MM       MM      MM        WMp,"            MM     `Mb MM    MM   ,M9                    #
#     MM      ,M'   MM        MM        MMmmdM9  ,M  `MM      MM       MMmmmmmmMM       ,gPMN.  jM"'      MM      MM MM    MMmmdM9                     #
#     MM      MV    MM.      ,MP        MM       AbmmmqMA     MM       MM      MM      ,M.  YMp.M'        MM     ,MP MM    MM  YM.                     #
#     MM     AW     `Mb.    ,dP'        MM      A'     VML    MM       MM      MM      8Mp   ,MMp         MM    ,dP' MM    MM   `Mb.                   #
#   .JMML.  ,M'       `"bmmd"'        .JMML.  .AMA.   .AMMA..JMML.   .JMML.  .JMML.    `YMbmm'``MMm.    .JMMmmmdP' .JMML..JMML. .JMM.                  #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Input/Output files and folders

UniProt_Data_file    = "_DATABASE_uniprotkb_AND_reviewed_true_2023_10_23.tsv"
UniProt_Data_columns = ["Sequence", "Catalytic activity"]

output_folder        = Path("Processed_Dataset/")
output_folder.mkdir(parents = True, exist_ok = True)

output_file_fasta      = output_folder / (Step_code + "SP_Sequences.fasta"       )
output_file_dataset_1  = output_folder / (Step_code + "Preprocessed_Dataset.csv" )
output_file_dataset_2  = output_folder / (Step_code + "Processed_Dataset.csv"    )






#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#  `7MM"""Mq.`7MM"""Mq. `7MM"""YMM  `7MM"""Mq.`7MM"""Mq.   .g8""8q.     .g8"""bgd `7MM"""YMM   .M"""bgd  .M"""bgd `7MMF'`7MN.   `7MF' .g8"""bgd        #
#    MM   `MM. MM   `MM.  MM    `7    MM   `MM. MM   `MM..dP'    `YM. .dP'     `M   MM    `7  ,MI    "Y ,MI    "Y   MM    MMN.    M .dP'     `M        #
#    MM   ,M9  MM   ,M9   MM   d      MM   ,M9  MM   ,M9 dM'      `MM dM'       `   MM   d    `MMb.     `MMb.       MM    M YMb   M dM'       `        #
#    MMmmdM9   MMmmdM9    MMmmMM      MMmmdM9   MMmmdM9  MM        MM MM            MMmmMM      `YMMNq.   `YMMNq.   MM    M  `MN. M MM                 #
#    MM        MM  YM.    MM   Y  ,   MM        MM  YM.  MM.      ,MP MM.           MM   Y  , .     `MM .     `MM   MM    M   `MM.M MM.    `7MMF'      #
#    MM        MM   `Mb.  MM     ,M   MM        MM   `Mb.`Mb.    ,dP' `Mb.     ,'   MM     ,M Mb     dM Mb     dM   MM    M     YMM `Mb.     MM        #
#  .JMML.    .JMML. .JMM.JMMmmmmMMM .JMML.    .JMML. .JMM. `"bmmd"'     `"bmmmd'  .JMMmmmmMMM P"Ybmmd"  P"Ybmmd"  .JMML..JML.    YM   `"bmmmdPY        #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# 



# Process one catalytic activity in one row.
# Extract the reaction, ChEBI numbers, Rhea numbers, and EC number from the catalytic activity string.
def get_rxn_info(one_catalytic_activity):

    # Extracting the Reaction after `` CATALYTIC ACTIVITY: ... Reaction= ``
    reaction = re.search(r"CATALYTIC ACTIVITY:.*Reaction=(.*?);", one_catalytic_activity)
    reaction = reaction.group(1) if reaction else ""

    # Extracting the numbers after ChEBI:CHEBI:
    chebi_numbers = re.findall(r"ChEBI:CHEBI:(\d+)", one_catalytic_activity)
    chebi_numbers = chebi_numbers if chebi_numbers else []

    # Extracting the string starting with Rhea:
    #rhea = re.findall(r"Rhea:RHEA(?:-COMP)?:(\d+)", one_catalytic_activity)
    rhea = re.findall(r"Rhea:RHEA(?!-COMP):(\d+)", one_catalytic_activity)
    rhea = rhea if rhea else []



    # Extracting the number after EC=
    ec = re.search(r"EC=(\d+\.\d+\.\d+\.\d+);", one_catalytic_activity)
    ec = ec.group(1) if ec else ""

    return [reaction, chebi_numbers, rhea, ec]




# Process one row of the dataframe.
# Extract the sequence and all catalytic activities from one row.
def preprocess_row(one_row):

    one_sequence    = one_row["Sequence"]
    activities_list = []

    if isinstance(one_row["Catalytic activity"], str):
        # Split by the pattern but without removing it
        activities      = re.split(r"(?=CATALYTIC ACTIVITY)", one_row["Catalytic activity"])
        activities      = [activity for activity in activities if activity]  # Removing any empty strings
        activities_data = [get_rxn_info(activity) for activity in activities]
        activities_list.extend(activities_data)

    return [one_sequence, activities_list]



if output_file_dataset_1.exists() == False:

    # Load the dataframe.
    df_data = pd.read_csv(filepath_or_buffer = UniProt_Data_file     , 
                        sep                =  '\t'                 , 
                        usecols            =  UniProt_Data_columns , 
                        )


    # Explode the dataframe into a list of lists
    all_data_rxn_info_list          = df_data.apply(preprocess_row, axis = 1).tolist()
    all_data_rxn_info_list_expanded = []

    for one_pair_seqs_rxns in all_data_rxn_info_list:
        one_sequence, reactions_info_list = one_pair_seqs_rxns
        for one_reaction_info in reactions_info_list:
            all_data_rxn_info_list_expanded.append([one_sequence, ] + one_reaction_info)


    # Print first 10 exploded data points to inspect
    print("Printing first 10 exploded data points to inspect...")
    for i in range(min(4, len(all_data_rxn_info_list_expanded))):
        print(all_data_rxn_info_list_expanded[i])

    print("\nSize of the preprocessed dataset: ", len(all_data_rxn_info_list_expanded))



    # Convert all_data_rxn_info_list_expanded to a dataframe
    df_all_data_rxn_info_list_expanded = pd.DataFrame(all_data_rxn_info_list_expanded, 
                                                      columns = [ "Sequence"    , 
                                                                  "Reaction"    , 
                                                                  "ChEBI_list"  , 
                                                                  "Rhea_list"   , 
                                                                  "EC"          , ])
    

    # Save the dataframe to a csv file
    df_all_data_rxn_info_list_expanded.to_csv(output_file_dataset_1, index = False)




else:
    


    # Define a function to convert string lists to actual lists
    def string_to_list(string):
        try:
            # Safely evaluate the string as a Python literal (list)
            return ast.literal_eval(string)
        except (ValueError, SyntaxError):
            # Return an empty list in case of an error (e.g., if the string is not a valid Python literal)
            return []

    # Specify the columns that need to be converted from string to list
    converters = {
        'ChEBI_list': string_to_list,
        'Rhea_list': string_to_list
    }

    # Reading the CSV file with the converters applied to the specified columns
    df_all_data_rxn_info_list_expanded = pd.read_csv(output_file_dataset_1, converters=converters)


beautiful_print(df_all_data_rxn_info_list_expanded[["Reaction"    , 
                                                    "ChEBI_list"  , 
                                                    "Rhea_list"   , 
                                                    "EC"          , ]])



all_data_rxn_info_list_expanded = df_all_data_rxn_info_list_expanded.values.tolist()






#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#  `7MMF'     A     `7MF'`7MM"""Mq. `7MMF'MMP""MM""YMM `7MM"""YMM        .M"""bgd `7MM"""Mq.     `7MM"""YMM  db       .M"""bgd MMP""MM""YMM  db        #
#    `MA     ,MA     ,V    MM   `MM.  MM  P'   MM   `7   MM    `7       ,MI    "Y   MM   `MM.      MM    `7 ;MM:     ,MI    "Y P'   MM   `7 ;MM:       #
#     VM:   ,VVM:   ,V     MM   ,M9   MM       MM        MM   d         `MMb.       MM   ,M9       MM   d  ,V^MM.    `MMb.          MM     ,V^MM.      #
#      MM.  M' MM.  M'     MMmmdM9    MM       MM        MMmmMM           `YMMNq.   MMmmdM9        MM""MM ,M  `MM      `YMMNq.      MM    ,M  `MM      #
#      `MM A'  `MM A'      MM  YM.    MM       MM        MM   Y  ,      .     `MM   MM             MM   Y AbmmmqMA   .     `MM      MM    AbmmmqMA     #
#       :MM;    :MM;       MM   `Mb.  MM       MM        MM     ,M      Mb     dM   MM             MM    A'     VML  Mb     dM      MM   A'     VML    #
#        VF      VF      .JMML. .JMM.JMML.   .JMML.    .JMMmmmmMMM      P"Ybmmd"  .JMML.         .JMML..AMA.   .AMMA.P"Ybmmd"     .JMML.AMA.   .AMMA.  #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Print the number of unique sequences and reactions
seqs_list = [one_data[0] for one_data in all_data_rxn_info_list_expanded]
rxn_list  = [one_data[1] for one_data in all_data_rxn_info_list_expanded]

seqs_set_list = list(set(seqs_list))
rxn_set_list  = list(set(rxn_list ))


print("\n\nInformation about the preprocessed dataset:")
print(    "len(seqs_set_list) : ", len(seqs_set_list ))
print(    "len(rxn_set_list)  : ", len(rxn_set_list  ))


# Write a fasta file with all the sequences
with open(output_file_fasta, "w") as f:
    for i in range(len(seqs_set_list)):
        f.write(">SP_" + str(i) + "\n")
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
# Get chebi_smiles_dict

import pickle
import requests
import xml.etree.ElementTree as ET

output_chebi_smiles_dict_path = output_folder / (Step_code + "chebi_smiles_dict.p")



# Function to convert a list of ChEBI codes to a dict mapping ChEBI codes to SMILES strings.
def chebi_list_to_smiles_dict(chebi_ids, smiles_dict = None, use_raw_text = True):

    """
    Convert a list of ChEBI codes to a dictionary mapping ChEBI codes to SMILES strings.
        :param chebi_list: - List of ChEBI codes (as integers or strings without "CHEBI:" prefix).
        :return:           - Dictionary {ChEBI code: SMILES string}
    """

    base_url           = "https://www.ebi.ac.uk/chebi/saveStructure.do?xml=true&chebiId="
    smiles_dict_exists = True


    # If a dictionary is not provided, create a new one
    if smiles_dict is None:
        smiles_dict        = {}
        smiles_dict_exists = False


    for index, chebi_id in enumerate(chebi_ids, start = 1):

        if smiles_dict_exists == False:
            print("Fetching SMILES for CHEBI:", chebi_id, "|", index, "out of", len(chebi_ids))

        # If the ChEBI code is already in the dictionary, skip it.
        if chebi_id in smiles_dict:
            continue
        

        print(f"Fetching SMILES not in dict for CHEBI: {chebi_id}")

        try:
            # Construct the URL for the API request
            url      = f"{base_url}{chebi_id}"
            response = requests.get(url)
            
            if response.status_code == 200:

                # Extract SMILES from raw text.
                if use_raw_text == True:

                    # Extract SMILES from raw text
                    content = response.text
                    smiles_matches = re.findall(r'<SMILES>(.*?)</SMILES>', content)

                    if smiles_matches:
                        if len(smiles_matches) > 1:
                            print("Warning: Multiple SMILES entries found. Using the first one.")
                        
                        smiles = smiles_matches[0]
                        smiles_dict[chebi_id] = smiles
                        print(f"Found SMILES: {chebi_id} : {smiles}")

                    else:
                        print(f"No SMILES found for CHEBI:{chebi_id}")
                        smiles_dict[chebi_id] = "No SMILES found"

                # Extract SMILES from XML.
                if use_raw_text == False:

                    # Preprocess the response content to escape special characters
                    preprocessed_content = re.sub(r"<(?!/?[a-zA-Z]+(?: [a-zA-Z]+=\"[^\"]*\")*/?>)", "&lt;", response.text)
                    preprocessed_content = re.sub(r"(?<!<)/>", "&gt;", preprocessed_content)

                    root = ET.fromstring(preprocessed_content)
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
            # If an error occurs, print the error message
            print(f"Error fetching SMILES for CHEBI:{chebi_id}: {e}")

            # Add the ChEBI code to the dictionary with a value of None
            smiles_dict[chebi_id] = None


    return smiles_dict



# Get a list of all ChEBI numbers from the ChEBI_list column in the dataframe.
all_chebi_numbers = df_all_data_rxn_info_list_expanded["ChEBI_list"].tolist()

# Flatten the list of lists into a list
all_chebi_numbers = [chebi_number for chebi_numbers in all_chebi_numbers for chebi_number in chebi_numbers] 

# Get a list of all unique ChEBI numbers
all_chebi_numbers = list(set(all_chebi_numbers))

# Print the number of unique ChEBI numbers
print("\n\nInformation about the processing dataset:")
print("len(all_chebi_numbers)        : ", len(all_chebi_numbers), " (number of unique ChEBI numbers)")




# For RHEA, get a list of lists, all_RHEA_numbers, and flatten it into a list
all_RHEA_numbers           = df_all_data_rxn_info_list_expanded["Rhea_list"].tolist()
all_RHEA_numbers_flattened = [rhea_number for rhea_numbers in all_RHEA_numbers for rhea_number in rhea_numbers]
all_unique_RHEA_list       = []
for rhea_list in all_RHEA_numbers:
    all_unique_RHEA_list.append(rhea_list) if rhea_list not in all_unique_RHEA_list else None
print("Number of all RHEA numbers     : ", len(set(all_RHEA_numbers_flattened))                     )
print("Number of all unique RHEA list : ", len(all_unique_RHEA_list)                                )
print("Number of non-empty RHEA list  : ", sum([1 for rhea_list in all_RHEA_numbers if rhea_list])  )





# Check if ChEBI-to-SMILES dictionary already exists
if (output_chebi_smiles_dict_path).exists():

    print("\n\nChEBI-to-SMILES dictionary already exists. Loading from file and update... ")

    with open(output_chebi_smiles_dict_path, "rb") as f:
        chebi_smiles_dict = pickle.load(f)

    # Update the dictionary.
    chebi_smiles_dict = chebi_list_to_smiles_dict(all_chebi_numbers, chebi_smiles_dict)

    # get all items in chebi_smiles_dict with None values or "No SMILES found"
    none_items_1 = {k: v for k, v in chebi_smiles_dict.items() if v == None}
    print("Number of **None** values in the dictionary: ", len(none_items_1))
    none_items_2 = {k: v for k, v in chebi_smiles_dict.items() if v == "No SMILES found"}
    print("Number of **No SMILES found** values in the dictionary: ", len(none_items_2))

    # Replace the values of None with "No SMILES found"
    chebi_smiles_dict = {k: "No SMILES found" if v == None else v for k, v in chebi_smiles_dict.items()}



    # Print the number of None values in the dictionary
    print("Number of **None** values in the dictionary: ", list(chebi_smiles_dict.values()).count(None))


else:

    print("\n\nChEBI-to-SMILES dictionary does not exist. Fetching from ChEBI API... ")

    # Get a dictionary mapping ChEBI numbers to SMILES strings
    chebi_smiles_dict = chebi_list_to_smiles_dict(all_chebi_numbers)

    # Print the number of None values in the dictionary
    print("Number of **None** values in the dictionary: ", list(chebi_smiles_dict.values()).count(None))



# Save the dictionary to a pickle file
with open(output_chebi_smiles_dict_path, "wb") as f:
    pickle.dump(chebi_smiles_dict, f)


# Print out the first 5 entries in the dictionary
print("Print First 10 entries in the chebi_smiles_dict:")
for key, value in list(chebi_smiles_dict.items())[:10]:
    print(key      , ":" , value       ) # Print the key and value
    print(type(key), ":" , type(value) ) # Print out the type of key and value





#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#   .g8"""bgd`7MM"""YMM MMP""MM""YMM  `7MM"""Mq. `7MMF'  `7MMF'7MM"""YMM       db                     .g8"""bgd`7MMF'  `7MMF'7MM"""YMM  7MM"""Yp,`7MMF'#
# .dP'     `M  MM    `7 P'   MM   `7    MM   `MM.  MM      MM   MM    `7      ;MM:         `MM.     .dP'     `M  MM      MM   MM    `7   MM    Yb  MM  #
# dM'       `  MM   d        MM         MM   ,M9   MM      MM   MM   d       ,V^MM.          `Mb.   dM'       `  MM      MM   MM   d     MM    dP  MM  #
# MM           MMmmMM        MM         MMmmdM9    MMmmmmmmMM   MMmmMM      ,M  `MM    MMMMMMMMMMMD MM           MMmmmmmmMM   MMmmMM     MM"""bg.  MM  #
# MM.    `7MMF'MM   Y  ,     MM         MM  YM.    MM      MM   MM   Y  ,   AbmmmqMA          ,M'   MM.          MM      MM   MM   Y  ,  MM    `Y  MM  #
# `Mb.     MM  MM     ,M     MM         MM   `Mb.  MM      MM   MM     ,M  A'     VML       .M'     `Mb.     ,'  MM      MM   MM     ,M  MM    ,9  MM  #
#   `"bmmmdPY.JMMmmmmMMM   .JMML.     .JMML. .JMM.JMML.  .JMML.JMMmmmmMMM.AMA.   .AMMA.               `"bmmmd' .JMML.  .JMML.JMMmmmmMMM JMMmmmd9 .JMML.#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Read "./_DATABASE_RHEA.tsv" and get the RHEA-to-SMILES dictionary


# Define the path to your TSV file
input_RHEA_file_path  = './_DATABASE_RHEA.tsv'
output_RHEA_info_dict = output_folder / (Step_code + "RHEA_info_dict.p")


# Read the equation and modify the coefficients and suffixes
def preprocess_rxn_eqn(equation):

    print("\nProcessing :", equation)

    # Split equation into reactants and products
    reactants, products = equation.split(" = ")
    reactants_list      = reactants.split(" + ")
    products_list       = products.split(" + ")

    # Initialize lists to hold identified coefficients and suffixes
    identified_coefficients = []
    identified_suffixes     = []

    # Step 1: Identify coefficients with n
    coefficient_pattern = re.compile(r"^\(?(n[\+\-]?\d*)\)?\s")
    #coefficient_pattern = re.compile(r"^\((n[\+\-]\d*)\) |^n ")
    for compound in reactants_list + products_list:
        coefficient_match = coefficient_pattern.match(compound)
        #print(coefficient_match)
        if coefficient_match:
            identified_coefficients.append(coefficient_match.group().strip())


    # Step 2: Identify suffixes with n that follow square brackets
    suffix_pattern      = re.compile(r"\[\S+\]\((n[\+\-]?\d*)\)")

    for compound in reactants_list + products_list:
        for match in suffix_pattern.finditer(compound):
            identified_suffixes.append(match.group(1))


    # Reporting identified coefficients and suffixes
    n_coefficients_exists = False
    n_suffixes_exists     = False

    if identified_coefficients:
        print("    Identified coefficients with 'n':", identified_coefficients)
        n_coefficients_exists = True
    else:
        print("    No coefficients contain 'n'")

    if identified_suffixes:
        print("    Identified suffix with 'n':", identified_suffixes)
        n_suffixes_exists = True
    else:
        print("    No suffix contain 'n'")

    # If no coefficients or suffixes contain 'n', return the original equation
    if not n_coefficients_exists and not n_suffixes_exists:
        return equation

    # Step 3: Determine n based on identified coefficients and suffixes
    min_n_value = 1
    for expression in identified_coefficients + identified_suffixes:
        if 'n' in expression:
            value = eval(expression.replace('n', '1'))
            min_n_value = min(min_n_value, value)

    n = 1 - min_n_value + 1  # Adjust n so the smallest (n-m) or (n) is 1
    print("    Determined n based on coefficients and suffixes:", n)


    # Step 4: Replace coefficients and suffixes with the computed actual number
    def replace_expression_suf(match):
        expr = match.group(1)
        return str(eval(expr.replace('n', str(n))))

    def update_compound_suf(compound):
        # This pattern captures the entire compound including its suffix to be modified.
        # Adjust this pattern if necessary to match the specific format of your compounds.
        compound_pattern = re.compile(r"(\[\S+\])\((n[\+\-]?\d*)\)")

        # Function to replace and move suffix to the front
        def replace_and_move(match):
            base_part = match.group(1)  # The compound part without the suffix
            expr = match.group(2)       # The expression part within the suffix
            evaluated_suffix = eval(expr.replace('n', str(n)))  # Evaluate the suffix expression
            return str(evaluated_suffix) + " " + base_part  # Prepend the evaluated suffix as a coefficient

        # Apply the replacement and moving process
        updated_compound = compound_pattern.sub(replace_and_move, compound)

        return updated_compound


    def replace_expression_coef(match):
        expr = match.group(1)
        try:
            # Evaluate the expression, replacing 'n' with its computed value
            return str(eval(expr.replace('n', str(n))))
        except Exception as e:
            print(f"Error evaluating expression: {expr}. Error: {e}")
            import pdb; pdb.set_trace()


    def update_compound_coef(compound):
        # Replace coefficients at the beginning
        compound = coefficient_pattern.sub(lambda m: replace_expression_coef(m) + " ", compound, 1)
        
        # Handle suffixes
        # Split the compound to separate the part within square brackets and the suffix
        parts = compound.split("[")
        if len(parts) > 1:
            base, suffix_part = parts[0], parts[1]
            # Replace the suffix expression
            suffix = suffix_pattern.sub(lambda m: replace_expression_coef(m), suffix_part)
            # Reassemble the compound
            compound = base + "[" + suffix
        return compound


    def decide_and_update(compound):
        # Check if "(n" is likely a coefficient (at the start or not followed by ']')
        if compound.startswith("n ") or "(n" in compound and not re.search(r"\[\S*\]\(n", compound):
            return update_compound_coef(compound)
        else:
            # If "(n" follows square brackets or doesn't start the compound, treat as suffix
            return update_compound_suf(compound)

    updated_reactants = [decide_and_update(comp) for comp in reactants_list]
    print("    Updated reactants:", updated_reactants)

    updated_products = [decide_and_update(comp) for comp in products_list]
    print("    Updated products:", updated_products)


    # Construct the new equation
    new_equation = " + ".join(updated_reactants) + " = " + " + ".join(updated_products)
    print("Processed  :", new_equation)
    return new_equation



# Function to process Equation column and extract coefficients
def process_rxn_eqn_basic_info(equation):
    # Split equation into reactants and products
    reactants, products = equation  .split(" = ")
    reactants_list      = reactants .split(" + ")
    products_list       = products  .split(" + ")

    # Process reactants and products to extract or assign coefficients
    def extract_coefficients(compound_list):
        coeffs = []
        for compound in compound_list:
            parts = compound.split(' ')
            # Check if the first part is purely numeric to determine if it's a coefficient
            if parts[0].isdigit():
                coeffs.append(int(parts[0]))
            else:
                coeffs.append(1)
        return coeffs


    reactant_coeffs = extract_coefficients(reactants_list)
    product_coeffs  = extract_coefficients(products_list )

    return reactant_coeffs, product_coeffs



# Function to process Equation column and further process product chebi_ids
def process_rxn_eqn_additional_info(equation):
    # Split the reaction into reactants and products
    reactants_str, products_str = equation     .split(" = ")
    reactants_list              = reactants_str.split(" + ")
    products_list               = products_str .split(" + ")
    
    reactants_list_abnormal = []
    products_list_abnormal  = []
    for c in reactants_list:
        if ( c.find("(in)" )  != -1 or 
             c.find("(out)")  != -1 or 
             c.find("(n)"  )  != -1 or 
             c.find("(n+1)")  != -1 or 
             c.find("(n-1)")  != -1 
            ):
            reactants_list_abnormal.append(1)
        else:
            reactants_list_abnormal.append(0)

    for c in products_list:
        if ( c.find("(in)" )  != -1 or 
             c.find("(out)")  != -1 or 
             c.find("(n)"  )  != -1 or 
             c.find("(n+1)")  != -1 or 
             c.find("(n-1)")  != -1 
            ):
            products_list_abnormal.append(1)
        else:
            products_list_abnormal.append(0)


    # Apply the function to reactants and products
    reactants_marks_in_out = copy.deepcopy(reactants_list_abnormal)
    products_marks_in_out  = copy.deepcopy(products_list_abnormal )
    

    # Prepare the list for mapping products to reactants
    product_to_reactant_map = []
    
    # Remove (in) or (out) and map products to reactants
    for product in products_list:
        cleaned_product = product.replace('(in)' , '') \
                                 .replace('(out)', '') \
                                 .replace('(n)'  , '') \
                                 .replace('(n-1)', '') \
                                 .replace('(n+1)', '').strip()
        if cleaned_product.split(' ')[0].isdigit():
            cleaned_product = cleaned_product.split(' ')[1]
        
        match_found = False
        for i, reactant in enumerate(reactants_list):
            cleaned_reactant = reactant.replace('(in)' , '') \
                                       .replace('(out)', '') \
                                       .replace('(n)'  , '') \
                                       .replace('(n-1)', '') \
                                       .replace('(n+1)', '').strip()
            if cleaned_reactant.split(' ')[0].isdigit():
                cleaned_reactant = cleaned_reactant.split(' ')[1]

            if cleaned_product == cleaned_reactant:
                product_to_reactant_map.append(i)
                match_found = True
                break

        if not match_found:
            if '(in)'  in product or \
               '(out)' in product or \
               '(n)'   in product or \
               '(n-1)' in product or \
               '(n+1)' in product :
                product_to_reactant_map.append(-2)
            else:
                product_to_reactant_map.append(-1)
    
    return [[reactants_marks_in_out, products_marks_in_out], product_to_reactant_map]



# Function to process ChEBI identifier column
def process_chebi_identifiers(chebi_str):
    chebi_pairs = chebi_str.split(";")
    chebi_ids = [int(chebi.split(":")[1]) for chebi in chebi_pairs]
    return chebi_ids


# Main function to process the TSV file and create the desired output
def Get_RHEA_info_dict(file_path):
    # Read the TSV file
    df = pd.read_csv(file_path, sep='\t')

    # Initialize an empty dict for the output
    output_dict = {}

    # Iterate through the dataframe
    len_rhea_file = len(df)
    for index, row in df.iterrows():
        print("\nProcessing row", index, "out of", len_rhea_file)

        preprocessed_eqn = preprocess_rxn_eqn(row['Equation'])

        # Process the equation to get reactant and product coefficients
        reactant_coeffs, product_coeffs = process_rxn_eqn_basic_info(preprocessed_eqn)
        # Process the ChEBI identifiers
        chebi_ids = process_chebi_identifiers(row['ChEBI identifier'])

        # Modify the product ChEBI IDs to account for (in) and (out) compounds
        rctt_chebi_ids = chebi_ids[ : len(reactant_coeffs)   ]
        prod_chebi_ids = chebi_ids[   len(reactant_coeffs) : ]
        prod_chebi_ids_it = iter(prod_chebi_ids)

        equation_additional_info = process_rxn_eqn_additional_info(preprocessed_eqn)
        product_to_reactant_map  = equation_additional_info[1]

        prod_chebi_ids_modified = []
        for one_info in product_to_reactant_map:
            if one_info   == -1:
                try:
                    prod_chebi_ids_modified.append(next(prod_chebi_ids_it))
                except:
                    print("Error: No more ChEBI IDs to iterate through.")
                    print(preprocessed_eqn)
                    print(equation_additional_info)
                    print(prod_chebi_ids)
                    # import pdb; pdb.set_trace() 
                    
            elif one_info == -2:
                print(preprocessed_eqn)
                prod_chebi_ids_modified.append(next(prod_chebi_ids_it))

            else:
                prod_chebi_ids_modified.append(rctt_chebi_ids[one_info])



        # ZX_!!! : This is a temporary fix. Need to find a better solution.
        if len(prod_chebi_ids_modified) != len(product_coeffs):
            print("Error: Length of modified ChEBI IDs does not match the length of product coefficients.")
            if reactant_coeffs == product_coeffs:
                prod_chebi_ids_modified = copy.deepcopy(rctt_chebi_ids)



        # Assign values to the output dict
        output_dict[row['Reaction identifier']] = [preprocessed_eqn                          , 
                                                   (reactant_coeffs, product_coeffs)         , 
                                                   rctt_chebi_ids + prod_chebi_ids_modified  , # modified chebi_ids
                                                   row['Enzyme class']                       , 
                                                  ]

    return output_dict


# Check if the dictionary already exists
if (output_RHEA_info_dict).exists() == False:
    print("\n\nRHEA info dictionary (RHEA_info_dict) does not exist. Fetching from RHEA TSV file... ")
    RHEA_info_dict = Get_RHEA_info_dict(input_RHEA_file_path)
    # Save the dictionary to a pickle file
    with open(output_RHEA_info_dict, "wb") as f:
        pickle.dump(RHEA_info_dict, f)
else:
    print("\n\nRHEA info dictionary (RHEA_info_dict) already exists. Loading from file... ")
    with open(output_RHEA_info_dict, "rb") as f:
        RHEA_info_dict = pickle.load(f)



# Print out the first 5 entries in the dictionary
print("Print First 10 entries in the RHEA_info_dict:")
for key, value in list(RHEA_info_dict.items())[:10]:
    print(key, ":", value)  # Print the key and value







#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#    .g8"""bgd `7MM"""YMM MMP""MM""YMM    `7MM"""Mq. `YMM'   `MP'`7MN.   `7MF'     .M"""bgd `7MMM.     ,MMF'`7MMF'`7MMF'      `7MM"""YMM   .M"""bgd    #
#  .dP'     `M   MM    `7 P'   MM   `7      MM   `MM.  VMb.  ,P    MMN.    M      ,MI    "Y   MMMb    dPMM    MM    MM          MM    `7  ,MI    "Y    #
#  dM'       `   MM   d        MM           MM   ,M9    `MM.M'     M YMb   M      `MMb.       M YM   ,M MM    MM    MM          MM   d    `MMb.        #
#  MM            MMmmMM        MM           MMmmdM9       MMb      M  `MN. M        `YMMNq.   M  Mb  M' MM    MM    MM          MMmmMM      `YMMNq.    #
#  MM.    `7MMF' MM   Y  ,     MM           MM  YM.     ,M'`Mb.    M   `MM.M      .     `MM   M  YM.P'  MM    MM    MM      ,   MM   Y  , .     `MM    #
#  `Mb.     MM   MM     ,M     MM           MM   `Mb.  ,P   `MM.   M     YMM      Mb     dM   M  `YM'   MM    MM    MM     ,M   MM     ,M Mb     dM    #
#    `"bmmmdPY .JMMmmmmMMM   .JMML.       .JMML. .JMM.MM:.  .:MMa.JML.    YM      P"Ybmmd"  .JML. `'  .JMML..JMML..JMMmmmmMMM .JMMmmmmMMM P"Ybmmd"     #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# 


all_Sequences_list = df_all_data_rxn_info_list_expanded["Sequence"  ].tolist()
all_Reactions_list = df_all_data_rxn_info_list_expanded["Reaction"  ].tolist()
all_ChEBI_list     = df_all_data_rxn_info_list_expanded["ChEBI_list"].tolist()
all_Rhea_list      = df_all_data_rxn_info_list_expanded["Rhea_list" ].tolist()
all_EC_list        = df_all_data_rxn_info_list_expanded["EC"        ].tolist()


assert print(len(all_Reactions_list)) == print(len(all_Rhea_list)) == print(len(all_EC_list)), "Size of the lists are not equal."

final_data_list    = []
unfound_RHEA_list  = []
unfound_ChEBI_list = []

# Print the first 5 entries in the lists
len_dataset = len(all_Sequences_list)
for idx, (seqn, rctn, cheb, rhea, ecno) in enumerate(zip(all_Sequences_list, all_Reactions_list, all_ChEBI_list, all_Rhea_list, all_EC_list)):
    print("\nProcessing", idx, "out of", len_dataset, "...")

    # if Rhea_list is empty, skip
    if not rhea:
        print("    Empty Rhea list. Skipping...")
        continue

    one_datapoint = [None, ] * 3
    for one_rhea in rhea:

        one_rhea_valid = True

        if one_rhea:
            one_rhea_id = "RHEA:" + one_rhea
            print("    Processing", one_rhea_id)

            if one_rhea_id not in RHEA_info_dict:
                unfound_RHEA_list.append(one_rhea_id)
                print(f"    Unfound RHEA: {one_rhea_id}")
                continue


            if one_rhea_id in RHEA_info_dict:
                one_datapoint_rctt_coef   = RHEA_info_dict[one_rhea_id][1][0]                                # e.g., [1, 1, 1]
                one_datapoint_prod_coef   = RHEA_info_dict[one_rhea_id][1][1]                                # e.g., [1, 1, 1]
                one_datapoint_rctt_cheb   = RHEA_info_dict[one_rhea_id][2][ : len(one_datapoint_rctt_coef) ] # e.g., [15377, 16452, 57392]
                one_datapoint_prod_cheb   = RHEA_info_dict[one_rhea_id][2][len(one_datapoint_rctt_coef) :  ] # e.g., [58853, 57287, 15378]

                one_datapoint_rctt_smls   = []
                one_datapoint_prod_smls   = []

                #--------------------------------------------------#
                # 
                for one_chebi in one_datapoint_rctt_cheb:
                    if str(one_chebi) in chebi_smiles_dict:
                        one_smiles = chebi_smiles_dict[str(one_chebi)]
                        if one_smiles != "No SMILES found":
                            one_datapoint_rctt_smls.append(one_smiles)
                        else:
                            unfound_ChEBI_list.append(one_chebi)
                            print(f"    Unknown SMILES for ChEBI: {one_chebi}", "Skipping...")
                            one_rhea_valid = False
                            break
                    else:
                        unfound_ChEBI_list.append(one_chebi)
                        print(f"    Unfound ChEBI: {one_chebi}", "Skipping...")
                        one_rhea_valid = False
                        break


                for one_chebi in one_datapoint_prod_cheb:
                    if str(one_chebi) in chebi_smiles_dict:
                        one_smiles = chebi_smiles_dict[str(one_chebi)]
                        if one_smiles != "No SMILES found":
                            one_datapoint_prod_smls.append(one_smiles)
                        else:
                            unfound_ChEBI_list.append(one_chebi)
                            print(f"    Unknown SMILES for ChEBI: {one_chebi}", "Skipping...")
                            one_rhea_valid = False
                            break
                    else:
                        unfound_ChEBI_list.append(one_chebi)
                        print(f"    Unfound ChEBI: {one_chebi}", "Skipping...")
                        one_rhea_valid = False
                        break


                if len(one_datapoint_rctt_coef) + len(one_datapoint_prod_coef) != len(one_datapoint_rctt_smls) + len(one_datapoint_prod_smls):
                    print(">>> Error Found: Number of coef's doesnt match with number of ChEBI's. Skipping...")
                    print(f"rctn: {rctn}")
                    print(f"cheb: {cheb}")
                    print("rctn:", RHEA_info_dict[one_rhea_id][0])
                    print("coef:", RHEA_info_dict[one_rhea_id][1])
                    print("cheb:", RHEA_info_dict[one_rhea_id][2])



                    one_rhea_valid = False



                if one_rhea_valid == False:
                    continue


                #--------------------------------------------------#
                # 
                one_datapoint_rctt_smls_w_coef = []
                one_datapoint_prod_smls_w_coef = []


                try:

                    for idx, one_smls in enumerate(one_datapoint_rctt_smls):
                        for _ in range(one_datapoint_rctt_coef[idx]):
                            one_datapoint_rctt_smls_w_coef.append(one_smls)
                    
                    for idx, one_smls in enumerate(one_datapoint_prod_smls):
                        for _ in range(one_datapoint_prod_coef[idx]):
                            one_datapoint_prod_smls_w_coef.append(one_smls)

                except Exception as e:

                    print(f">>> Error: {e}")
                   #print(f"seqn: {seqn}")
                    print(f"rctn: {rctn}")
                    print(f"cheb: {cheb}")
                    print(f"rhea: {rhea}")
                    print(f"ecno: {ecno}")

                    print(f"one_datapoint_rctt_smls: {one_datapoint_rctt_smls}")
                    print(f"one_datapoint_prod_smls: {one_datapoint_prod_smls}")
                    print(f"one_datapoint_rctt_coef: {one_datapoint_rctt_coef}")
                    print(f"one_datapoint_prod_coef: {one_datapoint_prod_coef}")


                    print(">>> Error Found: Unknown Error. ")
                    import pdb; pdb.set_trace()

                one_datapoint_smls_string = ".".join(one_datapoint_rctt_smls_w_coef) + ">>" + ".".join(one_datapoint_prod_smls_w_coef)
                print("    Processed SMILES reaction: ", "\n    ", one_datapoint_smls_string)

                one_datapoint[1] = one_datapoint_smls_string
                one_datapoint[0] = seqn
                one_datapoint[2] = ecno
                final_data_list.append(one_datapoint)



print("Number of ChEBI numbers not found in the ChEBI-to-SMILES dict : ", len(unfound_ChEBI_list) )
print("Number of RHEA numbers not found in the RHEA dict             : ", len(unfound_RHEA_list ) )




# Convert final_data_list to a dataframe
final_data_df = pd.DataFrame(final_data_list, columns = ["Sequence", "Reaction", "EC_number", ])
final_data_df.to_csv(output_file_dataset_2, index = True)







































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




