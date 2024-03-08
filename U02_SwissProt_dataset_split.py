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
import copy
import pandas as pd
import pickle

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

# 
Step_code = "U02_"

# 
output_folder        = Path("Processed_Dataset/")

SwissProt_split_dict_file = output_folder / ("U01_SwissProt_split_dict.p")
Processed_SwissProt_file  = output_folder / ("U00B_Processed_Dataset.csv")
SwissProt_fasta_file      = output_folder / ("U00B_SP_Sequences.fasta")


# Read the fasta file into a dictionary that maps the sequence id to the sequence.
SwissProt_dict = {}
with open(SwissProt_fasta_file, 'r') as file:
    for line in file:
        if line.startswith('>'):
            seq_id = line.strip().split(' ')[0][1:]
            SwissProt_dict[seq_id] = ''
        else:
            SwissProt_dict[seq_id] += line.strip()


# Read the processed dataset into a dataframe
# The csv file contains three columns: "Sequence", "Reaction", and "EC_number"
Processed_SwissProt_df = pd.read_csv(Processed_SwissProt_file, index_col = 0)
all_reaction_list     = Processed_SwissProt_df['Reaction'].tolist()
all_reaction_set_list = list(set(all_reaction_list))
print("Number of unique reactions: ", len(all_reaction_set_list))


# Get mapped_reaction
reaction_to_mapped_reaction_dict_path = output_folder / (Step_code + "reaction_to_mapped_reaction_dict.p")


if path.exists(reaction_to_mapped_reaction_dict_path):
    reaction_to_mapped_reaction_dict = pickle.load(open(reaction_to_mapped_reaction_dict_path, "rb"))

else:

    from rxnmapper import RXNMapper
    rxn_mapper = RXNMapper()

    reaction_modified_dict = {}
    reaction_modified_list = []


    for idx, one_rxn in enumerate(all_reaction_set_list):
        one_rxn_original = copy.deepcopy(one_rxn)
        one_rxn = one_rxn.replace("-*" , "C")
        one_rxn = one_rxn.replace("[*]", "C")
        one_rxn = one_rxn.replace("*"  , "C")

        reaction_modified_dict[one_rxn_original] = one_rxn
        reaction_modified_list.append(one_rxn)

        # print("\nMapping reaction: ", idx, " out of ", len(all_reaction_set_list))
        # try:
        #     rxn_mapper_result = rxn_mapper.get_attention_guided_atom_maps([one_rxn, ])
        #     mapped_rxn        = rxn_mapper_result[0]['mapped_rxn']
        #     reaction_to_mapped_reaction[one_rxn] = mapped_rxn
        #     print("Successfully mapped reaction: ", one_rxn)
        # except:
        #     print("Error mapping reaction: ", one_rxn)
        #     reaction_to_mapped_reaction[one_rxn] = None

    from rxnmapper import BatchedMapper
    rxn_mapper = BatchedMapper(batch_size = 32)

    rxn_map_results = list(rxn_mapper.map_reactions(reaction_modified_list))

    reaction_to_mapped_reaction_dict = {}
    for idx, (one_rxn_original, one_rxn_mapped) in enumerate(zip(all_reaction_set_list, rxn_map_results)):
        if one_rxn_mapped == ">>":
            print("Error mapping reaction: ", one_rxn_original)
            reaction_to_mapped_reaction_dict[one_rxn_original] = None
        reaction_to_mapped_reaction_dict[one_rxn_original] = one_rxn_mapped

    # Save the dictionary to a pickle file
    pickle.dump(reaction_to_mapped_reaction_dict, open(reaction_to_mapped_reaction_dict_path, "wb"))



# Read the split dictionary from the pickle file
# There are two keys in the dict SP_seqs_similar_to_UR50 and SP_seqs_dissimilar_to_UR50
SP_seqs_similar_to_UR50_split_dict = pd.read_pickle(SwissProt_split_dict_file)



# Extract the sequence IDs dissimilar to UniRef50
SP_seqs_dissimilar_to_UR50_clstr = SP_seqs_similar_to_UR50_split_dict["all_swissprot_clusters"]
print("SP_seqs_dissimilar_to_UR50_clstr      : " , SP_seqs_dissimilar_to_UR50_clstr[:3]  )
print("len(SP_seqs_dissimilar_to_UR50_clstr) : " , len(SP_seqs_dissimilar_to_UR50_clstr) )

# Randomly select 10K from the sequences dissimilar to UniRef50 using shuffle
import random
random.seed(42)
random.shuffle(SP_seqs_dissimilar_to_UR50_clstr)

num_clstr = len(SP_seqs_dissimilar_to_UR50_clstr)
clstr_for_valid_set = SP_seqs_dissimilar_to_UR50_clstr[                : (num_clstr // 6)     ] # values here adjust the number of datapoints in valid set.
clstr_for_test_set  = SP_seqs_dissimilar_to_UR50_clstr[ num_clstr // 6 : (num_clstr // 6) * 2 ] # values here adjust the number of datapoints in test  set.

SP_id_seqs_dissimilar_to_UR50_valid_set = [seq_id for clstr in clstr_for_valid_set for seq_id in clstr]
SP_id_seqs_dissimilar_to_UR50_test_set  = [seq_id for clstr in clstr_for_test_set  for seq_id in clstr]

# Extract the sequences dissimilar to UniRef50
SP_seqs_dissimilar_to_UR50_valid_set = [SwissProt_dict[seq_id.split(">")[1]] for seq_id in SP_id_seqs_dissimilar_to_UR50_valid_set]
SP_seqs_dissimilar_to_UR50_test_set  = [SwissProt_dict[seq_id.split(">")[1]] for seq_id in SP_id_seqs_dissimilar_to_UR50_test_set ]


# Now, split the dataframe into three sets: train, validation, and test, based on 
# the sequence in SP_seqs_dissimilar_to_UR50_vali_set and SP_seqs_dissimilar_to_UR50_test_set
# Save to three different csv files in the output folder with names
# U02_Processed_SwissProt_Dataset_train.csv , 
# U02_Processed_SwissProt_Dataset_valid.csv ,  
# U02_Processed_SwissProt_Dataset_test.csv  .

# Create a mapping from sequences to their respective categories (train, valid, test)
seq_to_category = {}
for seq in SP_seqs_dissimilar_to_UR50_valid_set:
    seq_to_category[seq] = 'valid'
for seq in SP_seqs_dissimilar_to_UR50_test_set:
    seq_to_category[seq] = 'test'

# Function to categorize each sequence in the dataframe
def categorize_sequence(row):
    seq = row['Sequence']
    return seq_to_category.get(seq, 'train')

# Apply the categorization
Processed_SwissProt_df['Category'] = Processed_SwissProt_df.apply(categorize_sequence, axis=1)

# Split the dataframe based on the category
train_df = Processed_SwissProt_df[Processed_SwissProt_df['Category'] == 'train']
valid_df = Processed_SwissProt_df[Processed_SwissProt_df['Category'] == 'valid']
test_df  = Processed_SwissProt_df[Processed_SwissProt_df['Category'] == 'test' ]

# Add a column to the dataframes that contains the mapped reactions for all three dataframes
train_df.loc[:, 'reaction_mapped'] = train_df['Reaction'].map(reaction_to_mapped_reaction_dict)
valid_df.loc[:, 'reaction_mapped'] = valid_df['Reaction'].map(reaction_to_mapped_reaction_dict)
test_df .loc[:, 'reaction_mapped'] = test_df ['Reaction'].map(reaction_to_mapped_reaction_dict)

# Remove the rows with None in the any columns
train_df = train_df.dropna()
valid_df = valid_df.dropna()
test_df  = test_df.dropna()

# Reset the index of the dataframes
train_df = train_df.reset_index(drop = True)
valid_df = valid_df.reset_index(drop = True)
test_df  = test_df.reset_index (drop = True)


# Save to CSV files
train_df.drop(columns=['Category']).rename(columns={"Sequence": "sequences", "Reaction": "reaction"}).to_csv(output_folder / "U02_Processed_SwissProt_Dataset_train.csv")
valid_df.drop(columns=['Category']).rename(columns={"Sequence": "sequences", "Reaction": "reaction"}).to_csv(output_folder / "U02_Processed_SwissProt_Dataset_valid.csv")
test_df .drop(columns=['Category']).rename(columns={"Sequence": "sequences", "Reaction": "reaction"}).to_csv(output_folder / "U02_Processed_SwissProt_Dataset_test.csv" )














































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




