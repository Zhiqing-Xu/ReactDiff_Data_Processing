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
                           'display.max_columns'   , 6    , 
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



import json
import numpy as np
from tqdm import tqdm



def write_train_sequences_to_fasta(data_dir, output_dir, output_file):
    """
    Write all training set sequences from the UniRefDataset to a new FASTA file.
    
    Parameters:
    - data_dir: The directory containing the dataset files.
    - output_file: The path to the output FASTA file.
    """
    # Load training indices
    with open(data_dir / 'splits.json', 'r') as f:
        train_indices = json.load(f)['train']
    
    # Load sequence offsets
    metadata = np.load(data_dir / 'lengths_and_offsets.npz')
    seq_offsets = metadata['seq_offsets']
    
    # Open the output FASTA file for writing
    with open(output_dir / output_file, 'w') as fasta_file:
        # Open the consensus.fasta file for reading
        with open(data_dir / 'consensus.fasta') as f:
            for idx in tqdm(train_indices, desc = "Writing sequences"):
                # Retrieve the sequence offset
                offset = seq_offsets[idx]
                # Seek to the position of the sequence in the file
                f.seek(offset)
                # Read the sequence
                consensus = f.readline().strip()
                # Write the sequence to the output file with the specified naming convention
                fasta_file.write(f">uniref50_seq_{idx}\n{consensus}\n")

# Example usage
data_dir    = Path('./_DATABASE_uniref50/')  # Adjust this to your data directory
output_dir  = Path('./Processed_Dataset/' )  # Adjust this to your desired output directory
output_file = 'U00_UniRef50_Sequences.fasta'  # Adjust this to your desired output file path
write_train_sequences_to_fasta(data_dir, output_dir, output_file)



































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




