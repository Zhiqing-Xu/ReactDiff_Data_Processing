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
import json
import numpy as np
import pandas as pd
from tqdm import tqdm


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
# `7MMF'     A     `7MF'`7MM"""Mq. `7MMF'MMP""MM""YMM `7MM"""YMM    `7MMF'   `7MF'`7MM"""Mq.     `7MM"""YMM  db       .M"""bgd MMP""MM""YMM  db        #
#   `MA     ,MA     ,V    MM   `MM.  MM  P'   MM   `7   MM    `7      MM       M    MM   `MM.      MM    `7 ;MM:     ,MI    "Y P'   MM   `7 ;MM:       #
#    VM:   ,VVM:   ,V     MM   ,M9   MM       MM        MM   d        MM       M    MM   ,M9       MM   d  ,V^MM.    `MMb.          MM     ,V^MM.      #
#     MM.  M' MM.  M'     MMmmdM9    MM       MM        MMmmMM        MM       M    MMmmdM9        MM""MM ,M  `MM      `YMMNq.      MM    ,M  `MM      #
#     `MM A'  `MM A'      MM  YM.    MM       MM        MM   Y  ,     MM       M    MM  YM.        MM   Y AbmmmqMA   .     `MM      MM    AbmmmqMA     #
#      :MM;    :MM;       MM   `Mb.  MM       MM        MM     ,M     YM.     ,M    MM   `Mb.      MM    A'     VML  Mb     dM      MM   A'     VML    #
#       VF      VF      .JMML. .JMM.JMML.   .JMML.    .JMMmmmmMMM      `bmmmmd"'  .JMML. .JMM.   .JMML..AMA.   .AMMA.P"Ybmmd"     .JMML.AMA.   .AMMA.  #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# 
Step_code = "U00A_"


# SwissProt sequences that are not similar to EvoDiff train sequences can be used as the validation set.
# The following function loads all EvoDiff train sequences from the UniRef50 dataset and write them to a new FASTA file.
# This step is preparation for splitting the SwissProt dataset into training and validation sets.

def write_train_sequences_to_fasta(data_dir, output_dir, output_file):

    """
    Write all **training** set sequences from the UniRefDataset to a new FASTA file.
    
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
                fasta_file.write(f">UR50_{idx}\n{consensus}\n")

    return


# Run the function to write the training sequences to a new FASTA file.
data_dir    = Path('./_DATABASE_uniref50/')  # Adjust this to your data directory
output_dir  = Path('./Processed_Dataset/' )  # Adjust this to your desired output directory
output_file = Step_code + 'UR50_Sequences.fasta'  # Adjust this to your desired output file path
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




