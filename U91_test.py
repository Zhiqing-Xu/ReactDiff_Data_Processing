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
#--------------------------------------------------#
import math
#--------------------------------------------------#

from typing   import Tuple, Dict
#--------------------------------------------------#
import torch
import torch.nn as nn
import torch.nn.functional as F


from ZX03_rxn_mpnn_args     import TrainArgs





#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Get X_TrainArgs


cmpd_process            = "d-MPNN"
MPNN_size               = "1280"


# Extra atom-level features.
cmpd_encodings_dict     = None                                                   # TODO Get dict from file.
X_cmpd_encodings_dim    = \
    len(list(cmpd_encodings_dict.values())[0]) if cmpd_encodings_dict != None else 0
    


# Compound side train arguments.
"""
'--features_path'
    - Whether to use MPNN only
    - comment out -> only MPNN
    - keep -> DONT EXCLUDE mol-level extra features
'--features_only'
    - Whether to use molecule-level features only
    - keep -> only mol-level extra features
    - comment out -> DONT EXCLUDE MPNN
'--undirected'
    - Whether to use directed message passing 
    - (only w/ bond_messages, ERROR w/ atom_messages)
'--no_atom_descriptor_scaling'
    - CANNOT set directly without providing extra features
'--no_bond_features_scaling'
    - CANNOT set directly without providing extra features
"""


arguments = ['--data_path'                                ,    ''                   ,
             '--dataset_type'                             ,    'regression'         ,
             '--hidden_size'                              ,    MPNN_size            ,
             '--depth'                                    ,    '3'                  ,

             '--atom_messages'                            , 
             '--reaction'                                 , 

             # DO NOT APPLY standardization on features
             '--no_features_scaling'                      , 

             # '--overwrite_default_atom_features'        , 
             # '--overwrite_default_bond_features'        , 

            ]


# Best Performance
if   cmpd_process == "u-MPNN-a": 
    pass

# Best for Testing
elif cmpd_process == "u-MPNN": 
    substruct_encoding_file = None

# Default `d-MPNN`
#   - No extra molecule/atom level features
#   - directed + bond message
elif cmpd_process == "d-MPNN": 
    arguments.remove('--atom_messages')
    substruct_encoding_file = None

# u-MPNN + extra molecule-level features
elif cmpd_process == "u-MPNN-a-m": 
    arguments.append('--features_path')

# Extra molecule-level features ONLY
elif cmpd_process == "Mol_Feat": 
    arguments.append('--features_path')
    arguments.append('')
    arguments.append('--features_only')
    substruct_encoding_file = None

# Exceptions
else:
    raise Exception(f"{cmpd_process} compound encoder type unknown!" )


# Use the "tap" package here, better than argparse
cmpd_TrainArgs = TrainArgs().parse_args(arguments)


# Additional Settings to be determined later (TODO)
# cmpd_TrainArgs.use_input_features = False

# End Configurations.




#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Read csv files 

import pandas as pd

train_set_path = "./data/U02_Processed_SwissProt_Dataset_train.csv"
valid_set_path = "./data/U02_Processed_SwissProt_Dataset_valid.csv"
test_set_path  = "./data/U02_Processed_SwissProt_Dataset_test.csv"

train_set = pd.read_csv(train_set_path)
valid_set = pd.read_csv(valid_set_path)
test_set  = pd.read_csv(test_set_path)

train_set_reaction_strings_list = train_set['reaction_smiles'].tolist()
valid_set_reaction_strings_list = valid_set['reaction_smiles'].tolist()
test_set_reaction_strings_list  = test_set ['reaction_smiles'].tolist()




















































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





















