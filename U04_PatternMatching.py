#!/usr/bin/env python
# coding: utf-8
#====================================================================================================#
# The following code ensures the code work properly in 
# MS VS, MS VS CODE and jupyter notebook on both Linux and Windows.

from rdkit import Chem
from rdkit.Chem import AllChem

#====================================================================================================#
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

import pickle

import numpy as np
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
# 







train_set_path = "./Processed_Dataset/U02_Processed_SwissProt_Dataset_train.csv"
valid_set_path = "./Processed_Dataset/U02_Processed_SwissProt_Dataset_valid.csv"
test_set_path  = "./Processed_Dataset/U02_Processed_SwissProt_Dataset_test.csv"

train_set = pd.read_csv(train_set_path)
valid_set = pd.read_csv(valid_set_path)
test_set  = pd.read_csv(test_set_path)

train_set_rxn_strs_list = train_set['Reaction'].tolist()
valid_set_rxn_strs_list = valid_set['Reaction'].tolist()
test_set_rxn_strs_list  = test_set ['Reaction'].tolist()



# Read reaction rules from ./RetroRules/knime-ready-rules_mnx-all-forward_ECOLI-iJO1366.csv
rxn_rules_fwd_path = "./RetroRules/knime-ready-rules_mnx-all-forward_ECOLI-iJO1366.csv"
rxn_rules_bwd_path = "./RetroRules/knime-ready-rules_mnx-all-reverse_ECOLI-iJO1366.csv"
rxn_rules_fwd_df   = pd.read_csv(rxn_rules_fwd_path)
rxn_rules_bwd_df   = pd.read_csv(rxn_rules_bwd_path)



# Make a dict that maps reaction rules to EC numbers, from the "Rule" column to the "EC" column
rule_fwd_to_EC_dict = dict(zip(rxn_rules_fwd_df['Rule'], rxn_rules_fwd_df['EC number']))
rule_bwd_to_EC_dict = dict(zip(rxn_rules_bwd_df['Rule'], rxn_rules_bwd_df['EC number']))

rule_fwd_list = rxn_rules_fwd_df['Rule'].tolist()
rule_bwd_list = rxn_rules_bwd_df['Rule'].tolist()




for one_rxn_str in train_set_rxn_strs_list:

    one_rxn_rctt = one_rxn_str.split(">>")[0]
    one_rxn_prod = one_rxn_str.split(">>")[1]

    one_rxn_rctt_list = one_rxn_rctt.split(".")
    one_rxn_prod_list = one_rxn_prod.split(".")

    print("Processing Reactants: ", one_rxn_rctt_list)

    for one_rule in rule_fwd_list:

        #one_rule_rctt_list = one_rule.split(">>")[0].split(".")
        
        r_fmt_rxn       = AllChem.ReactionFromSmarts(one_rule)
        r_fmt_rctt_list = [Chem.MolFromSmiles(rctt) for rctt in one_rxn_rctt_list]
        r_fmt_tctt_tupl = tuple(r_fmt_rctt_list)

        try:
            r_fmt_products_lists = r_fmt_rxn.RunReactants(r_fmt_tctt_tupl)
            if r_fmt_products_lists:
                print("Successfully simulate a reaction:", r_fmt_products_lists)
        except Exception:
            r_fmt_products_lists=[]
            # print("Failed to simulate a reaction:", r_fmt_products_lists)






















































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




