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


# 
Step_code = "U01_"

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#   `7MMF'      ,M'   .g8""8q.        `7MM"""Mq.    db   MMP""MM""YMM`7MMF'  `7MMF'     ,gM""bg         `7MM"""Yb. `7MMF'`7MM"""Mq.                    #
#     MM        MV  .dP'    `YM.        MM   `MM.  ;MM:  P'   MM   `7  MM      MM       8MI  ,8           MM    `Yb. MM    MM   `MM.                   #
#     MM       AW   dM'      `MM        MM   ,M9  ,V^MM.      MM       MM      MM        WMp,"            MM     `Mb MM    MM   ,M9                    #
#     MM      ,M'   MM        MM        MMmmdM9  ,M  `MM      MM       MMmmmmmmMM       ,gPMN.  jM"'      MM      MM MM    MMmmdM9                     #
#     MM      MV    MM.      ,MP        MM       AbmmmqMA     MM       MM      MM      ,M.  YMp.M'        MM     ,MP MM    MM  YM.                     #
#     MM     AW     `Mb.    ,dP'        MM      A'     VML    MM       MM      MM      8Mp   ,MMp         MM    ,dP' MM    MM   `Mb.                   #
#   .JMML.  ,M'       `"bmmd"'        .JMML.  .AMA.   .AMMA..JMML.   .JMML.  .JMML.    `YMbmm'``MMm.    .JMMmmmdP' .JMML..JMML. .JMM.                  #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# 
UniProt_Data_file    = "_DATABASE_uniprotkb_AND_reviewed_true_2023_10_23.tsv"
UniProt_Data_columns = ["Sequence", "Catalytic activity"]


output_folder        = Path("Processed_Dataset/")
output_folder.mkdir(parents = True, exist_ok = True)


output_file = output_folder / (Step_code + "SwissProt_split_dict.p")






#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#  `7MM"""Mq. `7MM"""YMM        db     `7MM"""Yb.         .g8"""bgd `7MMF'    `7MMF'   `7MF'.M"""bgd MMP""MM""YMM `7MM"""YMM  `7MM"""Mq.   .M"""bgd    #
#    MM   `MM.  MM    `7       ;MM:      MM    `Yb.     .dP'     `M   MM        MM       M ,MI    "Y P'   MM   `7   MM    `7    MM   `MM. ,MI    "Y    #
#    MM   ,M9   MM   d        ,V^MM.     MM     `Mb     dM'       `   MM        MM       M `MMb.          MM        MM   d      MM   ,M9  `MMb.        #
#    MMmmdM9    MMmmMM       ,M  `MM     MM      MM     MM            MM        MM       M   `YMMNq.      MM        MMmmMM      MMmmdM9     `YMMNq.    #
#    MM  YM.    MM   Y  ,    AbmmmqMA    MM     ,MP     MM.           MM      , MM       M .     `MM      MM        MM   Y  ,   MM  YM.   .     `MM    #
#    MM   `Mb.  MM     ,M   A'     VML   MM    ,dP'     `Mb.     ,'   MM     ,M YM.     ,M Mb     dM      MM        MM     ,M   MM   `Mb. Mb     dM    #
#  .JMML. .JMM.JMMmmmmMMM .AMA.   .AMMA.JMMmmmdP'         `"bmmmd'  .JMMmmmmMMM  `bmmmmd"' P"Ybmmd"     .JMML.    .JMMmmmmMMM .JMML. .JMM.P"Ybmmd"     #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Load the dataset
clusters            = {}    # Dictionary to hold clusters
current_cluster_id  = None
swissprot_sequences = set()
uniref50_sequences  = set()

SwissProt_seq_tag   = ["SwissProt_seq_" , "SP_"   , ][1]
UniRef50_seq_tag    = ["uniref50_seq_"  , "UR50_" , ][1]



# Parse the .clstr file to organize sequences by cluster
with open('./Processed_Dataset/clustered.fasta.clstr', 'r') as file:

    # Get the number of lines in the file.
    line_num = "80M"
    print("The number of lines in the file: ", line_num)


    for line_idx, line in enumerate(file):

        if line_idx % 1000000 == 0:
            print(f"Processed {int(line_idx / 1000000)}M lines out of {line_num}.")

        if line.startswith('>Cluster'):
            current_cluster_id = line.strip().split(' ')[1]
            clusters[current_cluster_id] = []

        else:
            #seq_id = line.split('>')[1].split('...')[0]
            seq_id = line.split('\t')[1].split('...')[0]
            clusters[current_cluster_id].append(seq_id)

            if   SwissProt_seq_tag in seq_id:
                swissprot_sequences.add(seq_id)

            elif UniRef50_seq_tag in seq_id:
                uniref50_sequences.add(seq_id)



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#     .M"""bgd `7MM"""Mq.`7MMF'      `7MMF'MMP""MM""YMM       .M"""bgd `7MM"""YMM    .g8""8q.    .M"""bgd                                              #
#    ,MI    "Y   MM   `MM. MM          MM  P'   MM   `7      ,MI    "Y   MM    `7  .dP'    `YM. ,MI    "Y                                              #
#    `MMb.       MM   ,M9  MM          MM       MM           `MMb.       MM   d    dM'      `MM `MMb.                                                  #
#      `YMMNq.   MMmmdM9   MM          MM       MM             `YMMNq.   MMmmMM    MM        MM   `YMMNq.                                              #
#    .     `MM   MM        MM      ,   MM       MM           .     `MM   MM   Y  , MM.    M ,MP .     `MM                                              #
#    Mb     dM   MM        MM     ,M   MM       MM           Mb     dM   MM     ,M `Mb.   AMMP' Mb     dM                                              #
#    P"Ybmmd"  .JMML.    .JMMmmmmMMM .JMML.   .JMML.         P"Ybmmd"  .JMMmmmmMMM   `"bmmdPVM  P"Ybmmd"                                               #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Identify SwissProt sequences clustered with UniRef50 sequences
import copy
swissprot_dissimilar   = set()
swissprot_similar      = set()
all_swissprot_clusters = []

len_clusters = len(clusters) #37M
print("Number of clusters: ", len_clusters)

for idx, (clstr_k, clstr_v)  in enumerate(clusters.items()):
    
    if idx % 1000000 == 0:
        print("\ncluster_k: ", clstr_k, "\ncluster_v: ", clstr_v[:5])
        print("Processing cluster ", int(idx / 1000000), "M out of ", "37M", ".")


    if any(SwissProt_seq_tag in seq for seq in clstr_v) and any(UniRef50_seq_tag in seq for seq in clstr_v):
    # If the cluster contains both SwissProt and UniRef50 sequences
        pass
    
    elif all(SwissProt_seq_tag in seq for seq in clstr_v):
    # If the cluster contains only SwissProt sequences
        all_swissprot_clusters.append(clstr_v)

        for seq in clstr_v:
            swissprot_dissimilar.add(seq)


# SwissProt sequences not clustered with UniRef50 sequences
swissprot_similar = swissprot_sequences - swissprot_dissimilar



# Reporting
print(f"Total SwissProt sequences                        : {len(swissprot_sequences  )}"  ) # 
print(f"Total UniRef50 sequences                         : {len(uniref50_sequences   )}"  ) # 
print(f"SwissProt sequences similar to UniRef50          : {len(swissprot_similar    )}"  ) # 
print(f"SwissProt sequences with no UniRef50 similarity  : {len(swissprot_dissimilar )}"  ) # 58K


# Save swissprot_similar and swissprot_dissimilar to a pickle file which contains a dictionary.
split_dict = { "SP_seqs_similar_to_UR50"    : swissprot_similar      , 
               "SP_seqs_dissimilar_to_UR50" : swissprot_dissimilar   , 
               "all_swissprot_clusters"     : all_swissprot_clusters , }

pd.to_pickle(split_dict, output_file)





























































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




