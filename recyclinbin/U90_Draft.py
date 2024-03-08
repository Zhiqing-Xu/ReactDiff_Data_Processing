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
from random import shuffle



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

# 
output_folder = Path("Processed_Dataset/")
cluster_file  = output_folder / "U00_SwissProt_Sequences_CDhit_clstr.sorted"
fasta_file    = output_folder / "U00_SwissProt_Sequences.fasta"



def read_clusters(filename):
    clusters = {}
    seq_to_cluster = {}
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('>Cluster'):
                cluster_id = line.strip().split()[1]
                clusters[cluster_id] = []
            else:
                part = line.split(',')[1].strip()
                seq_id = part.split('...')[0].strip()
                is_representative = '*' in part
                clusters[cluster_id].append((seq_id, is_representative))
                seq_to_cluster[seq_id] = cluster_id

    # Debugging: Print some entries from seq_to_cluster
    print(list(seq_to_cluster.items())[:5])  # Print first 5 entries

    return clusters, seq_to_cluster



def split_clusters(clusters, split_ratio=0.95):
    all_representatives = [(cluster, seq_id) for cluster, members in clusters.items() for seq_id, is_rep in members if is_rep]
    shuffle(all_representatives)  # Randomize to ensure a fair split
    
    split_point = int(len(all_representatives) * split_ratio)
    train_clusters = set(cluster for cluster, _ in all_representatives[:split_point])
    test_clusters = set(cluster for cluster, _ in all_representatives[split_point:])

    # Debugging: Print cluster assignments
    print(f"Number of clusters in training set: {len(train_clusters)}")
    print(f"Number of clusters in test set: {len(test_clusters)}")

    return train_clusters, test_clusters


def read_fasta_sequences(filename, seq_to_cluster, train_clusters, test_clusters):
    seqs_train_cdhit = []
    seqs_test_cdhit = []
    with open(filename, 'r') as file:
        current_seq = []
        for line in file:
            if line.startswith('>'):
                # Process the previous sequence
                if current_seq:
                    # Correcting the format of the sequence ID to include '>'
                    seq_id = '>' + current_seq[0][1:].strip().split()[0]
                    cluster_id = seq_to_cluster.get(seq_id)
                    if cluster_id is None:
                        print(f"Warning: Sequence ID {seq_id} not found in cluster mapping.")
                    elif cluster_id in train_clusters:
                        seqs_train_cdhit.append(''.join(current_seq))
                    elif cluster_id in test_clusters:
                        seqs_test_cdhit.append(''.join(current_seq))
                    else:
                        print(f"Warning: Cluster ID {cluster_id} for sequence {seq_id} not found in train or test sets.")
                current_seq = [line]
            else:
                current_seq.append(line)

        # Handle the last sequence
        if current_seq:
            seq_id = '>' + current_seq[0][1:].strip().split()[0]
            cluster_id = seq_to_cluster.get(seq_id)
            if cluster_id in train_clusters:
                seqs_train_cdhit.append(''.join(current_seq))
            elif cluster_id in test_clusters:
                seqs_test_cdhit.append(''.join(current_seq))

    return seqs_train_cdhit, seqs_test_cdhit


# Read clusters and create a sequence to cluster mapping
clusters, seq_to_cluster = read_clusters(cluster_file)

# Debugging: Print number of sequences in clusters and in mapping
print(f"Total number of clusters: {len(clusters)}")
print(f"Total number of sequences in mapping: {len(seq_to_cluster)}")

# Split clusters into training and test sets
train_clusters, test_clusters = split_clusters(clusters, split_ratio=0.95)

# Read sequences from the fasta file and split into training and test sets
seqs_train_cdhit, seqs_test_cdhit = read_fasta_sequences(fasta_file, seq_to_cluster, train_clusters, test_clusters)

# Print results
print(f"Number of sequences in training set: {len(seqs_train_cdhit)}")
print(f"Number of sequences in test set: {len(seqs_test_cdhit)}")

# Now seqs_train_cdhit and seqs_test_cdhit contain the sequences








# Create a function to write a CSV file from a list of sequences
def create_csv_from_sequences(filename, sequences):
    # Extract only the sequence data, excluding the sequence names
    sequence_data = ["".join(seq.split("\n")[1:]) for seq in sequences]

    # Create a DataFrame with the specified columns and data
    df = pd.DataFrame({
        'sequences': sequence_data,
        'reaction': 'CCO>>CCO',
        'reaction_mapped': '[CH3:1][CH2:2][OH:3]>>[CH3:1][CH2:2][OH:3]'
    })

    # Reset index to start from 1 instead of 0
    df.index += 1
    df.reset_index(inplace=True)

    # Rename the index column to match the specified format
    df.rename(columns={'index': ''}, inplace=True)

    # Write the DataFrame to a CSV file
    df.to_csv(filename, index=False)

# File paths for the CSV files
train_csv_file = 'processed_SwissProt_dataset_sim_train.csv'
test_csv_file  = 'processed_SwissProt_dataset_sim_test.csv'

# Create CSV files for training and test sets
create_csv_from_sequences(train_csv_file, seqs_train_cdhit)
create_csv_from_sequences(test_csv_file, seqs_test_cdhit)












































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




