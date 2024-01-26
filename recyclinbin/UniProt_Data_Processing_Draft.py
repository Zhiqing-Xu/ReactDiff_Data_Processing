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




from google.colab import drive
drive.mount('/content/drive')

import pandas as pd

columns = ["Sequence", "Catalytic activity"]



df = pd.read_csv('/content/drive/MyDrive/uniprotkb_AND_reviewed_true_2023_10_23.tsv', sep='\t',usecols=columns)

df.head(2)

"""# "sequence,[reation,[chebi],[Rhea],EC],[reation,[chebi],[Rhea],EC]"
"""

import pandas as pd
import re


def extract_data(activity):

    # Extracting the Reaction after CATALYTIC ACTIVITY: Reaction=
    reaction = re.search(r"CATALYTIC ACTIVITY:.*Reaction=(.*?);", activity)
    reaction = reaction.group(1) if reaction else " "

    # Extracting the numbers after ChEBI:CHEBI:
    chebi_numbers = re.findall(r"ChEBI:CHEBI:(\d+)", activity)
    chebi_numbers = chebi_numbers if chebi_numbers else [" "]

    # Extracting the string starting with Rhea:
    rhea = re.findall(r"Rhea:RHEA(?:-COMP)?:(\d+)", activity)
    rhea = rhea if rhea else [" "]

    # Extracting the number after EC=
    ec = re.search(r"EC=(\d+\.\d+\.\d+\.\d+);", activity)
    ec = ec.group(1) if ec else " "

    return [reaction, chebi_numbers, rhea, ec]

def process_row(row):
    sequence = row["Sequence"]
    activities_list = []
    if isinstance(row["Catalytic activity"], str):
        # Split by the pattern but without removing it
        activities = re.split(r"(?=CATALYTIC ACTIVITY: Reaction=)", row["Catalytic activity"])
        activities = [activity for activity in activities if activity]  # Removing any empty strings
        activities_data = [extract_data(activity) for activity in activities]
        activities_list.extend(activities_data)
    return [sequence, activities_list]

result = df.apply(process_row, axis=1).tolist()

# Printing the result
for item in result:
    print(item)

exploded_data = []

for item in result:
    sequence, reactions = item
    if reactions:
        for reaction in reactions:
            exploded_data.append([sequence, reaction])
    else:
        exploded_data.append([sequence, []])

# Print first 10 exploded data points to inspect
for i in range(min(10, len(exploded_data))):
    print(exploded_data[i])