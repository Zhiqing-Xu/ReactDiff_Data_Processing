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



import re


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
    suffix_pattern = re.compile(r"\[\S+\]\((n[\+\-]?\d*)\)$")

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


# Example usage
eqn_list = ["[phosphate](n) + UTP = [phosphate](n+1) + UDP"                                                           ,
            "(R)-carnitine + H(+) + NADPH + O2 = (3R)-3-hydroxy-4-oxobutanoate + H2O + NADP(+) + trimethylamine"      ,
            "[(3R)-hydroxybutanoate](n) + H2O = (3R)-hydroxybutanoate dimer + [(3R)-hydroxybutanoate](n-2) + H(+)"    ,
            "[(3R)-hydroxybutanoate](n) + H2O = (3R)-hydroxybutanoate trimer + [(3R)-hydroxybutanoate](n-3) + H(+)"   ,
            "[(3R)-hydroxybutanoate](n) + H2O = (3R)-hydroxybutanoate tetramer + [(3R)-hydroxybutanoate](n-4) + H(+)" ,
            "[(3R)-hydroxybutanoate](n) + H2O = (3R)-hydroxybutanoate pentamer + [(3R)-hydroxybutanoate](n-5) + H(+)" ,

            ]


for one_eqn in eqn_list:
    preprocess_rxn_eqn(one_eqn)