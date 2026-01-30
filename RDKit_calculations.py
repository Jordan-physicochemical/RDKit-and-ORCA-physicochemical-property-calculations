# This script calculates many physicochemical properties for compounds input as a SMILES list.
# Python version: 3.12.12
# RDKit version: 2025.3.3
# Date: 2025-10-14
# Author: Jordan Campbell
# Adapted from google collab and github resources:
# https://colab.research.google.com/github/schwallergroup/ai4chem_course/blob/main/notebooks/01%20-%20Basics/01d_rdkit_basics.ipynb#scrollTo=4v2cbSCiRl4q
# https://github.com/schwallergroup/ai4chem_course/blob/main/notebooks/01%20-%20Basics/01d_rdkit_basics.ipynb

#!/usr/bin/env python3

# Install rdkit
!pip install rdkit

## Calculate (almost) all rdkit properties for 190 fungicides
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors.MoleculeDescriptors import MolecularDescriptorCalculator
import pandas as pd

# Upload Smiles_fungicide_list.csv
from google.colab import files
uploaded = files.upload()

df_new = pd.read_csv('Smiles_fungicide_list.csv')
df_new.head()

# Rename columns for clarity
df_new.columns = ['SMILES', 'Name']

# Remove rows with missing or non-string SMILES
df_new = df_new[df_new['SMILES'].apply(lambda x: isinstance(x, str) and x.strip() != '')].copy()

# Convert SMILES to RDKit Mol objects
df_new['Molecule'] = df_new['SMILES'].apply(Chem.MolFromSmiles)

# Filter out invalid molecules
df_new = df_new[df_new['Molecule'].notnull()].reset_index(drop=True)

# Get all descriptor names, excluding those causing errors
descriptor_names = [name for name in dir(Descriptors)
                    if not name.startswith('_') and callable(getattr(Descriptors, name))
                    and name not in ['setupAUTOCorrDescriptors', 'PropertyFunctor']] # Exclude problematic descriptors

# Create descriptor calculator
calculator = MolecularDescriptorCalculator(descriptor_names)

# Calculate descriptors
properties = df_new['Molecule'].apply(calculator.CalcDescriptors)
df_properties = pd.DataFrame(properties.tolist(), columns=descriptor_names)

# Combine with fungicide names and SMILES
df_final_new = pd.concat([df_new[['Name', 'SMILES']], df_properties], axis=1)

# Display the final descriptor table
display(df_final_new)

# Create and download csv file
df_final_new.to_csv('rdkit_fungicide_properties.csv', index=False)
files.download('rdkit_fungicide_properties.csv')

df_final_new.to_csv('20251014_rdkit_fungicide_properties_round10.csv', index=False)
files.download('20251014_rdkit_fungicide_properties_round10.csv')
