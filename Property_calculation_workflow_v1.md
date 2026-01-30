# Workflow used to calculate physicochemical properties of 217 single-site monomeric fungicides with known modes of action. 
#### 384 Properties were calculated using RDKit, and another 11 quantum properties were calculated using ORCA.

## RDKit Calculations
RDKit Molecular descriptors were calculated from SMILES codes using RDKit (version 2025.3.3). 
SMILES strings were imported from CSV files and converted to RDKit molecule objects using Chem.MolFromSmiles. 
Entries with missing or invalid SMILES were excluded. 
A comprehensive set of RDKit molecular descriptors was computed using rdkit.Chem.Descriptors and MolecularDescriptorCalculator. 
All callable descriptor functions were included except setupAUTOCorrDescriptors and PropertyFunctor, which were excluded due to incompatibility with batch calculation. 
Descriptor values were calculated for each molecule and exported together with compound names and SMILES identifiers as a CSV table for downstream analysis. 
Python (3.12.12) script for RDKit (2025.3.3) workflow is shown in "rdkit_calculations.py".

## ORCA Calculations
Two computational approaches were employed for ORCA (6.1) calculations, depending on molecular complexity. 
### Method 1 (Direct optimization, n=213): 
Geometry optimizations and single-point energy calculations were performed using the RI-B2PLYP double-hybrid density functional with Grimme's D3BJ dispersion correction and the def2-TZVP basis set. 
The RIJCOSX approximation was used for computational efficiency with the def2/J auxiliary basis for Coulomb fitting and def2-TZVP/C for correlation fitting. 
Tight SCF convergence criteria (TightSCF) and tight optimization thresholds (TightOpt) were applied. 
The input keyword line was: 
“! RIJCOSX RI-B2PLYP D3BJ def2-TZVP def2/J def2-TZVP/C TightSCF TightOpt Opt”
### Method 2 (Two-step approach, n=4): 
For four structurally complex fungicides that encountered convergence issues with Method 1 (fluoxapiprolin, validamycin, fenpicoxamid, iminoctadine), a two-step protocol was implemented. 
First, geometries were optimized using the computationally efficient B97-3c composite method with tight convergence criteria: 
“ ! B97-3c Opt TightSCF TightOpt” 
Subsequently, single-point energy calculations were performed on the B97-3c-optimized geometries using the same RI-B2PLYP level of theory as Method 1: 
“ ! RIJCOSX RI-B2PLYP D3BJ def2-TZVP def2/J def2-TZVP/C TightSCF”

### Dipole moment calculations: 
Dipole moments were calculated at the SCF level for all compounds to ensure consistency. 
For two compounds (dodine and mandipropamid), dipole moment calculations were inadvertently omitted from the initial Method 1 because they failed to achieve optimization convergence. 
For these compounds, additional single-point calculations were performed using the final optimized geometries from Method 1 with explicit dipole moment requests to maintain consistency with the remaining dataset.
Extracted using RStudio script 
