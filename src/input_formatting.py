#! /usr/bin/env python3
# -*- coding: utf-8 -*-
## Dependencies
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import os
import glob
import re

# Function to process TMT-MS2 output (MSfragger protein.tsv files) for further TPP analysis using Gaussian processing
#############################################################################################################################################################
"""
Author: Johannes F. Hevler [jfhevler@stanford.edu]

Function to format and normalize protein intensities followed by min-max intensity scaling. 
"""

def msfragger_to_gp(file_dir, output_dir, qupm_filter = 1,
                    na_filter = 'Yes', temp_gradient = [37, 41, 45, 50, 54, 57, 58, 62, 68, 72],
                    control_condition='vehicle', compound_conc={'vehicle': 0, 'ngi1': 10, 'kifunensine': 10, '3fax':50, '2ff': 100}):
    
    # Check whether required naming scheme for files in directoy is correct:
    # Define the directory and pattern
    directory = file_dir
    pattern = r'^protein_[a-zA-Z0-9]+_[0-9]+\.tsv$'
    for file_name in os.listdir(directory):
        if file_name.endswith('.tsv'):
            if re.match(pattern, file_name):
                print(f"{file_name}: File naming matches the naming scheme")
            else:
                raise ValueError(f"{file_name}: File naming does not match the naming scheme! Please adhere to following scheme: protein_*condition*_*replicate*.tsv (e.g. protein_kifunensine_1.tsv)")
            
    # Load files from directory
    files = glob.glob(file_dir+'/*.tsv')
    tpptr_dfs = pd.concat([pd.read_csv(fp,sep = "\t").assign(Filename=os.path.basename(fp)) for fp in files])
    
    # remove contaminations
    tpptr_dfs = tpptr_dfs[~tpptr_dfs['Protein'].str.contains('contam_')]

    # For proteins without Gene name set Protein ID as Gene name
    tpptr_dfs['Gene'] = np.where(tpptr_dfs['Gene'].isna(), tpptr_dfs['Protein ID'], tpptr_dfs['Gene'])

    # Add columns with information about tested condition and replicate number
    tpptr_dfs['Condition'] = tpptr_dfs['Filename'].astype(str).str.extract('(.*\_(.*?)\_.*)', expand=True)[1]
    tpptr_dfs['Sample'] = tpptr_dfs['Filename'].astype(str).str.extract('(.*protein\_(.*?)\.tsv.*)', expand=True)[1]

    # reduce size of dataframe by extracting columns of interest
    tpptr_dfs = tpptr_dfs.iloc[:, [1,3,11,14,17]+list(range(23,36))]
    
    # rename columns to fit column naming of input table for Bioconductor TPP R packages 
    new_col_names = ["accession", "gene_name", "qupm", "qssm", "MS1_intensity", "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "filename", "condition", "sample"]
    tpptr_dfs.columns = new_col_names

    # remove proteins without quantified MS1 and abundance in T1 + T2 + T3 + T4 + T5
    tpptr_dfs = (tpptr_dfs[~(tpptr_dfs['T1'] == 0) & ~(tpptr_dfs['T2'] == 0) & ~(tpptr_dfs['T3'] == 0) & ~(tpptr_dfs['T4'] == 0) & ~(tpptr_dfs['T5'] == 0)])
    
    # remove proteins that are below the set qssm filter
    if qupm_filter != None:
        tpptr_dfs = tpptr_dfs[tpptr_dfs['qupm'] > qupm_filter]
    
    else:
        pass     

    # Identify rows with the same Genename but different AccessionCode for each condition
    mask = tpptr_dfs.groupby(['gene_name', 'sample'])['accession'].transform('nunique') > 1

    # Add a suffix to the GeneName to make it unique within each group
    tpptr_dfs.loc[mask, 'gene_name'] += '_' + tpptr_dfs[mask].groupby(['gene_name', 'sample']).cumcount().astype(str)

    # Reset the index to avoid duplicate labels
    tpptr_dfs.reset_index(drop=True, inplace=True)

    # return table with Uniprot Accession, Gene name and Intensity
    tpptr_dfs = tpptr_dfs[["accession", "gene_name", "qupm", "qssm", "condition", "sample", "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10"]]

    # log2 transform intensity values to perform a median normalization
    for column in tpptr_dfs[["T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10"]]:
            tpptr_dfs[column] = np.log2(tpptr_dfs[column]+1e-8)

    # replace -inf values (log2(0)) with NaN
    tpptr_dfs.replace([np.inf, -np.inf, 0], np.nan, inplace=True)

    if na_filter == 'Yes':
        tpptr_dfs = tpptr_dfs.dropna()
    
    else:
        pass 

    # melt dataframe to perform median normalization
    col_names = tpptr_dfs.iloc[:, 6:16].columns
    tpptr_dfs = tpptr_dfs.melt(id_vars=["accession", "gene_name", "qupm", "qssm", "condition", "sample"], value_vars=col_names)

    # remove NaN values (only removes single values not the whole protein)
    tpptr_dfs = tpptr_dfs.dropna()
    
    # median normalization
    tpptr_dfs['median_sample'] = tpptr_dfs.groupby(['sample','variable'])['value'].transform('median')
    tpptr_dfs['median_all'] = tpptr_dfs.groupby(['variable'])['value'].transform('median')
    tpptr_dfs['norm_value'] = tpptr_dfs['value'] - (tpptr_dfs['median_sample'] - tpptr_dfs['median_all'])
    
      
    # Calculate scaled intensities
    
    ## Scale the data between 0 and 1 within each group
    def minmax_scale(group):
        min_val = group['norm_value'].min()
        max_val = group['norm_value'].max()
        group['scaled_intensity'] = (group['norm_value'] - min_val) / (max_val - min_val)
        return group

    # Group by 'Protein', 'condition', and 'replicate', then apply scaling
    tpptr_dfs = tpptr_dfs.groupby(['accession', 'sample']).apply(minmax_scale).reset_index(drop=True)

    # Map Temperatures
    temp_range = temp_gradient
    mapping_scheme_normalized = {
    'T1': temp_range[0],
    'T2': temp_range[1],
    'T3': temp_range[2],
    'T4': temp_range[3],
    'T5': temp_range[4],
    'T6': temp_range[5],
    'T7': temp_range[6],
    'T8': temp_range[7],
    'T9': temp_range[8],
    'T10': temp_range[9]
    }
    
    tpptr_dfs['temperature'] = tpptr_dfs['variable'].apply(lambda x: mapping_scheme_normalized.get(x, 'Other'))

    # Restore the original shape of the DataFrame   
    tpptr_dfs_restored = tpptr_dfs.pivot_table(index=["accession", "gene_name", "qupm", "qssm", "condition", "sample"], columns="temperature", values="scaled_intensity").reset_index()

    # Restore the original shape of the DataFrame but with anti log2 normalized values
    tpptr_dfs_undo_log2 = tpptr_dfs.copy()
    tpptr_dfs_undo_log2['norm_value_exp'] = tpptr_dfs_undo_log2["norm_value"].rpow(2)
    tpptr_dfs_restored_log2 = tpptr_dfs_undo_log2.pivot_table(index=["accession", "gene_name", "qupm", "qssm", "condition", "sample"], columns="temperature", values="norm_value_exp").reset_index()
    tpptr_dfs_restored_log2.columns = [str(col) if isinstance(col, (int, float)) else col for col in tpptr_dfs_restored_log2.columns]
        
    # Calculate relative abundances
    lowest_temp_index = tpptr_dfs_restored_log2.columns.get_loc('37')

    # Convert columns to numeric
    for col in tpptr_dfs_restored_log2.columns[lowest_temp_index:lowest_temp_index+10]:
        tpptr_dfs_restored_log2[col] = pd.to_numeric(tpptr_dfs_restored_log2[col], errors='coerce')

    # Iterate over each TMT channel column and compute relative abundances
    for col in tpptr_dfs_restored_log2.columns[lowest_temp_index:lowest_temp_index+10]:  # Skip the lowest temperature column
        tpptr_dfs_restored_log2['fc_'+col] = tpptr_dfs_restored_log2[col] / tpptr_dfs_restored_log2['37']        
    tpptr_dfs_restored_log2 = tpptr_dfs_restored_log2.iloc[:, [0,1,2,3,4,5]+list(range(16,26))]
    
    # add compound concentration information
    def combine_columns(row):
        return str(row['accession']) + '-(' + row['gene_name'] + ')'
    
    # add a unique Protein_ID column by combining Uniprot accession and Genename
    tpptr_dfs_restored['uniqueID'] = tpptr_dfs_restored.apply(combine_columns, axis=1)
    tpptr_dfs_restored_log2['uniqueID'] = tpptr_dfs_restored_log2.apply(combine_columns, axis=1)
    
    tpptr_dfs_restored['compoundConcentration'] = tpptr_dfs_restored['condition'].apply(lambda x: compound_conc.get(x, 'Other'))
    tpptr_dfs_restored_log2['compoundConcentration'] = tpptr_dfs_restored_log2['condition'].apply(lambda x: compound_conc.get(x, 'Other'))

    # format dataframe to fit the GP input and save dataframes for each sample
    drug_conditions = tpptr_dfs_restored[tpptr_dfs_restored['condition'] != control_condition]['condition'].unique()  

    for condition in drug_conditions:
                    
        # drug dataframe
        drug_df = tpptr_dfs_restored[(tpptr_dfs_restored['condition'] == condition) & (~tpptr_dfs_restored['sample'].str.contains(control_condition))]
        drug_df['replicate_ids'] = drug_df.groupby('uniqueID').transform('size')
        drug_df2 = tpptr_dfs_restored_log2[(tpptr_dfs_restored_log2['condition'] == condition) & (~tpptr_dfs_restored_log2['sample'].str.contains(control_condition))]
        drug_df2['replicate_ids'] = drug_df2.groupby('uniqueID').transform('size')  
            
        # control dataframe
        control_df = tpptr_dfs_restored[(tpptr_dfs_restored['condition'] == control_condition)]
        control_df['replicate_ids'] = control_df.groupby('uniqueID').transform('size')

        control_df2 = tpptr_dfs_restored_log2[(tpptr_dfs_restored_log2['condition'] == control_condition)]
        control_df2['replicate_ids'] = control_df2.groupby('uniqueID').transform('size')  

        # combine dataframes
        conditions_df = pd.concat([control_df, drug_df])
        conditions_df2 = pd.concat([control_df2, drug_df2])

        # add condition information
        conditions_df['condition'] = condition
        conditions_df['replicate'] = conditions_df['sample'].apply(lambda x: x.split('_')[-1])

        conditions_df2['condition'] = condition
        conditions_df2['replicate'] = conditions_df2['sample'].apply(lambda x: x.split('_')[-1])

        # extract columns of interest and re-order columns
        conditions_df = conditions_df.iloc[:, [16,2,4,19,17,18]+list(range(6,16))]
        conditions_df2 = conditions_df2.iloc[:, [16,2,4,19,17,18]+list(range(6,16))]
        new_col_names = ["uniqueID", "uniquePeptideMatches", "condition", "replicate","compoundConcentration","replicate_ids", "37", "41", "45", "50", "54", "57", "58", "62", "68", "72"]
        conditions_df.columns = new_col_names
        conditions_df2.columns = new_col_names

        # melt dataframe
        col_names = conditions_df.iloc[:, 6:16].columns
        conditions_df = conditions_df.melt(id_vars=["uniqueID", "uniquePeptideMatches", "condition", "replicate", "compoundConcentration", "replicate_ids"], value_vars=col_names)
        conditions_df2 = conditions_df2.melt(id_vars=["uniqueID", "uniquePeptideMatches", "condition", "replicate", "compoundConcentration", "replicate_ids"], value_vars=col_names)
        
        # change column names
        new_col_names = ["uniqueID", "uniquePeptideMatches", "condition", "replicate","compoundConcentration","replicate_ids", "x", "y"]
        conditions_df.columns = new_col_names
        conditions_df2.columns = new_col_names

        # change order of columns
        conditions_df = conditions_df[["condition", "uniqueID", "y", "x", "compoundConcentration", "replicate", "uniquePeptideMatches", "replicate_ids"]]
        conditions_df2 = conditions_df2[["condition", "uniqueID", "y", "x", "compoundConcentration", "replicate", "uniquePeptideMatches", "replicate_ids"]]
        
        # filter for proteins with a certain number of datapoints
        conditions_df['datapoints'] = conditions_df.groupby('uniqueID').transform('size')
        conditions_df2['datapoints'] = conditions_df2.groupby('uniqueID').transform('size')
        
        # Add max and min number of qupm
        conditions_df['max_qupm'] = conditions_df.groupby(['uniqueID'])['uniquePeptideMatches'].transform('max')
        conditions_df['min_qupm'] = conditions_df.groupby(['uniqueID'])['uniquePeptideMatches'].transform('min') 
        conditions_df2['max_qupm'] = conditions_df2.groupby(['uniqueID'])['uniquePeptideMatches'].transform('max')
        conditions_df2['min_qupm'] = conditions_df2.groupby(['uniqueID'])['uniquePeptideMatches'].transform('min')   
            
        # For compoundConcentration == 0, indicate that its the control
        if control_condition == 'vehicle':

            conditions_df['condition'] = np.where(conditions_df['compoundConcentration'] == 0, 'vehicle', conditions_df['condition'])
            conditions_df2['condition'] = np.where(conditions_df2['compoundConcentration'] == 0, 'vehicle', conditions_df2['condition'])
        else:
            control_comp_conc = compound_conc[control_condition]
            conditions_df['condition'] = np.where(conditions_df['compoundConcentration'] == control_comp_conc, control_condition, conditions_df['condition'])
            conditions_df2['condition'] = np.where(conditions_df2['compoundConcentration'] == control_comp_conc, control_condition, conditions_df2['condition'])   

        # save dataframe
        conditions_df.to_csv(f'{output_dir}/gpmelt_{condition}_scaled_intensities.csv', index=False)
        conditions_df2.to_csv(f'{output_dir}/gpmelt_{condition}_fc_intensities.csv', index=False)

        # inform the user that formatted datframes have been saved
        print(f"gpmelt_{condition}_scaled_intensities.csv and gpmelt_{condition}_fc_intensities.csv has been saved to {output_dir}!")    

    return(conditions_df)

"""Process and normalize intensities for GP analysis of TPP data.

MSfragger output tables are formatted and contaminations as well as proteins with low numbers of quantified peptides (the default minimum requirement is two quantified peptide matches) are removed. Next, protein abundances are min-max scaled for each replicate . Output is used to generate input files for GP regression of melting curves, as well as [slimTPCA](https://github.com/wangjun258/Slim_TPCA) and [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) analysis to assess changed thermal stability, protein-protein interactions, complex dynamics as well as protein abundance. 

Parameters
----------
file_dir : :class:'string'
Path to directory with MSfragger TMT10-MS2 protein.tsv files.

output_dir : :class:'string'
Path to directory where proccessed files are saved.

qupm_filter : :class:'numeric'
Threshold for removing proteins based on the number of quantified peptide matches. Default=1, removes proteins with only one quantified peptide.

na_filter : :class:'numeric'
Threshold for removing proteins with missing values. If set to "Yes", proteins with missing values are removed from analysis. 

temp_gradient : :class:'string'
A string indicating the applied temperature gradient for the TPP experiment. Default=37, 41, 45, 50, 54, 57, 58, 62, 68, 72 degrees celcius.

control_condition : :class:'string'
Name of the control condition, e.g. "vehicle".

compound_conc : :class:'list'
A list indicating the tested drug conditions and respective concentrations.

Returns
-------
tpptr_dfs_restored : :class:'pandas.DataFrame'
Dataframe that can be used as input to generate input files for gaussian processing modeling of melting curves, [RTPCA](https://www.bioconductor.org/packages/release/bioc/html/Rtpca.html)), and [slimTPCA](https://github.com/wangjun258/Slim_TPCA) analysis.

Usage
-----
In order to work, file naming needs to adhere to this format: protein_*condition*_*replicate*.tsv (e.g. protein_kifunensine_1.tsv)
msfragger_to_gp(file_dir = r'D:\Projects\TPP_altered_glycosylation\Data_analysis\EXP2\Data_analysis_12fr\Data_formatted_for_analysis\MsFragger_output', output_dir=r'D:\Projects\TPP_altered_glycosylation\Data_analysis\EXP2\Data_analysis_12fr\gp_modeling\data\min_max_scaling', qupm_filter = 1)
"""
