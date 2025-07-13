#! /usr/bin/env python3
from copy import deepcopy

# Import util functions
from utils import update_checklist
from utils import compute_sum_extreme
from utils import benjamini_hochberg

def compute_likelihood_statistics(result_dict_full, result_dict_null, parameters):

    # Update the checklist
        
    tasks = [
    ("1. Build and fit full model", True),
    ("2. Create a joint model and null dataset", True),
    ("3. Evaluate and predict models", True),
    ("4. Build and fit null model", True),
    ("5. Compute likelihood ratio test statistics", False),
    ("6. Combine and create result files", False)
    ]

    update_checklist(tasks)

    # exclude poorly fitted ids from statistic calculation
    exclude_poor_fits = parameters['exclude_poor_fits']
    
    # Get dataframes with mll and lr values
    mll_df_full = result_dict_full['full_mll_values']
    lr_df_full = result_dict_full['lr_values_full_vs_joint']
    lr_df_null = result_dict_null['lr_values_full_vs_joint']

    # Filter data to only contain proteins with "good" fit (mll > 0.5)
    mll_df_full_filt = mll_df_full[mll_df_full['mll'] > 0.5]
    mll_df_full_filt['count'] = mll_df_full_filt.groupby('protein')['protein'].transform('size')
    mll_df_full_filt = mll_df_full_filt[mll_df_full_filt['count'] > 1]
    mll_df_full_filt = mll_df_full_filt.drop_duplicates('protein')
    
    if exclude_poor_fits == True:
        # Create a clean protein column for the null dataset
        lr_df_null['protein_clean'] = lr_df_null['protein'].str.replace(r'_sample_\d+', '', regex=True)
        # Filter lr dataframes to only contain proteins with "good" fits
        lr_df_full = lr_df_full[lr_df_full['protein'].isin(mll_df_full_filt['protein'])]
        lr_df_null = lr_df_null[lr_df_null['protein_clean'].isin(mll_df_full_filt['protein'])]
    else:
        None

    # Create a new dataframe to copy p-Values to
    p_val = lr_df_full.copy()
    
    # Define the null distribution approximation of the statistic (use LR values of null dataset (sampled from joint model))
    null_distribution = lr_df_null['lr'].tolist()
    
    # Compute the size of the null distribution approximation
    p_val['size_null'] = len(null_distribution)
    
    # Compare LR values of real dataset to LR of null dataset (sampled from joined model and calculate the number of values as extreme in the null distribution approximation
    p_val['n_lr_values_extreme'] = p_val['lr'].apply(compute_sum_extreme, approximated_likelihood_ratio=null_distribution)
    
    # Calculate p-Value
    p_val['pValue'] = (p_val['n_lr_values_extreme'] + 1) / (p_val['size_null'] + 1)
    
    # Adjust p-Value using Benjamini-Hochberg
    p_val['adj_pValue'] = benjamini_hochberg(p_val['pValue'])

    # copy results to result dict
    # Add predictions to result dict
    likelihood_statistics_result_dict = deepcopy(result_dict_full)
    likelihood_statistics_result_dict.update({
            "gp_likelihood_statistics_df" : p_val})
    
    # Update checklist
    tasks[4] = ("5. Compute likelihood ratio test statistics", True)
    update_checklist(tasks)
    return(likelihood_statistics_result_dict) 