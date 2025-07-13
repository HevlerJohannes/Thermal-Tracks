#! /usr/bin/env python3

import numpy as np
import torch
from matplotlib import pyplot as plt
from IPython.display import Markdown, display, clear_output
import pandas as pd
from copy import deepcopy
import seaborn as sns

#############################
########### UTILS ###########
#############################

def plot_posterior_distributions(ax, prot, plot_df, conditions):
    colors = ['#99d8c9', '#2ca25f']
    color_areas = ['#ff7f00', '#6a3d9a']

    # Plot the joint fitted curve
    ax.plot('x', 'y', label='joint', 
            data=plot_df[(plot_df['type'] == 'fitted') & (plot_df['condition'] == 'joint')],
            color='black', linestyle='--', alpha=0.5)

    # Plot the fitted curves and confidence intervals for each condition
    for cond, color, color_area in zip(conditions, colors, color_areas):
        ax.plot('x', 'y', label=f'fitted - {cond}', 
                data=plot_df[(plot_df['type'] == 'fitted') & (plot_df['condition'] == cond)], 
                color=color, linestyle='-')
        ax.fill_between(x=plot_df[(plot_df['type'] == 'fitted') & (plot_df['condition'] == cond)]['x'],
                        y1=plot_df[(plot_df['type'] == 'fitted') & (plot_df['condition'] == cond)]['conf_upper'],
                        y2=plot_df[(plot_df['type'] == 'fitted') & (plot_df['condition'] == cond)]['conf_lower'],
                        color=color_area, alpha=0.1, label=f'confidence interval - {cond}')
        # Plot measured data
        ax.scatter('x', 'y', label=f'measured - {cond}', 
                   data=plot_df[(plot_df['type'] == 'measured') & (plot_df['condition'] == cond)], 
                   color=color, s=5)

    # Add plot titles and labels
    ax.set_title(f'Posterior distribution - {prot}', fontsize=10, weight='bold')
    ax.set_xlabel('Temperature [ºC]', fontsize=10, weight='bold')
    ax.set_ylabel('Scaled intensity', fontsize=10, weight='bold')
    ax.legend(loc='lower left', fontsize='small')

# Function to plot covariance matrix directly into an axis
def plot_covariance_matrix_sns(ax, covariance_matrix, title="Covariance Matrix"):
    sns.heatmap(covariance_matrix.numpy(), cmap="coolwarm", annot=False, fmt=".1f", square=True, ax=ax)
    ax.set_title(title, fontsize=10, weight='bold')
    ax.set_xlabel("Test Points", fontsize=10, weight='bold')
    ax.set_ylabel("Test Points", fontsize=10, weight='bold')

def plot_covariance_matrix_imshow(ax, covariance_matrix, title="Covariance Matrix"):
    im = ax.imshow(covariance_matrix.numpy(), cmap="coolwarm", vmin=0, vmax=1)
    ax.set_title(title, fontsize=10, weight='bold')
    ax.set_xlabel("Test Points", fontsize=10, weight='bold')
    ax.set_ylabel("Test Points", fontsize=10, weight='bold')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Covariance Value', fontsize=10, weight='bold')
    
    # Adjust colorbar tick label sizes
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(8)


# Function to normalize the covariance matrix
def normalize_covariance_matrix(covariance_matrix):
    min_val = torch.min(covariance_matrix)
    max_val = torch.max(covariance_matrix)
    normalized_matrix = (covariance_matrix - min_val) / (max_val - min_val)
    return normalized_matrix

# Convert tensors to DataFrames
def cov_matrices_to_df(covariance_matrices):

    # prepare data for ataFrame Creation
    covariance_dfs = [(key, value.numpy()) for key, value in covariance_matrices.items()]

    # create a DataFrame with two columns: 'Key' and 'Tensor'
    covariance_df = pd.DataFrame(covariance_dfs, columns=['key', 'covariance'])
    
    # extract uniqueID and condition:
    covariance_df[['uniqueID', 'condition']] = covariance_df['key'].apply(lambda x: pd.Series(x.rsplit('_', 1)))
    
    # normalize covariance
    covariance_df['normalized_covariance'] = covariance_df['covariance'].apply(lambda x: normalize_covariance_matrix(torch.tensor(x)).numpy())
    
    # select relevant columns for final DataFrame:
    covariance_norm_df = covariance_df[['uniqueID', 'condition', 'normalized_covariance']]
    
    # return dataframe
    return covariance_norm_df

def combine_mll_dfs(mll_full, mll_joint):

    # Pivot the DataFrame
    df_full = mll_full.copy()
    df_pivot_full = df_full.pivot(index='protein', columns='condition', values='mll')

    # Reset the index to make 'protein' a column
    df_pivot_full = df_pivot_full.reset_index()

    # Remove the multi-level column name
    df_pivot_full.columns.name = None

    # Get mlls for joint model
    df_joint = mll_joint.copy()
    
    # Map mll of joint model to mll of full model
    df_pivot_full.loc[:, 'joint'] = df_pivot_full.loc[:,'protein'].map(df_joint.drop_duplicates('protein').set_index('protein')['mll'])
    
    # return combined dataframe
    return df_pivot_full

# compute likelihood ratios (LR)
def compute_likelihood_ratio(mll_full_df, mll_joint_df):
    mll_full = mll_full_df.copy()
    mll_full['mll_sum'] = mll_full.groupby(['protein'])['mll'].transform('sum')
    mll_full = mll_full.drop_duplicates(subset=['protein'])[['protein', 'mll_sum']]
    LR_df = pd.DataFrame({'protein' : mll_full['protein']})
    LR_df["mll_full"] = mll_full['mll_sum']
    LR_df.loc[:, "mll_joint"] = LR_df.loc[:,'protein'].map(mll_joint_df.drop_duplicates('protein').set_index('protein')['mll'])
    LR_df["lr"] = -2 * (LR_df['mll_joint']-LR_df['mll_full'])
    return(LR_df)

# compare LR to approximated null distribution and identify extreme values
'Adopted from: LeSueur, Cecile, Magnus Rattray, and Mikhail Savitski, BioRxiv (2024)'
def compute_sum_extreme(experimental_likelihood_ratio, approximated_likelihood_ratio):
    return sum(1 for ratio in approximated_likelihood_ratio if ratio >= experimental_likelihood_ratio)

def print_bold(string):
    display(Markdown(f"**{string}**"))

def benjamini_hochberg(p_values):
    """
    Implements the Benjamini-Hochberg procedure for controlling the false discovery rate.
    
    Parameters:
    p_values (array-like): Array of p-values
    
    Returns:
    numpy array: Adjusted p-values
    """
    p_values = np.array(p_values)
    n = len(p_values)
    sorted_indices = np.argsort(p_values)
    sorted_p_values = p_values[sorted_indices]
    adjusted_p_values = np.zeros(n)
    
    for i in range(n):
        adjusted_p_values[i] = sorted_p_values[i] * n / (i + 1)
    
    # Ensure the p-values are monotonic increasing
    adjusted_p_values = np.minimum.accumulate(adjusted_p_values[::-1])[::-1]
    
    # Place the adjusted p-values back into the original order
    adjusted_p_values_final = np.empty_like(adjusted_p_values)
    adjusted_p_values_final[sorted_indices] = adjusted_p_values
    
    # Cap the adjusted p-values at 1
    adjusted_p_values_final = np.minimum(adjusted_p_values_final, 1)
    
    return adjusted_p_values_final

# Function to update and display the checklist

def update_checklist(tasks, live_update=None, plot_data=None):
    clear_output(wait=True)
    
    # Display the checklist
    checklist = ""
    for task_info in tasks:
        if len(task_info) == 3:
            task, done, subprocesses = task_info
            checklist += f"- [{'x' if done else ' '}] {task}\n"
            for subtask, subdone in subprocesses:
                checklist += f"    - [{'x' if subdone else ' '}] {subtask}\n"
                if live_update and live_update['task'] == task and live_update['subtask'] == subtask:
                    checklist += f"        {live_update['message']}\n"
        elif len(task_info) == 2:
            task, done = task_info
            checklist += f"- [{'x' if done else ' '}] {task}\n"
            if live_update and live_update['task'] == task:
                checklist += f"    {live_update['message']}\n"
        else:
            raise ValueError(f"Task tuple must have either 2 or 3 elements, but got {len(task_info)}: {task_info}")
                
    display(Markdown(checklist))
    
    # Render the plot
    if plot_data and plot_data['task'] == live_update['task']:
        plt.figure(figsize=(6, 4))
        plt.plot(plot_data['x'], plot_data['y'], label=plot_data['label'])
        plt.xlabel(plot_data['xlabel'])
        plt.ylabel(plot_data['ylabel'])
        plt.title(plot_data['title'])
        plt.legend()
        plt.show()            