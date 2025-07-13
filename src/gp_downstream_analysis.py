#! /usr/bin/env python3
import pickle
import numpy as np
import pandas as pd
from copy import deepcopy
import matplotlib.pyplot as plt
import torch
import os
import gc
from IPython.display import Markdown, display, clear_output
import multiprocessing as mp
import matplotlib
import warnings
from matplotlib.colors import ListedColormap, Normalize, to_rgba
from matplotlib.cm import ScalarMappable
import seaborn as sns

# Ignore all warnings
warnings.filterwarnings("ignore") 

# Set figure style
rc = {'axes.spines.right':False, 'axes.spines.top':False, 'axes.edgecolor':'.2', 'svg.fonttype':'none'}
sns.set_theme(style='ticks', context='paper', font='arial', rc=rc)

####################################################################
#                  Import gp result pickle file                    #
####################################################################
def load_selected_keys_from_pickle(file_path, selected_keys):
    
    selected_data = {}
    
    with open(file_path, 'rb') as file:
        while True:
            try:
                # Load the next object from the pickle file
                data = pickle.load(file)
                for key in selected_keys:
                    if key in data:
                        selected_data[key] = data[key]
            except EOFError:
                # End of file reached
                break
            except KeyError:
                # Continue if the key is not found in the current data object
                continue

    return selected_data

# Usage
def import_gp_result_dict(parameters):

    # Get path for gp result pickle file
    result_file_path = parameters['result_dir']
    for file in os.listdir(result_file_path):
        if file.endswith(".pkl"):
            pickle_file_path = (os.path.join(result_file_path, file))    
    
    keys_to_load = ['full_mll_values', 'proteins','conditions', 'output_path', 'joint_mll_values','gp_result_df', 'gp_full_covariance', 'gp_joint_covariance', 'gp_likelihood_statistics_df']
    selected_data = load_selected_keys_from_pickle(pickle_file_path, keys_to_load)

    return(selected_data)


####################################################################
#                Generate result plots (Effect size)               #
####################################################################
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
def plot_covariance_matrix_from_df(ax, covariance_matrix, title="Covariance Matrix"):
    # Extract the normalized_covariance matrix from the DataFrame
    if not covariance_matrix.empty:
        # Assuming the first element is representative
        normalized_covariance = covariance_matrix['normalized_covariance'].iloc[0]

        # Convert the list of lists to a NumPy array
        if isinstance(normalized_covariance, list):
            normalized_covariance = np.array(normalized_covariance)
    im = ax.imshow(normalized_covariance, cmap="coolwarm", vmin=0, vmax=1)
    ax.set_title(title, fontsize=10, weight='bold')
    ax.set_xlabel("Test Points", fontsize=10, weight='bold')
    ax.set_ylabel("Test Points", fontsize=10, weight='bold')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Covariance Value', fontsize=10, weight='bold')
    
    # Adjust colorbar tick label sizes
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(10)

# function to generate result plot
def generate_result_plot(prot, prediction_result_df, covariance_matrices, conditions, control, perturbation, plot_dir_path):
    
    prot_pred_df = prediction_result_df[prediction_result_df['uniqueID'] == prot]
    prot_covariance_matrix = covariance_matrices[covariance_matrices['uniqueID'] == prot]

    # Create figure
    fig, axs = plt.subplots(2, 2, figsize=(15, 12))
    try:

        # Plot the posterior distribution
        plot_posterior_distributions(axs[1, 1], prot, prot_pred_df, conditions)
        # Plot the covariance matrices
        plot_covariance_matrix_from_df(axs[0, 0], prot_covariance_matrix[prot_covariance_matrix['condition'] == control], title=f"{prot} - {control}")
        plot_covariance_matrix_from_df(axs[0, 1], prot_covariance_matrix[prot_covariance_matrix['condition'] == perturbation], title=f"{prot} - {perturbation}")
        plot_covariance_matrix_from_df(axs[1, 0], prot_covariance_matrix[prot_covariance_matrix['condition'] == 'joint'], title=f"{prot} - joint")
        
        # Save the figure directly as a PDF
        fig.savefig(f'{plot_dir_path}/result_plots_{prot}.pdf', format='pdf')
        #print(f"Saved plot for: {prot} at {plot_dir_path}/result_plots_{prot}.pdf")  # Debugging statement
    except Exception as e:
        print(f"Error generating plot for {prot}: {e}")
    finally:
        plt.clf()
        plt.close(fig)
        # Explicitly delete large objects
        del prot_pred_df, prot_covariance_matrix
        gc.collect()

def generate_result_plots(result_dict, parameters):
    matplotlib.use('Agg') 
    # Print checklist
    tasks = [("Create plots for model comparison", False)]  
    update_checklist(tasks)

    prediction_result_df = result_dict['gp_result_df'].copy()
    covariance_matrices_full = cov_matrices_to_df(result_dict['gp_full_covariance'].copy())
    covariance_matrices_joint = cov_matrices_to_df(result_dict['gp_joint_covariance'].copy())
    covariance_matrices = pd.concat([covariance_matrices_full, covariance_matrices_joint])
    proteins2test = result_dict['proteins'].copy()
    conditions = result_dict['conditions'].copy()
    output_path = parameters['result_dir']
    control = parameters['control_condition']
    perturbation = parameters['perturbation']

    # Create directory for plots
    # Define the name of the new directory
    plot_dir = "result_plots"
    # Create the full path for the new directory
    plot_dir_path = os.path.join(output_path, plot_dir)
    # Create the directory if it doesn't exist
    if not os.path.exists(plot_dir_path):
        os.makedirs(plot_dir_path)
       
    iteration_counter = 0
    
    for prot in proteins2test:
        iteration_counter += 1
        generate_result_plot(prot, prediction_result_df, covariance_matrices, conditions, control, perturbation, plot_dir_path)
        # Update checklist        
        live_update_message = f"(Plots generated: {iteration_counter}/{len(proteins2test)})"
        update_checklist(tasks, live_update={"task": "Create plots for model comparison", "message": live_update_message})
            
    live_update_message = f"(All plots generated and saved!)"
    update_checklist(tasks, live_update={"task": "Create plots for model comparison", "message": live_update_message})
    tasks = [("Create plots for model comparison", True)] 
    update_checklist(tasks)
     

####################################################################
#          Calculate area between curves (Effect size)             #
####################################################################
def calculate_difference(control_lower, control_upper, perturbation_lower, perturbation_upper, control_y, perturbation_y):
    differences = np.zeros(control_y.size, dtype=np.float64)
    
    for i in range(control_y.size):
        difference = 0.0

        # Calculate signed difference based on the given conditions
        if (
            (perturbation_lower[i] > control_lower[i] and perturbation_lower[i] < control_upper[i]) or
            (perturbation_upper[i] > control_lower[i] and perturbation_upper[i] < control_upper[i]) or
            (control_lower[i] > perturbation_lower[i] and control_lower[i] < perturbation_upper[i]) or
            (control_upper[i] > perturbation_lower[i] and control_upper[i] < perturbation_upper[i])
        ):
            difference = 0.0
        else:
            difference = perturbation_y[i] - control_y[i]  # Calculate signed difference

        differences[i] = difference  # Store the difference
    
    return differences

def calculate_differences(group, control, perturbation, strict_auc_calc = [False, True]):
    control_data = group[group['condition'] == control]
    perturbation_data = group[group['condition'] == perturbation]

    if strict_auc_calc == False:
        control_lower = control_data['conf_lower'].values
        control_upper = control_data['conf_upper'].values
        perturbation_lower = perturbation_data['conf_lower'].values
        perturbation_upper = perturbation_data['conf_upper'].values
    else:
        control_lower = control_data['conflik_lower'].values
        control_upper = control_data['conflik_upper'].values
        perturbation_lower = perturbation_data['conflik_lower'].values
        perturbation_upper = perturbation_data['conflik_upper'].values

    control_y = control_data['y'].values
    perturbation_y = perturbation_data['y'].values

    # Calculate differences using the CPU function
    differences = calculate_difference(control_lower, control_upper, perturbation_lower, perturbation_upper, control_y, perturbation_y)

    # Calculate the sum of differences
    total_difference = np.sum(differences)
    total_absolute_difference = np.sum(np.abs(differences))

    # Return both the individual differences and the total summed difference
    return differences, total_difference, total_absolute_difference

def compute_effect_size(result_dict, parameters, strict_auc = False):

    assert strict_auc == False or strict_auc == True, 'strict_auc needs to be boolean'
    
    # Get some parameters
    perturbation= parameters['perturbation']
    control  = parameters['control_condition']

    # Get gp result dataframes
    predictions_df = result_dict['gp_result_df'].copy()
    fits_df = predictions_df[predictions_df['type'] != 'measured']
    pval_stats_df = result_dict['gp_likelihood_statistics_df'].copy()

    protein_groups = fits_df.groupby('uniqueID')
    
    # Create an empty list to store the results
    auc_results = []
    
    # Loop over each protein group
    for protein, group in protein_groups:
        differences, total_diff, total_absolute_difference = calculate_differences(group, control, perturbation, strict_auc_calc = strict_auc)

        auc_results.append({
            'protein': protein,
            'effect_size': total_diff,
            'absolute_effect_size': total_absolute_difference,
            'auc': differences.tolist()  # Store differences as a list
        })
    # Convert results to a DataFrame
    auc_results_df = pd.DataFrame(auc_results)

    # Map effect size to LR stats df
    pval_stats_df.loc[:,'effect_size'] = pval_stats_df.loc[:,'protein'].map(auc_results_df.drop_duplicates('protein').set_index('protein')['effect_size'])
    pval_stats_df.loc[:,'absolute_effect_size'] = pval_stats_df.loc[:,'protein'].map(auc_results_df.drop_duplicates('protein').set_index('protein')['absolute_effect_size'])

    # Save results to result dict
    auc_result_dict = deepcopy(result_dict)
    auc_result_dict.update({
            "effect_size_df" : auc_results_df,
            "gp_likelihood_statistics_df" : pval_stats_df})
    
    return auc_result_dict

####################################################################
#                   Create and update a checklist                  #
####################################################################
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

####################################################################
#                Plotting Effect size vs pValue                    #
####################################################################

def cluster_dataframe_column(df, column_name, new_column_name='cluster_group', zero_range=1e-6):
    """
    Clusters the values in a specified column of a DataFrame into 11 groups:
    - Five groups for positive values.
    - One group for values around zero (no effect).
    - Five groups for negative values.

    :param df: Input DataFrame.
    :param column_name: Name of the column to be clustered.
    :param new_column_name: Name of the new column to store the groups. Default is 'group'.
    :param zero_range: Threshold for defining values around zero. Default is 1e-6.
    :return: DataFrame with a new column containing the assigned groups.
    """
    # Extract column values
    values = df[column_name]
    
    # Separate positive, negative, and zero values
    pos_values = values[values > zero_range]
    neg_values = values[values < -zero_range]

    # Handle edge cases where max or min values might not exist
    max_pos = max(pos_values) if not pos_values.empty else zero_range
    min_neg = min(neg_values) if not neg_values.empty else -zero_range

    # Calculate bin edges for positive and negative values
    pos_bins = np.linspace(0, max_pos, 6)  # 5 bins for positive values
    neg_bins = np.linspace(min_neg, 0, 6)  # 5 bins for negative values

    # Combine bins: negative bins + zero bin + positive bins
    bins = np.concatenate((neg_bins[:-1], [-zero_range, zero_range], pos_bins[1:]))

    # Define custom labels
    labels = [
        'destabilization_5', 'destabilization_4', 'destabilization_3', 
        'destabilization_2', 'destabilization_1', 'no effect',
        'stabilization_1', 'stabilization_2', 'stabilization_3', 
        'stabilization_4', 'stabilization_5'
    ]

    # Apply clustering to the column
    df[new_column_name] = pd.cut(values, bins=bins, labels=labels, include_lowest=True)
    
    return df

# Creating a circular Manhattan plot with annotations for surface and N-glycosylated proteins
def plot_manhattan_nglyco(result_dict, parameters, annotation_dict = {}, pvalue=['adj_pValue', 'pValue'], figsize = [15,15], fill_to_end = [True, False]):

    if pvalue not in ['adj_pValue', 'pValue']:
        raise ValueError("p-value must be 'adj_pValue' or 'pValue'")
    
    # Get data
    gp_stats_df = result_dict['gp_likelihood_statistics_df'].copy()
    glyco_proteins= annotation_dict['glyco_proteins'].copy()
    cs_proteins = annotation_dict['cs_proteins'].copy()

    # Get information about condition
    perturbation = parameters['perturbation']
    control = parameters['control_condition']
    
    # Format data
    gp_stats_df['accession'] = gp_stats_df["protein"].str.split('-').str[0]
    gp_stats_df['gene_name'] = gp_stats_df["protein"].str.split('-').str[1]
    gp_stats_df['gene_name'] = gp_stats_df["gene_name"].str.extract('.*\((.*?)\).*', expand = True)
        
    # Add annotations
    gp_stats_df.loc[:,'glycan_type'] = gp_stats_df.loc[:,'accession'].map(glyco_proteins.drop_duplicates('accession').set_index('accession')['glycan_type'])
    gp_stats_df.loc[:,'glycan_type'] = np.where(gp_stats_df.loc[:,'glycan_type'].isna(), 'none', gp_stats_df.loc[:,'glycan_type'])
    gp_stats_df.loc[:,'cell_surface'] = np.where(gp_stats_df.loc[:,'accession'].isin(cs_proteins['ID_link']), 'yes', 'no')  
    gp_stats_df['protein_type'] = np.where((gp_stats_df['glycan_type'] == 'nglycan') & (gp_stats_df['cell_surface'] == 'yes'), 'cs_nglyco', 'other')
    gp_stats_df['protein_type'] = np.where((gp_stats_df['glycan_type'] == 'nglycan') & (gp_stats_df['cell_surface'] == 'no'), 'nglyco', gp_stats_df['protein_type'])
    
    # Prepare data for the plot
    ## Transform p-values to -log10(p-value)
    log10_p_values = -np.log10(gp_stats_df[pvalue])
    
    ## Cluster effect sizes into groups
    gp_stats_df = cluster_dataframe_column(gp_stats_df, column_name='effect_size')
    
    # Transform p-values to -log10(p-value)
    gp_stats_df['-log10_p_values'] = -np.log10(gp_stats_df[pvalue])
    
    # Set up figure and axis
    # Create a polar plot
    fig, ax = plt.subplots(figsize=(figsize[0], figsize[1]), subplot_kw={'projection': 'polar'})

    # Define background colors for each group (adjust colors if needed)
    background_colors = ['#053061', '#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#bdbdbd', 
                        '#fddbc7', '#f4a582', '#d6604d', '#b2182b', '#67001f']
    
    # Rotate the plot 45 degrees to the left
    ax.set_theta_offset(np.pi / 2)
    # Convert background colors to RGBA with alpha = 0.4
    background_colors_rgba = [to_rgba(color, alpha=0.4) for color in background_colors]

    # Create a colormap from the background colors with alpha
    cmap = ListedColormap(background_colors_rgba)
    norm = Normalize(vmin=0, vmax=len(background_colors_rgba) - 1)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    # Original group sizes
    group_sizes = [0.1, 0.2, 0.3, 0.5, 0.8, 0.4, 0.8, 0.5, 0.3, 0.2, 0.1]

    # Normalize the group sizes
    total_sum = sum(group_sizes)
    normalized_group_sizes = [x / total_sum for x in group_sizes]

    # Verify that the sum of normalized group sizes is 1
    assert np.isclose(sum(normalized_group_sizes), 1), "Group sizes must sum up to 1."
    
    # Calculate angle ranges based on normalized sizes
    start_angle = 0
    gap = np.pi / 180  # make gap (1 degree in radians)
    angle_ranges = []

    # Calculate angle range for each group based on normalized sizes
    for size in normalized_group_sizes:
        group_angle = (2 * np.pi * size) - gap
        end_angle = start_angle + group_angle
        angle_ranges.append((start_angle, end_angle))
        start_angle = end_angle + gap  # Update start angle for next group

    # Fill background for each group
    for i, (start_angle, end_angle) in enumerate(angle_ranges):
        theta = np.linspace(start_angle, end_angle, 100)
        r = np.full_like(theta, max(gp_stats_df['-log10_p_values']) * 1.2)  # Slightly larger radius to cover the area
        ax.fill_between(theta, 0, r, color=background_colors[i % len(background_colors)], alpha=0.4, zorder=0)

    # Define color mapping function
    def get_color(row):
        if row['-log10_p_values'] > -np.log10(0.05):
            if row['protein_type'] == 'cs_nglyco':
                return '#7570b3'  # Significant cell surface glycoproteins
            elif row['protein_type'] == 'nglyco':
                return '#1b9e77'  # Significant glycoproteins
            else:
                return 'black'  # Significant other proteins
        else:
            return '#bababa'  # Non-significant proteins

    # Define alpha setting function
    def get_alpha(row):
        if row['protein_type'] in ['cs_nglyco', 'nglyco'] and row['-log10_p_values'] > -np.log10(0.05):
            return 1.0  # Full opacity for cell surface proteins and glycoproteins
        else:
            return 0.3  # Lower opacity for other proteins

    # Define scatter size function
    def get_size(row):
        if row['protein_type'] in ['cs_nglyco', 'nglyco'] and row['-log10_p_values'] > -np.log10(0.05):
            return 72  # Bigger scatter size for significant cell surface- and glyco proteins
        else:
            return 36  # Standard scatter size for other proteins    

    # Create a circular plot with each group in its designated section, colored by protein type and significance
    for i, label in enumerate(gp_stats_df['cluster_group'].cat.categories):
        group_mask = gp_stats_df['cluster_group'] == label
        group_r = gp_stats_df.loc[group_mask, '-log10_p_values']
        group_size = group_r.size  # Get the number of points in the current group
        

        # Check if group size is 1, center the point if true
        if group_size < 4:
            # Set angles evenly spaced within the angle range for the group
            group_theta = np.linspace(angle_ranges[i][0], angle_ranges[i][1], group_size + 2)[1:-1]
        else:
            # Calculate angles dynamically for the current group based on its size
            group_theta = np.linspace(angle_ranges[i][0], angle_ranges[i][1], group_size, endpoint=False)
    
    
        
        # Determine colors and alphas based on log_p_values and protein type
        colors = gp_stats_df.loc[group_mask].apply(get_color, axis=1)
        alphas = gp_stats_df.loc[group_mask].apply(get_alpha, axis=1)
        scatter_size = gp_stats_df.loc[group_mask].apply(get_size, axis=1)
        
        if group_size == 0:
            continue
        else:
            # Determine zorder based on color
            if any(color in ['#7570b3', '#1b9e77'] for color in colors):
                zorder = 2
            else:
                zorder = 1  # Default zorder

            # Plot with variable alpha and zorder
            ax.scatter(group_theta, group_r, c=colors, s = scatter_size, edgecolors='face', alpha=alphas, label=label, zorder=zorder)

    # Set title and labels
    ax.set_title(f'Protein thermal stability changes {perturbation} vs. {control}', va='bottom', fontsize=12, weight='bold')
    ax.set_xticklabels([])  # remove x-tick labels for a cleaner plot
    ax.set_yticklabels([])  # remove y-tick labels for a cleaner plot

    # Set y-axis limits
    if fill_to_end == True:
        ax.set_ylim(0, max(log10_p_values) * 1.2)
    else:
        None

    # Adding significance threshold line after setting axis limits
    significance_threshold = -np.log10(0.05)
    theta = np.linspace(0, 2 * np.pi, 100)
    r = np.full_like(theta, significance_threshold)
    ax.plot(theta, r, color='black', linestyle='dashed', linewidth=1, zorder=2)

    # Adding custom legend
    pval_label = 'adj. p-value' if pvalue == 'adj_pValue' else 'p-value'
    legend_handles = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#7570b3', markersize=12, label=f'significant cell surface glycoprotein ({pval_label} < 0.05)', alpha=1),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#1b9e77', markersize=12, label=f'significant glycoprotein ({pval_label} < 0.05)', alpha=1),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=8, label=f'significant other protein ({pval_label} < 0.05)', alpha=0.3),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#bababa', markersize=8, label=f'non-significant protein ({pval_label} > 0.05)', alpha=0.3),
        plt.Line2D([0], [0], color='black', linestyle='dashed', linewidth=1, label=f'Significance threshold ({pval_label} = 0.05)')
    ]

    ax.legend(handles=legend_handles, loc='lower right', bbox_to_anchor=(1.3, -0.1))

    # Add the colorbar legend for stabilization gradient
    cbar = fig.colorbar(sm, ax=ax, orientation='horizontal', pad=0.05, fraction=0.02)
    cbar.set_label('Stabilization Gradient\nmost destabilized \u2194 most stabilized', rotation=0, labelpad=5, fontsize=8, weight='bold')
    cbar.set_ticks([])

    # Display the plot
    output_path = parameters['result_dir']
    plt.savefig(f'{output_path}/thermal_stability_changes_{perturbation}_vs_{control}_with_annotations.pdf', dpi = 400)
    plt.savefig(f'{output_path}/thermal_stability_changes_{perturbation}_vs_{control}_with_annotations.svg')
    print(f'A circular Manhattan plot has been created: {output_path}/thermal_stability_changes_{perturbation}_vs_{control}_with_annotations')
    plt.show()
    plt.close()
    
    # Save results to result dict
    plot_result_dict = deepcopy(result_dict)
    plot_result_dict.update({
            "gp_likelihood_statistics_df" : gp_stats_df})
    
    return plot_result_dict

# Creating a circular Manhattan plot
def plot_manhattan(result_dict, parameters, pvalue=['adj_pValue', 'pValue'], figsize = [15,15], fill_to_end = [True, False]):

    if pvalue not in ['adj_pValue', 'pValue']:
        raise ValueError("p-value must be 'adj_pValue' or 'pValue'")
    
    # Get data
    gp_stats_df = result_dict['gp_likelihood_statistics_df'].copy()
    
    # Get information about condition
    perturbation = parameters['perturbation']
    control = parameters['control_condition']
    
    # Format data
    gp_stats_df['accession'] = gp_stats_df["protein"].str.split('-').str[0]
    gp_stats_df['gene_name'] = gp_stats_df["protein"].str.split('-').str[1]
        
    # Prepare data for the plot
    ## Transform p-values to -log10(p-value)
    log10_p_values = -np.log10(gp_stats_df[pvalue])
    
    ## Cluster effect sizes into groups
    gp_stats_df = cluster_dataframe_column(gp_stats_df, column_name='effect_size')
    
    ## Transform p-values to -log10(p-value)
    gp_stats_df['-log10_p_values'] = -np.log10(gp_stats_df[pvalue])
    
    # Set up figure and axis
    # Create a polar plot
    fig, ax = plt.subplots(figsize=(figsize[0], figsize[1]), subplot_kw={'projection': 'polar'})

    # Define background colors for each group (adjust colors if needed)
    background_colors = ['#053061', '#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#bdbdbd', 
                        '#fddbc7', '#f4a582', '#d6604d', '#b2182b', '#67001f']
    
    # Rotate the plot 45 degrees to the left
    ax.set_theta_offset(np.pi / 2)
    # Convert background colors to RGBA with alpha = 0.4
    background_colors_rgba = [to_rgba(color, alpha=0.4) for color in background_colors]

    # Create a colormap from the background colors with alpha
    cmap = ListedColormap(background_colors_rgba)
    norm = Normalize(vmin=0, vmax=len(background_colors_rgba) - 1)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    # Original group sizes
    group_sizes = [0.1, 0.2, 0.3, 0.5, 0.8, 0.4, 0.8, 0.5, 0.3, 0.2, 0.1]

    # Normalize the group sizes
    total_sum = sum(group_sizes)
    normalized_group_sizes = [x / total_sum for x in group_sizes]

    # Verify that the sum of normalized group sizes is 1
    assert np.isclose(sum(normalized_group_sizes), 1), "Group sizes must sum up to 1."
    
    # Calculate angle ranges based on normalized sizes
    start_angle = 0
    gap = np.pi / 180  # example gap (1 degree in radians)
    angle_ranges = []

    # Calculate angle range for each group based on normalized sizes
    for size in normalized_group_sizes:
        group_angle = (2 * np.pi * size) - gap
        end_angle = start_angle + group_angle
        angle_ranges.append((start_angle, end_angle))
        start_angle = end_angle + gap  # Update start angle for next group

    # Fill background for each group
    for i, (start_angle, end_angle) in enumerate(angle_ranges):
        theta = np.linspace(start_angle, end_angle, 100)
        r = np.full_like(theta, max(gp_stats_df['-log10_p_values']) * 1.2)  # Slightly larger radius to cover the area
        ax.fill_between(theta, 0, r, color=background_colors[i % len(background_colors)], alpha=0.4, zorder=0)

    # Define color mapping function
    def get_color(row):
        if row['-log10_p_values'] >= -np.log10(0.05):
            return 'black'  # Significant proteins
            
        else:
            return '#bababa'  # Non-significant proteins

    # Define alpha setting function
    def get_alpha(row):
        if  row['-log10_p_values'] >= -np.log10(0.05):
            return 1.0  # Full opacity for significant proteins
        else:
            return 0.3  # Lower opacity for other proteins

    # Create a circular plot with each group in its designated section, colored by protein type and significance
    for i, label in enumerate(gp_stats_df['cluster_group'].cat.categories):
        group_mask = gp_stats_df['cluster_group'] == label
        group_r = gp_stats_df.loc[group_mask, '-log10_p_values']
        group_size = group_r.size  # Get the number of points in the current group

        # Check if group size is one, center the point if true
        if group_size < 4:
            # Set angles evenly spaced within the angle range for the group
            group_theta = np.linspace(angle_ranges[i][0], angle_ranges[i][1], group_size + 2)[1:-1]
        else:
            # Calculate angles dynamically for the current group based on its size
            group_theta = np.linspace(angle_ranges[i][0], angle_ranges[i][1], group_size, endpoint=False)
    
    
        
        # Determine colors and alphas based on log_p_values and protein type
        colors = gp_stats_df.loc[group_mask].apply(get_color, axis=1)
        alphas = gp_stats_df.loc[group_mask].apply(get_alpha, axis=1)
        
        if group_size == 0:
            continue
        else:
            # Plot with variable alpha
            ax.scatter(group_theta, group_r, c=colors, edgecolors='face', alpha=alphas, label=label)

    # Set title and labels
    ax.set_title(f'Protein thermal stability changes {perturbation} vs. {control}', va='bottom', fontsize=12, weight='bold')
    ax.set_xticklabels([])
    ax.set_yticklabels([])  # Remove y-tick labels for a cleaner plot

    # Set y-axis limits (defines whether background color goes until edge of the circle or keeps a small gap)
    if fill_to_end == True:
        ax.set_ylim(0, max(log10_p_values) * 1.2)
    else:
        None

    # Adding significance threshold line after setting axis limits
    significance_threshold = -np.log10(0.05)
    theta = np.linspace(0, 2 * np.pi, 100)
    r = np.full_like(theta, significance_threshold)
    ax.plot(theta, r, color='black', linestyle='dashed', linewidth=1, zorder=2)

    # Adding custom legend
    pval_label = 'adj. p-value' if pvalue == 'adj_pValue' else 'p-value'
    legend_handles = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='black', markersize=10, label=f'significant protein ({pval_label} < 0.05)', alpha=1),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#bababa', markersize=10, label=f'non-significant protein ({pval_label} > 0.05)', alpha=0.3),
        plt.Line2D([0], [0], color='black', linestyle='dashed', linewidth=1, label=f'Significance threshold ({pval_label} = 0.05)')
    ]

    ax.legend(handles=legend_handles, loc='lower right', bbox_to_anchor=(1.3, -0.1))

    # Add the colorbar legend for stabilization gradient
    cbar = fig.colorbar(sm, ax=ax, orientation='horizontal', pad=0.05, fraction=0.02)
    cbar.set_label('Stabilization Gradient\nmost destabilized \u2194 most stabilized', rotation=0, labelpad=5, fontsize=8, weight='bold')
    cbar.set_ticks([])

    # Display the plot
    output_path = parameters['result_dir']
    plt.savefig(f'{output_path}/thermal_stability_changes_{perturbation}_vs_{control}.pdf', dpi = 400)
    plt.savefig(f'{output_path}/thermal_stability_changes_{perturbation}_vs_{control}.svg')
    print(f'A circular Manhattan plot has been created: {output_path}/thermal_stability_changes_{perturbation}_vs_{control}')
    plt.show()
    plt.close()

    # Save results to result dict
    plot_result_dict = deepcopy(result_dict)
    plot_result_dict.update({
            "gp_likelihood_statistics_df" : gp_stats_df})

    return plot_result_dict

####################################################################
#                     Rerun LR statistics                          #
####################################################################

def run_lr_test(mll_full, lr_full, lr_null, mll_filter = None, conditions = ['condition_a', 'condition_b']):
    
    # Get proteins with good fits
    if mll_filter != None:
        proteins_filt = mll_full[(mll_full[conditions[0]] < mll_filter) | (mll_full[conditions[1]] < mll_filter)]['protein']
        lr_df_full = lr_full[~lr_full['protein'].isin(proteins_filt)]
        lr_df_null = lr_null[~lr_null['protein'].isin(proteins_filt)].copy()
        #lr_df_null = lr_null.copy()
        
    else:
        lr_df_full = lr_full.copy()
        lr_df_null = lr_null.copy()

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
        
    return(p_val) 

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