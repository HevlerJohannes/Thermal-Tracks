#! /usr/bin/env python3
import os
import numpy as np
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import pickle
from copy import deepcopy
import torch
import gc

# Import util functions
from utils import update_checklist
from utils import normalize_covariance_matrix
from utils import plot_posterior_distributions
from utils import plot_covariance_matrix_imshow
from utils import combine_mll_dfs


# Combine and save result files
def prepare_gp_results(result_dict_full, result_dict_null, parameters):

    plot_models = parameters['create_plots']
    if plot_models == True:

        # Update the checklist
        tasks = [
            ("1. Build and fit full model", True),
            ("2. Create a joint model and null dataset", True),
            ("3. Evaluate and predict models", True),
            ("4. Build and fit null model", True),
            ("5. Compute likelihood ratio test statistics", True),
            ("6. Combine and create result files", False, [
                ("Create plots for model comparison", False),
                ("Combine and save results", False)
            ])
            ]
    else:
        tasks = [
            ("1. Build and fit full model", True),
            ("2. Create a joint model and null dataset", True),
            ("3. Evaluate and predict models", True),
            ("4. Build and fit null model", True),
            ("5. Compute likelihood ratio test statistics", True),
            ("6. Combine and create result files", False, [
                ("Combine and save results", False)
            ])
            ]
        
    update_checklist(tasks)

    # Get predictions, covariance and protein ids
    prediction_result_df = result_dict_full['gp_result_df'].copy()
    covariance_matrices_full = result_dict_full['gp_full_covariance'].copy()
    covariance_matrices_joint = result_dict_full['gp_joint_covariance'].copy()

    proteins2test = result_dict_full['proteins'].copy()
    conditions = result_dict_full['conditions'].copy()
    output_path = parameters['result_dir']
    control = parameters['control_condition']
    perturbation = parameters['perturbation']

           
    if plot_models == True: # To do add code to plot covariance matrix!!!

        # Create directory for plots
        # Define the name of the new directory
        plot_dir = "result_plots"

        # Create the full path for the new directory
        plot_dir_path = os.path.join(output_path, plot_dir)

        # Create the directory if it doesn't exist
        if not os.path.exists(plot_dir_path):
            os.makedirs(plot_dir_path)

        # Initialize an iteration counter
        iteration_counter = 0
                        
        # Create plots as a single pdf page
        for prot in proteins2test:
            
            with PdfPages(f'{plot_dir_path}/result_plots_{prot}.pdf') as pdf:
                        
                # Increment the iteration counter
                iteration_counter += 1
                plot_df = prediction_result_df[prediction_result_df['uniqueID'] == prot]
                
                covariance_matrix_control = covariance_matrices_full[f"{prot}_{control}"]
                normalized_covariance_matrix_control = normalize_covariance_matrix(covariance_matrix_control)
                covariance_matrix_perturbation = covariance_matrices_full[f"{prot}_{perturbation}"]
                normalized_covariance_matrix_perturbation = normalize_covariance_matrix(covariance_matrix_perturbation)
                covariance_matrix_joint = covariance_matrices_joint[f"{prot}_joint"]
                normalized_covariance_matrix_joint = normalize_covariance_matrix(covariance_matrix_joint)

                fig, axs = plt.subplots(2, 2, figsize=(10, 8))

                # Plot the posterior distribution
                plot_posterior_distributions(axs[1, 1], prot, plot_df, conditions)

                # Plot the covariance matrices
                plot_covariance_matrix_imshow(axs[0, 0], normalized_covariance_matrix_control, title=f"{prot} - {control}")
                plot_covariance_matrix_imshow(axs[0, 1], normalized_covariance_matrix_perturbation, title=f"{prot} - {perturbation}")
                plot_covariance_matrix_imshow(axs[1, 0], normalized_covariance_matrix_joint, title=f"{prot} - joint")

                plt.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)
                # Explicitly delete large objects
                del plot_df, normalized_covariance_matrix_control, normalized_covariance_matrix_perturbation, normalized_covariance_matrix_joint 
                gc.collect()
                
                # Update checklist        
                live_update_message = f"(Plots generated: {iteration_counter}/{len(proteins2test)})"
                update_checklist(tasks, live_update={"task": "6. Combine and create result files", "subtask": "Create plots for model comparison", "message": live_update_message})
            
        live_update_message = f"(All plots generated and saved!)"
        update_checklist(tasks, live_update={"task": "6. Combine and create result files", "subtask": "Create plots for model comparison", "message": live_update_message})
        tasks[5][2][0] = ("Create plots for model comparison", True)
        update_checklist(tasks)

    # Save created models and likelihoods (full and joint)
    torch.save(result_dict_full["full_model_list"].state_dict(), f'{output_path}/gp_full_model_{parameters["control_condition"]}_{parameters["perturbation"]}.pth')
    torch.save(result_dict_full["full_likelihood_list"].state_dict(), f'{output_path}/gp_full_likelihood_{parameters["control_condition"]}_{parameters["perturbation"]}.pth')
    torch.save(result_dict_full["joint_model_list"].state_dict(), f'{output_path}/gp_joint_model_{parameters["control_condition"]}_{parameters["perturbation"]}.pth')
    torch.save(result_dict_full["joint_likelihood_list"].state_dict(), f'{output_path}/gp_joint_likelihood_{parameters["control_condition"]}_{parameters["perturbation"]}.pth')


    # Save result dictionaries as pickle file
    # List of keys you want to remove
    keys_to_remove = ['full_model_list', 'full_likelihood_list', 'full_state_dict_list', 'joint_model_list', 'joint_likelihood_list', 'joint_state_dict_list'] 
    # Remove the keys from the dictionary
    result_dict_filt = deepcopy(result_dict_full)
    for key in keys_to_remove:
        result_dict_filt.pop(key, None)

    with open(f'{output_path}/combined_gp_results_{parameters["control_condition"]}_{parameters["perturbation"]}.pkl', 'wb') as file:
            pickle.dump(result_dict_filt, file)

    # Save single result dataframes for easy access
    
    # Full model
    keys_full = ['exactgp_input', 'full_fit_parameters_df', 'lr_values_full_vs_joint','gp_result_df', 'gp_likelihood_statistics_df']
    file_names_full = ['model_input', 'loss_full_model', 'lr_values_full_model','predictions', 'lr_test_statistics']

    for i, j in zip(keys_full, file_names_full):
        result_dict_full[i].to_csv(f'{output_path}/{j}.csv', index = False)

    # combine mll values full model
    mll_values_combined_full = combine_mll_dfs(result_dict_full['full_mll_values'].copy(), result_dict_full['joint_mll_values'].copy()) 
    mll_values_combined_full.to_csv(f'{output_path}/mll_values_full_model.csv', index = False)
    # Null model
    keys_null = ['exactgp_input', 'full_fit_parameters_df', 'lr_values_full_vs_joint']
    file_names_null = ['sampled_null_dataset', 'loss_null_model', 'lr_values_null_model']

    for i, j in zip(keys_null, file_names_null):
        result_dict_null[i].to_csv(f'{output_path}/{j}.csv', index = False) 

    # combine mll values null model
    mll_values_combined_null = combine_mll_dfs(result_dict_null['full_mll_values'].copy(), result_dict_null['joint_mll_values'].copy())        
    mll_values_combined_null.to_csv(f'{output_path}/mll_values_null_model.csv', index = False)

    # Update checklist
    if plot_models == True:
        tasks[5][2][1] = ("Combine and save results", True)
        update_checklist(tasks)
        
        update_checklist(tasks)
    else:
        tasks[5][2][0] = ("Combine and save results", True)
        update_checklist(tasks)
       
    tasks[5] = ("6. Combine and create result files", True)
    update_checklist(tasks)             