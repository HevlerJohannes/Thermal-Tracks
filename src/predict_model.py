#! /usr/bin/env python3
import numpy as np
import pandas as pd
import torch
import gpytorch
import pandas as pd
from copy import deepcopy


# Import util functions
from utils import update_checklist, is_notebook


def predict_and_evaluate_batched(result_dict, parameters, null_dataset=[True, False]):

    # Update the checklist
    if is_notebook():
        if null_dataset == False:
            tasks = [
                ("1. Build and fit full model", True),
                ("2. Create a joint model and null dataset", True),
                ("3. Evaluate and predict models", False, [
                    ("Full model", False),
                    ("Joint model", False),
                    ("Combine and save predictions", False)
                ]),
                ("4. Build and fit null model", False),
                ("5. Compute likelihood ratio test statistics", False),
                ("6. Combine and create result files", False)
            ]
        else:
            tasks = [
                ("1. Build and fit full model", True),
                ("2. Create a joint model and null dataset", True),
                ("3. Evaluate and predict full and joint models", True),
                ("4. Build and fit null model", False, [
                    ("Train null model", True),
                    ("Create a joint model", True),
                    ("Evaluate and predict null model", False)
                ]),
                ("5. Compute likelihood ratio test statistics", False),
                ("6. Combine and create result files", False)
            ]
    
        update_checklist(tasks)

    # Get model and likelihood lists for full and joint models
    full_model_list = result_dict['full_model_list']  # List of (model, likelihood, metadata) tuples
    metadata_list = result_dict['full_model_metadata']  # Flat list of metadata dicts
    joint_model = result_dict['joint_model_list']
    joint_likelihood = result_dict['joint_likelihood_list']

    # Get additional model information
    n_total_tasks = result_dict['n_total_tasks']
    proteins2test = result_dict['proteins']
    n_prot = result_dict['n_proteins']
    conds = result_dict['conditions']
    n_cond = result_dict['n_conditions']
    tpptr_df = result_dict['exactgp_input'].copy()
    
    # Get output and perturbation information
    output_path = str(parameters['result_dir'])
    pert = parameters['perturbation']
    
    # Get device
    device_str = result_dict.get('device', 'cpu')
    device = torch.device(device_str)
    dtype = parameters['dtype']

    # Define number of prediction points
    n_predict = parameters['n_predictions']

    ############################## MODEL EVALUATE and PREDICT ##############################
    # Predictions for the full model
    
    # Set all models into eval mode
    for model, likelihood, _ in full_model_list:
        model.eval()
        likelihood.eval()

    # Make predictions (on selected test points - here a grid from minimum to maximum temperature)
    min_x = tpptr_df['x'].min()
    max_x = tpptr_df['x'].max()
    
    # Prepare prediction storage
    predictions_full_mean_list = []
    predictions_full_conf_lower_list = []
    predictions_full_conf_upper_list = []
    fconf_full_lower_list = []
    fconf_full_upper_list = []
    covariance_matrices = {}
    
    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        # Specify test point locations
        test_x = torch.linspace(min_x, max_x, n_predict).to(dtype = dtype, device=device)
        
        # Iterate through each batched model
        for model, likelihood, meta_chunk in full_model_list:
            # Forward pass for all tasks in this batch
            f_pred = model(test_x)  # Batched output: [num_tasks, n_predict]
            
            # Get predictions with likelihood
            pred_dist = likelihood(f_pred)
            
            # Extract predictions for each task in the batch
            for task_idx, meta in enumerate(meta_chunk):
                protein = meta['protein']
                condition = meta['condition']
                
                # Extract this task's predictions
                if model.num_tasks > 1:
                    # Batched model
                    task_mean = f_pred.mean[task_idx].cpu()
                    task_conf_lower = f_pred.confidence_region()[0][task_idx].cpu()
                    task_conf_upper = f_pred.confidence_region()[1][task_idx].cpu()
                    
                    pred_mean = pred_dist.mean[task_idx].cpu()
                    pred_conf_lower = pred_dist.confidence_region()[0][task_idx].cpu()
                    pred_conf_upper = pred_dist.confidence_region()[1][task_idx].cpu()
                    
                    # Store covariance matrix
                    covar_matrix = f_pred.covariance_matrix[task_idx].cpu()
                else:
                    # Single task model
                    task_mean = f_pred.mean.cpu()
                    task_conf_lower = f_pred.confidence_region()[0].cpu()
                    task_conf_upper = f_pred.confidence_region()[1].cpu()
                    
                    pred_mean = pred_dist.mean.cpu()
                    pred_conf_lower = pred_dist.confidence_region()[0].cpu()
                    pred_conf_upper = pred_dist.confidence_region()[1].cpu()
                    
                    # Store covariance matrix
                    covar_matrix = f_pred.covariance_matrix.cpu()
                
                predictions_full_mean_list.append(pred_mean.numpy())
                predictions_full_conf_lower_list.append(pred_conf_lower.numpy())
                predictions_full_conf_upper_list.append(pred_conf_upper.numpy())
                fconf_full_lower_list.append(task_conf_lower.numpy())
                fconf_full_upper_list.append(task_conf_upper.numpy())
                
                # Store covariance matrix with unique key
                key = f"{protein}_{condition}"
                covariance_matrices[key] = covar_matrix

    # Create dataframe with fits (proteins)
    result_rows = []
    pred_idx = 0  # Track which prediction we're using

    # Build dataframe using actual metadata order
    for model, likelihood, meta_chunk in full_model_list:
        for task_idx, meta in enumerate(meta_chunk):
            protein = meta['protein']
            condition = meta['condition']
            
            # Add rows for each test point using the correct metadata
            for i in range(len(test_x)):
                result_rows.append({
                    'uniqueID': protein,
                    'condition': condition,
                    'x': test_x.cpu().numpy()[i],
                    'y': predictions_full_mean_list[pred_idx][i],
                    'conf_lower': fconf_full_lower_list[pred_idx][i],
                    'conf_upper': fconf_full_upper_list[pred_idx][i],
                    'conflik_lower': predictions_full_conf_lower_list[pred_idx][i],
                    'conflik_upper': predictions_full_conf_upper_list[pred_idx][i],
                    'type': 'fitted'
                })
            pred_idx += 1

    result_full_df = pd.DataFrame(result_rows)
    
    # Update the checklist
    if is_notebook():
        live_update_message = f"Full model predictions complete"
        update_checklist(tasks, live_update={"task": "3. Evaluate and predict models",
                                            "subtask": "Full model","message": live_update_message}) 
        # Update the checklist
        if null_dataset == False:
            tasks[2][2][0] = ("Full model", True)
            update_checklist(tasks)
    
    ############################## JOINT MODEL PREDICTIONS ##############################
    
    # Predictions for the joint model
    predictions_joint_mean = []
    predictions_joint_conf_lower = []
    predictions_joint_conf_upper = []
    predictions_joint_conflik_lower = []
    predictions_joint_conflik_upper = []
    covariance_matrices_joint = {}

    for i, (combined_model, combined_likelihood) in enumerate(zip(joint_model.models, joint_likelihood.likelihoods)):
        protein = proteins2test[i]
        condition = 'joint'
        
        # Set submodel into eval mode
        combined_model.eval()
        combined_likelihood.eval()
        
        with torch.no_grad(), gpytorch.settings.fast_pred_var():
            # GP posterior predictions
            f_pred_joint = combined_model(test_x)
            predictions_joint_conf_lower.append(f_pred_joint.confidence_region()[0].cpu().detach().numpy())
            predictions_joint_conf_upper.append(f_pred_joint.confidence_region()[1].cpu().detach().numpy())
            
            # Likelihood predictions (with observation noise)
            pred_joint_lik = combined_likelihood(f_pred_joint)
            predictions_joint_mean.append(pred_joint_lik.mean.cpu().detach().numpy())
            predictions_joint_conflik_lower.append(pred_joint_lik.confidence_region()[0].cpu().detach().numpy())
            predictions_joint_conflik_upper.append(pred_joint_lik.confidence_region()[1].cpu().detach().numpy())

            # Calculate the covariance matrix
            covariance_matrix = f_pred_joint.covariance_matrix.cpu().detach()
        
            # Create a unique key for the covariance matrix
            key = f"{protein}_{condition}"
            covariance_matrices_joint[key] = covariance_matrix

    # Dataframe with fits for joint model 
    result_joint_df = pd.DataFrame({
        'uniqueID': np.repeat(proteins2test, len(test_x)),
        'condition': np.tile(np.repeat('joint', len(test_x)), n_prot),
        'type': np.tile(np.repeat('fitted', len(test_x)), n_prot),
        'x': np.tile(test_x.cpu().numpy(), n_prot),
        'y': np.array(predictions_joint_mean).flatten(),
        'conf_lower': np.array(predictions_joint_conf_lower).flatten(),
        'conf_upper': np.array(predictions_joint_conf_upper).flatten(),
        'conflik_lower': np.array(predictions_joint_conflik_lower).flatten(),
        'conflik_upper': np.array(predictions_joint_conflik_upper).flatten()
    })

    # Update the checklist
    if is_notebook():
        live_update_message = f"Joint model predictions complete"
        update_checklist(tasks, live_update={"task": "3. Evaluate and predict models",
                                            "subtask": "Joint model","message": live_update_message}) 
        # Update the checklist
        if null_dataset == False:
            tasks[2][2][1] = ("Joint model", True)
            update_checklist(tasks)

    # Merge dataframes
    result_df = pd.concat([result_full_df, result_joint_df])

    ############################## EXTRACT TRAINING DATA ##############################
    
    # Make dataframe of training input from batched models
    train_data_list = []
    
    for model, likelihood, meta_chunk in full_model_list:
        train_x = model.train_inputs[0].cpu()  # Same x for all tasks in batch
        train_y = model.train_targets.cpu()  # [num_tasks, n_points]
        
        # Extract for each task
        for task_idx, meta in enumerate(meta_chunk):
            protein = meta['protein']
            condition = meta['condition']
            
            if model.num_tasks > 1:
                # Batched model - extract this task's targets
                y_values = train_y[task_idx].numpy()
            else:
                # Single task model
                y_values = train_y.numpy()
            
            x_values = train_x.numpy()
            
            for x_val, y_val in zip(x_values, y_values):
                train_data_list.append({
                    'uniqueID': protein,
                    'condition': condition,
                    'x': x_val,
                    'y': y_val,
                    'type': 'measured'
                })
    
    # Add input to predicton dataframe
    inputs_df = pd.DataFrame(train_data_list)
    # Make sure values are numeric
    inputs_df['x'] = inputs_df['x'].astype(float)
    inputs_df['y'] = inputs_df['y'].astype(float)
    
    prediction_result_df = pd.concat([inputs_df, result_df])

    # Add predictions to result dict
    prediction_result_dict = deepcopy(result_dict)
    prediction_result_dict.update({
        "gp_result_df": prediction_result_df,
        "gp_full_covariance": covariance_matrices,
        "gp_joint_covariance": covariance_matrices_joint
    })
    
    # Update the checklist
    if is_notebook():
        if null_dataset == False:
            tasks[2][2][2] = ("Combine and save predictions", True)
            update_checklist(tasks)
            tasks[2] = ("3. Evaluate and predict models", True, tasks[2][2])
            update_checklist(tasks)

        else:
            # Mark the generation of a joint model as complete
            tasks[3][2][2] = ("Evaluate and predict null model", True)
            update_checklist(tasks)
            tasks[3] = ("4. Build and fit null model", True, tasks[3][2])
            update_checklist(tasks)
    
    # Save updates to summery file        
    if null_dataset == False:        
        with open(f'{output_path}/{pert}_Thermal_Tracks_summary.txt', 'a') as f:
            f.write("\n" + "="*70+"\n")
            f.write("PREDICTION SUMMARY\n")
            f.write("="*70+"\n")
            f.write(f"  Full model predictions: {len(result_full_df)}\n")
            f.write(f"  Joint model predictions: {len(result_joint_df)}\n")
            f.write(f"  Training data points: {len(inputs_df)}\n")
            f.write(f"  Total output rows: {len(prediction_result_df)}\n")
            f.write("="*70 + "\n")
    
    return prediction_result_dict
  