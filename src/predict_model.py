#! /usr/bin/env python3
import numpy as np
import pandas as pd
import torch
import gpytorch
import pandas as pd
from copy import deepcopy


# Import util functions
from utils import update_checklist


def predict_and_evaluate(result_dict, parameters, null_dataset = [True, False]):

    # Update the checklist
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

    # get model and lik lists for full and joint models
    full_model = result_dict['full_model_list']
    full_likelihood = result_dict['full_likelihood_list']
    joint_model = result_dict['joint_model_list']
    joint_likelihood = result_dict['joint_likelihood_list']

    # get additional model information
    n_models = result_dict['n_full_models']
    proteins2test = result_dict['proteins']
    n_prot = result_dict['n_proteins']
    conds = result_dict['conditions']
    n_cond = result_dict['n_conditions']
    tpptr_df = result_dict['exactgp_input'].copy()
    
    # define number of prediction points
    n_predict = parameters['n_predictions']
    
    ############################## MODEL EVALUATE and PREDICT ##############################
    # Predictions for the full model
    # Set model into eval mode
    full_model.eval()
    full_likelihood.eval()

    # Make predictions (on selected test points - here a grid from minimum to maximum temperature)
    min_x = tpptr_df['x'].min()
    max_x = tpptr_df['x'].max()

    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        # specify test point locations
        test_x = torch.linspace(min_x, max_x, n_predict).double()
        test_x_list = [test_x] * n_models # use the same x for all proteins and conditions
        # Predictions for all outcomes as a list
        flist_full = full_model(*test_x_list) # posterior process in each model
        # Compute the covariance matrix
        covar_full_list = [f.covariance_matrix for f in flist_full]
        #fmean_full = torch.stack([f.mean for f in flist_full]) # should be the same as predictions_mean
        fconf_full_lower = torch.stack([f.confidence_region()[0] for f in flist_full])
        fconf_full_upper = torch.stack([f.confidence_region()[1] for f in flist_full])
        predictions_full = full_likelihood(*full_model(*test_x_list))
        predictions_full_mean = torch.stack([pred.mean for pred in predictions_full]) # array with fitted mean per protein and condition
        predictions_full_conf_lower = torch.stack([pred.confidence_region()[0] for pred in predictions_full]) # array with lower confidence region per peptide and condition
        predictions_full_conf_upper = torch.stack([pred.confidence_region()[1] for pred in predictions_full]) # array with upper confidence region per peptide and condition
    
    # data frame with fits (proteins)
    result_full_df = pd.DataFrame({'uniqueID' : np.repeat(proteins2test, n_cond*len(test_x)),
                        'condition' : np.tile(np.repeat(conds, len(test_x)), n_prot),
                        'x' : np.tile(test_x, n_prot * n_cond)})
    result_full_df['y'] = torch.flatten(predictions_full_mean)
    result_full_df['conf_lower'] = torch.flatten(fconf_full_lower)
    result_full_df['conf_upper'] = torch.flatten(fconf_full_upper)
    result_full_df['conflik_lower'] = torch.flatten(predictions_full_conf_lower)
    result_full_df['conflik_upper'] = torch.flatten(predictions_full_conf_upper)
    result_full_df['type'] = "fitted"

    # Store the covariance matrices
    covariance_matrices = {}
    for i, protein in enumerate(proteins2test):
        for j, condition in enumerate(conds):
            key = f"{protein}_{condition}"
            covariance_matrices[key] = covar_full_list[i * len(conds) + j]

    
    # Update the checklist
    if null_dataset == False:
        tasks[2][2][0] = (f"Full model", True)
        update_checklist(tasks)
    else:
        None
    
    # Predictions for the joint model
    predictions_joint_mean = []
    predictions_joint_conf_lower = []
    predictions_joint_conf_upper = []
    covariance_matrices_joint = {}

    for i, (combined_model, combined_likelihood) in enumerate(zip(joint_model.models, joint_likelihood.likelihoods)):
        protein = proteins2test[i]
        condition = 'joint'
        # Set submodel into eval mode
        combined_model.eval()
        combined_likelihood.eval()
        with torch.no_grad(), gpytorch.settings.fast_pred_var():
            f_pred_joint = combined_model(test_x)  # Add unsqueeze to add a dimension to test_x
            predictions_joint_mean.append(f_pred_joint.mean.detach().numpy())
            predictions_joint_conf_lower.append(f_pred_joint.confidence_region()[0].detach().numpy())
            predictions_joint_conf_upper.append(f_pred_joint.confidence_region()[1].detach().numpy())

            # Calculate the covariance matrix
            covariance_matrix = f_pred_joint.covariance_matrix.detach()
        
            # Create a unique key for the covariance matrix
            key = f"{protein}_{condition}"
            covariance_matrices_joint[key] = covariance_matrix


    # data frame with fits for joint model 
    result_joint_df = pd.DataFrame({'uniqueID' : np.repeat(proteins2test, len(test_x)),
                        'condition' : np.tile(np.repeat('joint', len(test_x)), n_prot),
                        'type' : np.tile(np.repeat('fitted', len(test_x)), n_prot),
                        'x' : np.tile(test_x, n_prot),
                        'y' : np.array(predictions_joint_mean).flatten(),
                        'conf_lower' : np.array(predictions_joint_conf_lower).flatten(),
                        'conf_upper' : np.array(predictions_joint_conf_upper).flatten()})

    # Update the checklist
    if null_dataset == False:
        tasks[2][2][1] = (f"Joint model", True)
        update_checklist(tasks)
    else:
        None

    # merge dataframes
    result_df = pd.concat([result_full_df, result_joint_df])

    # Make data frame of training input
    x_train = [(torch.stack(submodel.train_inputs).flatten()) for submodel in full_model.models]
    n_xtrain_cond = [len(x) for x in x_train]
    y_train = [submodel.train_targets for submodel in full_model.models]
    [len(y) for y in y_train]
    inputs_df =  pd.DataFrame({'uniqueID' : np.repeat(np.repeat(proteins2test, n_cond), n_xtrain_cond),
                            'condition' : np.repeat(np.tile(conds, n_prot), n_xtrain_cond),
                            'x' : torch.cat(x_train),
                            'y': torch.cat(y_train)})
    inputs_df['type'] = "measured"
    prediction_result_df = pd.concat([inputs_df, result_df])

    # Add predictions to result dict
    prediction_result_dict = deepcopy(result_dict)
    prediction_result_dict.update({
            "gp_result_df" : prediction_result_df,
            "gp_full_covariance" : covariance_matrices,
            "gp_joint_covariance" : covariance_matrices_joint})
    
    # Update the checklist
    if null_dataset == False:
        tasks[2][2][2] = (f"Combine and save predictions", True)
        update_checklist(tasks)
        tasks[2] = ("3. Evaluate and predict models", True, tasks[2][2])
        update_checklist(tasks)
    else:
        # Mark the generation of a joint model as complete
        tasks[3][2][2] = ("Evaluate and predict null model", True)
        update_checklist(tasks)
        tasks[3] = ("4. Build and fit null model", True, tasks[3][2])
        update_checklist(tasks)           
    
    return prediction_result_dict
            

  