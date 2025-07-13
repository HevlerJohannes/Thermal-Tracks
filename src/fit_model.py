#! /usr/bin/env python3
import pandas as pd
import os
import numpy as np
import csv
from matplotlib import pyplot as plt
import torch
import gpytorch
from gpytorch.mlls import SumMarginalLogLikelihood

# Import util functions
from utils import update_checklist

# Import gp modules
from build_model import building_exactgp_model

####################################################################
#                    Fit Multitask Exact GP                        #
####################################################################

def train_model(tpptr_df_input, parameters, null_dataset = [True, False]):
    
    if null_dataset == False:
         
        tasks = [
        ("1. Build and fit full model", False, [
            ("Prepare input for GP process", False),
            ("Build models", False),
            ("Initialize models", False),
            ("Train models", False),
            ("Compute mll for full models", False),
            ("Save models and training results", False)
        ]),
        ("2. Creating a joint model and null dataset", False),
        ("3. Evaluate and predict models", False),
        ("4. Build and fit null model", False),
        ("5. Compute likelihood ratio test statistics", False),
        ("6. Combine and create result files", False)]

    else:
         
         tasks = [
        ("1. Build and fit full model", True),
        ("2. Create a joint model and null dataset", True),
        ("3. Evaluate and predict full and joint models", True),
        ("4. Build and fit null model", False, [
            ("Train null model", False),
            ("Create a joint model", False),
            ("Evaluate and predict null model", False)
        ]),
        ("5. Compute likelihood ratio test statistics", False),
        ("6. Combine and create result files", False)
        ]
         
    # Update the checklist
    update_checklist(tasks)
    
    ############################## MODEL BUILDING ##############################
    # Prepare input for gp process
    tpptr_df = tpptr_df_input.copy()
    pert = parameters['perturbation']
    control = parameters['control_condition']
    output_path = str(parameters['result_dir'])
    proteins2test = tpptr_df['uniqueID'].unique()
    conds = tpptr_df['condition'].unique()
    
    # Update the checklist
    if null_dataset == False:
         
        tasks[0][2][0] = ("Prepare input for GP process", True)
        update_checklist(tasks)
    
    else:
         None
    
    # Parameters for training
    lengthscale_prior = parameters['lengthscale_prior']
    lengthscale_minconstraint = parameters['lengthscale_minconstraint']
    lengthscale_mult = parameters['lengthscale_mult']
    n_iterations = parameters['training_iterations']
    
    ############################## MODEL BUILDING - FULL MODEL ##############################
    model_list, l_list, metadata_list, n_cond, n_prot, n_models = building_exactgp_model(tpptr_df, proteins2test, conds, lengthscale_prior, lengthscale_minconstraint, lengthscale_mult, mean = gpytorch.means.ZeroMean())
    
    # Update the checklist
    if null_dataset == False:
         
        tasks[0][2][1] = ("Build models", True)
        update_checklist(tasks)
    
    else:
         None
    
    full_model = gpytorch.models.IndependentModelList(*model_list)
    full_likelihood = gpytorch.likelihoods.LikelihoodList(*l_list)
    
    ############################## MODEL TRAINING ##############################
    # Fit the full model

    full_model.train()
    full_likelihood.train()

    # Update the checklist
    if null_dataset == False:
         
        tasks[0][2][2] = ("Initialize model", True)
        update_checklist(tasks)

        tasks[0][2][3] = ("Train models", False)
        update_checklist(tasks)

    else:
         None

    # Optimize the models using the Adam optimizer
    optimizer = torch.optim.Adam([
            {'params': full_model.parameters()},  
        ], lr=parameters['learningRate'], amsgrad = parameters['amsgrad'])

    sum_mll = SumMarginalLogLikelihood(full_likelihood, full_model)

    # Train
    training_iterations = n_iterations
    LossValues = []
    training_info = []
    for i in range(training_iterations):
        optimizer.zero_grad()
        output = full_model(*full_model.train_inputs)
        loss = -sum_mll(output, full_model.train_targets)
        loss.backward()
        #print('Training iteration %d/%d - Loss: %.6f' % (i + 1, training_iterations, loss.item()), end='\r')
        LossValues.append(loss.item())
        optimizer.step()

        # Live update for each iteration
        if null_dataset == False:
        
            live_update_message = f"(Training iteration {i + 1}/{n_iterations} - Loss: {loss.item():.6f})"
            plot_data = {
                "task": "1. Build and fit full model",
                "subtask": "Train models",
                "x": list(range(1, len(LossValues) + 1)),
                "y": LossValues,
                "label": "Loss",
                "xlabel": "Iteration",
                "ylabel": "Loss",
                "title": "Training Loss Over Iterations"
            }
            update_checklist(tasks, live_update={"task": "1. Build and fit full model", "subtask": "Train models", "message": live_update_message}, plot_data=plot_data)

        else:
            live_update_message = f"(Training iteration {i + 1}/{n_iterations} - Loss: {loss.item():.6f})"
            update_checklist(
                tasks,
                live_update={
                    "task": "4. Build and fit null model",
                    "subtask": "Train null model",
                    "message": live_update_message
                }
            )
            

        if i > 1 and abs(LossValues[i] - LossValues[i-1]) <= 0.0000011:
                break

        # Get loss, lengthscale, and noise for each protein
        for idx, m in enumerate(full_model.models):
            lengthscale = m.covar_module.lengthscale.item()
            noise = full_likelihood.likelihoods[idx].noise_covar.noise.mean().item()  # Use mean if noise is a tensor
            training_info.append([i + 1, loss.item(), idx + 1, lengthscale, noise])

    # Update the checklist
    if null_dataset == False:         
        tasks[0][2][3] = ("Train model", True)  # Mark training as complete
        update_checklist(tasks)
    else:
         None        

    # Save loss, lengthscale, and noise for each protein in dataframe
    training_info_df = pd.DataFrame(training_info, columns=['iteration', 'loss', 'model', 'lengthscale', 'noise'])
    n_rows = len(training_info_df)
    training_info_df['protein'] = np.tile(proteins2test, n_rows // len(proteins2test) + 1)[:n_rows] 
    training_info_df['condition'] = np.tile(np.repeat(conds, n_prot), n_rows // (n_prot * n_cond) + 1)[:n_rows]
    training_info_df = training_info_df[['protein', 'condition', 'iteration', 'loss', 'model', 'lengthscale', 'noise']]   
    
    # generate convergence plot
    if null_dataset == False:
        plt.plot(range(len(LossValues)), LossValues)
        plt.xlabel('Number of iterations', fontsize=12)
        plt.ylabel('Loss', fontsize=12)
        plt.title('Full Model training - '+str(pert), fontsize=14)
        plt.savefig(output_path + '/loss_full_model_gp_'+str(control)+'_'+str(pert)+'.pdf')
        plt.close()
    else:
        plt.plot(range(len(LossValues)), LossValues)
        plt.xlabel('Number of iterations', fontsize=12)
        plt.ylabel('Loss', fontsize=12)
        plt.title('Null Model training - '+str(pert), fontsize=14)
        plt.savefig(output_path + '/loss_null_model_gp_'+str(control)+'_'+str(pert)+'.pdf')
        plt.close()
         

    # save Loss
    if null_dataset == False:
         # Save the DataFrame to a CSV file
         output_file_path = f'{output_path}/{pert}_training_info_df.csv'
         training_info_df.to_csv(output_file_path, index=False)
         with open(f'{output_path}/{pert}_loss_full_model.csv', mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Iteration', 'Loss'])
            for idx, loss_value in enumerate(LossValues):
                writer.writerow([idx + 1, loss_value])
    else:
         with open(f'{output_path}/{pert}_loss_null_model.csv', mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Iteration', 'Loss'])
            for idx, loss_value in enumerate(LossValues):
                writer.writerow([idx + 1, loss_value])
         
                   

    # Compute the mll for full model
    mll_values_full_model = [sub_mll(output, target).item() for sub_mll, output, target in zip(sum_mll.mlls, output, full_model.train_targets)]
    mll_values_full_model_df = pd.DataFrame({'protein' : np.repeat(proteins2test, n_cond), 'condition' : np.tile(conds, n_prot), 'mll' : mll_values_full_model})

    # Update the checklist
    if null_dataset == False:
         
        tasks[0][2][4] = ("Compute mll for full model", True)
        update_checklist(tasks)
    
    # save parameters for all models after training
    list_state_dict_full = []
    for submodel in full_model.models:
            list_state_dict_full.append(submodel.state_dict())

    # Create a dict with model dict and training information   
    train_results_dict = {
            "full_model_list" : full_model,
            "full_likelihood_list" : full_likelihood,
            "full_state_dict_list" : list_state_dict_full,
            "full_fit_parameters_df" : training_info_df,
            "full_mll_values" : mll_values_full_model_df,
            "modeled_proteins_list" : proteins2test,
            "full_model_order" : metadata_list,
            "n_full_models" : n_models,
            "proteins" : proteins2test,
            "n_proteins": n_prot,
            "conditions" : conds,
            "n_conditions": n_cond,
            "exactgp_input" : tpptr_df,
            "output_path" : output_path
            }
    if null_dataset == False:

        # Update the checklist
        tasks[0][2][5] = ("Save models and training results", True)
        update_checklist(tasks)

        # Mark the training process as complete
        tasks[0] = ("1. Build and fit full model", True, tasks[0][2])
        update_checklist(tasks)
    
    else:       
        # Mark the generation of a joint model as complete
        tasks[3][2][0] = ("Train null model", True)
        update_checklist(tasks)  
         
        
    return train_results_dict