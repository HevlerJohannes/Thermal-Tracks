#! /usr/bin/env python3
import pandas as pd
import numpy as np
import torch
import gpytorch
from copy import deepcopy

# Import util functions
from utils import update_checklist
from utils import compute_likelihood_ratio

# Import gp modules

def define_joint_model(result_dict, parameters, null_dataset = [True, False]):

    # Update the checklist
    if null_dataset == False:
    
        tasks = [
        ("1. Build and fit full model", True),
        ("2. Create a joint model and null dataset", False, [
            (f"Combine training input and targets for {parameters['control_condition']} & {parameters['perturbation']}", False),
            ("Compute mll for joint model", False),
            ("Sample from the joint model (sampling of true negatives, null dataset)", False),
            ("Save joint model, and null dataset", False)
        ]),
        ("3. Evaluate and predict models", False),
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
            ("Create a joint model", False),
            ("Evaluate and predict null model", False)
        ]),
        ("5. Compute likelihood ratio test statistics", False),
        ("6. Combine and create result files", False)
        ]
            
    # Display the updated checklist
    update_checklist(tasks)    
    
    # Load data from training
    full_model = deepcopy(result_dict['full_model_list'])
    metadata_list = deepcopy(result_dict['full_model_order'])
    mll_values_full_model_df = deepcopy(result_dict['full_mll_values'])
    conds = deepcopy(result_dict['conditions'])
    list_state_dict = deepcopy(result_dict["full_state_dict_list"])
   
    # Load the parameters of the trained model
    for submodels, sub_state_dict in zip(full_model.models, list_state_dict):
        submodels.load_state_dict(sub_state_dict)

    ############################## Joint model ##############################
    # Create a joint model by combining training input and targets of the two conditions for each protein by joining model of the control and drug condition for respective proteins

    joint_models = []
    joint_likelihoods = []
    combined_inputs = []
    combined_targets = []
    mll_values = []
    uniqueID_list = []
    join_control_df = pd.DataFrame(columns=['protein_model_join_a', 'protein_model_join_b', 'condition_model_join_a', 'condition_model_join_b'])
    all_samples = []
    
    # Loop through model list and join models (control and drug) for each protein 
    for i in range(0, len(full_model.models), 2):
        
        # Get protein id for model merge
        protein_id = metadata_list[i][0]

        # Create dataframe to verify that models were correctly joint
        join_control_df.loc[i, ['protein_model_join_a']] = metadata_list[i][0]
        join_control_df.loc[i, ['protein_model_join_b']] = metadata_list[i+1][0]
        join_control_df.loc[i, ['condition_model_join_a']] = metadata_list[i][1]
        join_control_df.loc[i, ['condition_model_join_b']] = metadata_list[i+1][1]
        
        assert metadata_list[i][0] == metadata_list[i+1][0], 'Failing to join the correct protein models!'
        
        # Retrieve the models, load dict
        model_i = full_model.models[i]
        state_dict_model_i = deepcopy(model_i.state_dict())
        model_i_plus1 = full_model.models[i+1]
        state_dict_model_i_plus1 = deepcopy(model_i_plus1.state_dict())
            
            
        # Combine training data from both models
        combined_train_input = torch.cat([model_i.train_inputs[0], model_i_plus1.train_inputs[0]], dim=0)
        combined_train_target = torch.cat([model_i.train_targets, model_i_plus1.train_targets], dim=0)

        # Verify that the combined data sizes are as expected
        assert combined_train_input.size(0) == model_i.train_inputs[0].size(0) + model_i_plus1.train_inputs[0].size(0), 'Size of training input sizes is incorrect - Joining models stopped!'
        assert combined_train_target.size(0) == model_i.train_targets.size(0) + model_i_plus1.train_targets.size(0), 'Size of training targets are incorrect - Joining models stopped!'

        # Create a new model and likelihood for the joint model
        combined_likelihood = gpytorch.likelihoods.GaussianLikelihood()
        combined_model = ExactGPModel(combined_train_input, combined_train_target, combined_likelihood, lengthscale_prior=None, lengthscale_minconstraint=None)

        # Combine state dictionaries and average all parameters       
        if parameters['alpha_joint'] != None:
            alpha = parameters['alpha_joint']
            combined_state_dict = {}
            for key in state_dict_model_i.keys():
                combined_state_dict[key] = (state_dict_model_i[key] + (1-alpha) * state_dict_model_i_plus1[key])

        else:
            combined_state_dict = {}
            for key in state_dict_model_i.keys():
                combined_state_dict[key] = (state_dict_model_i[key] + state_dict_model_i_plus1[key]) / 2

        # Load dict for joint model
        combined_model.load_state_dict(combined_state_dict)

        joint_models.append(combined_model)
        joint_likelihoods.append(combined_likelihood)
        combined_inputs.append(combined_train_input)
        combined_targets.append(combined_train_target)
        uniqueID_list.append(protein_id)

        # Update progress of joining
        if null_dataset == False:

            live_update_message = f"(Models joint: {i + 1}/{(len(full_model.models))})"
            update_checklist(tasks, live_update={"task": "2. Create a joint model and null dataset",
                                                "subtask": f"Combine training input and targets for {parameters['control_condition']} & {parameters['perturbation']}",
                                                "message": live_update_message})
        else:
            live_update_message = f"(Models joint: {i + 1}/{(len(full_model.models))})"
            update_checklist(tasks, live_update={"task": "4. Build and fit null model",
                                                "subtask": "Create a joint model",
                                                "message": live_update_message})
            
        # Compute mll for joint model
        mll = gpytorch.mlls.ExactMarginalLogLikelihood(combined_likelihood, combined_model)
        output = combined_model(*combined_model.train_inputs)
        mll_value = mll(output, combined_model.train_targets).item()
        mll_values.append(mll_value)

        # Sample from the joint model (sampling of true negatives) - Samples are collected at the same temperatures that were used as inputs for training the model
        if null_dataset == False:
        
            # Initialize a dataframe to store all samples (dataframe with protein information and original training inputs)
            prior_joint_df = pd.DataFrame()

            # Get original training inputs for a protein across drug and control condition 
            test_x_joint_model = deepcopy(combined_model.train_inputs)
            test_x_joint_model = tuple(x.to(torch.float64) for x in test_x_joint_model)

            # Ensure the model parameters are of type torch.float64
            for param in combined_model.parameters():
                param.data = param.data.to(torch.float64)

            for param in combined_likelihood.parameters():
                param.data = param.data.to(torch.float64)

            with gpytorch.settings.prior_mode(True):
                combined_model.eval()
                combined_model.eval()
                prior_joint_model = combined_likelihood(combined_model(*test_x_joint_model))
                prior_joint = prior_joint_model.sample().detach().numpy()
                
                
            # Get information about the sampled condition
            # Initialize the previous temperature to compare
            temp = combined_model.train_inputs
            temp_flat = temp[0]
            output_conds = []
            current_cond_index = 0
            previous_temp = temp_flat[0]

            # Iterate over the temp list starting from the second element
            for i in range(len(temp_flat)):
                t = temp_flat[i]
                # Check if the current temperature is lower than the previous temperature
                if i > 0 and t < previous_temp:
                    # Switch to the next condition in the conds list
                    current_cond_index = (current_cond_index + 1) % len(conds)
                
                # Append the current condition to the output list
                output_conds.append(conds[current_cond_index])

                # Update the previous temperature
                previous_temp = t

            # Convert the list to a NumPy array
            output_conds_array = np.array(output_conds)

            # make final dataframe
            prior_joint_df.loc[:, ('uniqueID')] = np.repeat(protein_id, len(temp_flat))
            prior_joint_df.loc[:, ('x')] = temp_flat
            prior_joint_df.loc[:, ('y')] = prior_joint
            prior_joint_df.loc[:, ('condition')] = output_conds_array
            prior_joint_df.loc[:, ('model')] = 'joint'
            prior_joint_df.loc[:, ('comparison')] = f'{conds[0]} - {conds[1]}'
            prior_joint_df.loc[:, ('data')] = 'joint model priors'
                
            all_samples.append(prior_joint_df)

            # create dataframe of sampled data                
            all_samples_df = pd.concat(all_samples)
        
        else:
            None
        
    # Convert the joint model and likelihood lists to IndependentModelList and LikelihoodList
    joint_model_list = gpytorch.models.IndependentModelList(*joint_models)
    joint_likelihood_list = gpytorch.likelihoods.LikelihoodList(*joint_likelihoods)
    
    # Update the checklist
    if null_dataset == False:

        tasks[1][2][0] = (f"Combine training input and targets for {parameters['control_condition']} & {parameters['perturbation']}", True)
        update_checklist(tasks)
        tasks[1][2][1] = (f"Compute mll for joint model", True)
        update_checklist(tasks)
        tasks[1][2][2] = (f"Sample from the joint model (sampling of true negatives, null dataset)", True)
        update_checklist(tasks)

    else:
        None 

    # Save mll of joined model in dataframe
    mll_values_joint_model_df = pd.DataFrame({'protein': uniqueID_list,
    'condition': 'joint', 'mll': mll_values})

    # Save parameters for joint model
    list_state_dict_joint = []
    for submodel in joint_model_list.models:
            list_state_dict_joint.append(submodel.state_dict())    

    # Compute likelihood ratios for full vs joint model
    lr_values_full_vs_joint = compute_likelihood_ratio(mll_full_df = mll_values_full_model_df, mll_joint_df = mll_values_joint_model_df)

    # Create a dict with model dict and training information
    join_model_result_dict = deepcopy(result_dict)
    if null_dataset == False:
        
        join_model_result_dict.update({
                "joint_model_list" : joint_model_list,
                "joint_likelihood_list" : joint_likelihood_list,
                "joint_state_dict_list" : list_state_dict_joint,
                "joint_mll_values" : mll_values_joint_model_df,
                "lr_values_full_vs_joint": lr_values_full_vs_joint,
                "model_join_control_df" : join_control_df,
                "sampled_prior_df": all_samples_df
                })
        # Mark the generation of a joint model and null dataset as complete
        tasks[1][2][3] = (f"Save joint model, and null dataset", True)
        update_checklist(tasks)
        tasks[1] = ("2. Create a joint model and null dataset", True, tasks[1][2])
        update_checklist(tasks)
    
    else:
        
        join_model_result_dict.update({
                "joint_model_list" : joint_model_list,
                "joint_likelihood_list" : joint_likelihood_list,
                "joint_state_dict_list" : list_state_dict_joint,
                "joint_mll_values" : mll_values_joint_model_df,
                "lr_values_full_vs_joint": lr_values_full_vs_joint,
                "model_join_control_df" : join_control_df,
                })
        
        # Mark the generation of a joint model as complete
        tasks[3][2][1] = ("Create a joint model", True)
        update_checklist(tasks)    
    
    return join_model_result_dict

# ExactGPModel
class ExactGPModel(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood, mean = gpytorch.means.ZeroMean(), lengthscale_prior=None, lengthscale_minconstraint = ['min', 'mean', 'median', 'max', None]):
        super(ExactGPModel, self).__init__(train_x, train_y, likelihood)
        #self.mean_module = gpytorch.means.ConstantMean() # posterior mean need not be zero, prior has high influence only for low noise values - e.g. actual noise as fixed noise
        self.mean_module = mean
        #self.covar_module = gpytorch.kernels.ScaleKernel(gpytorch.kernels.RBFKernel()) 
        
        if not lengthscale_minconstraint is None:
            train_x_values = train_x.unique()
            DistVec = train_x_values[1:len(train_x_values)]-train_x_values[0:len(train_x_values)-1] # compute the distance between each pair of consecutive temperatures
            
            if lengthscale_minconstraint == 'min':
                Constt = min(DistVec) # we use the minium distance between temperatures as a lower limit for the lengthscale
            elif lengthscale_minconstraint == 'mean':
                Constt = torch.mean(DistVec) # we use the mean distance between temperatures as a lower limit for the lengthscale
            elif lengthscale_minconstraint == 'median':
                Constt = torch.median(DistVec) # we use the median distance between temperatures as a lower limit for the lengthscale
            elif lengthscale_minconstraint == 'max':
                Constt = max(DistVec) # we use the max distance between temperature as a lower limit for the lengthscale
            
            lengthscale_constraint = gpytorch.constraints.GreaterThan(Constt)
        else:
            lengthscale_constraint = None # there is no lower limit for the lengthscale
    
    
        self.covar_module = gpytorch.kernels.RBFKernel(lengthscale_prior = lengthscale_prior, lengthscale_constraint = lengthscale_constraint)
            

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)