#! /usr/bin/env python3
import pandas as pd
import numpy as np
import torch
import gpytorch
from copy import deepcopy

# Import util functions
from utils import update_checklist, is_notebook, compute_likelihood_ratio

class ExactGPModel(gpytorch.models.ExactGP):
    """
    Single-task ExactGP model for joint models (non-batched)
    """
    def __init__(self, train_x, train_y, likelihood, scale_kernel=False, 
                 mean=gpytorch.means.ZeroMean(), 
                 lengthscale_prior=None, 
                 lengthscale_minconstraint=None, 
                 lengthscale_mult=1.0):
        super(ExactGPModel, self).__init__(train_x, train_y, likelihood)

        self.mean_module = mean
        
        # Calculate lengthscale constraint
        lengthscale_constraint = None
        if lengthscale_minconstraint is not None:
            train_x_values = train_x.unique()
            if len(train_x_values) > 1:
                DistVec = train_x_values[1:] - train_x_values[:-1]
                
                if lengthscale_minconstraint == 'min':
                    Constt = lengthscale_mult * torch.min(DistVec)
                elif lengthscale_minconstraint == 'mean':
                    Constt = lengthscale_mult * torch.mean(DistVec)
                elif lengthscale_minconstraint == 'median':
                    Constt = lengthscale_mult * torch.median(DistVec)
                elif lengthscale_minconstraint == 'max':
                    Constt = lengthscale_mult * torch.max(DistVec)
                
                lengthscale_constraint = gpytorch.constraints.GreaterThan(Constt)
        
        # Covariance module
        if scale_kernel:
            self.covar_module = gpytorch.kernels.ScaleKernel(
                gpytorch.kernels.RBFKernel(
                    lengthscale_prior=lengthscale_prior,
                    lengthscale_constraint=lengthscale_constraint
                )
            )
        else:
            self.covar_module = gpytorch.kernels.RBFKernel(
                lengthscale_prior=lengthscale_prior,
                lengthscale_constraint=lengthscale_constraint
            )

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


def define_joint_model_batched(result_dict, parameters, null_dataset=[True, False]):

    # Update the checklist
    if is_notebook():
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

    # Load data from training - UPDATED for batched model format
    full_model_list = deepcopy(result_dict['full_model_list'])  # List of (model, likelihood, metadata) tuples
    metadata_list = deepcopy(result_dict['full_model_metadata'])  # Flat list of metadata dicts
    mll_values_full_model_df = deepcopy(result_dict['full_mll_values'])
    conds = deepcopy(result_dict['conditions'])
    list_state_dict = deepcopy(result_dict["full_state_dict_list"])
    scale_kernel = parameters.get('ScaledKernel', False)
    output_path = str(parameters['result_dir'])
    pert = parameters['perturbation']
    
    # Detect device from the trained models
    device_str = result_dict.get('device', 'cpu')
    device = torch.device(device_str)
    dtype = parameters['dtype']
    # Load the parameters of the trained model
    for (model, likelihood, _), state_dict in zip(full_model_list, list_state_dict):
        model.load_state_dict(state_dict, strict=False)

    ############################## Joint model ##############################
    # Create a joint model by combining training input and targets of the two conditions for each protein
    
    # First, organize metadata by protein and condition
    protein_task_map = {}
    
    for meta in metadata_list:
        protein = meta['protein']
        condition = meta['condition']
        
        if protein not in protein_task_map:
            protein_task_map[protein] = {}
        
        protein_task_map[protein][condition] = (
            meta['chunk_index'], 
            meta['task_index'], 
            meta
        )
    
    # Verify all proteins have both conditions
    proteins_to_process = []
    for protein, cond_dict in protein_task_map.items():
        if len(cond_dict) == len(conds):
            proteins_to_process.append(protein)
        else:
            print(f"Warning: Protein {protein} missing conditions, skipping...")
    
    joint_models = []
    joint_likelihoods = []
    combined_inputs = []
    combined_targets = []
    mll_values = []
    uniqueID_list = []
    join_control_df = pd.DataFrame(columns=['protein_model_join_a', 'protein_model_join_b', 
                                            'condition_model_join_a', 'condition_model_join_b'])
    all_samples = []

    # Loop through proteins and create joint models
    for idx, protein in enumerate(proteins_to_process):
        
        # Get metadata for both conditions
        cond_data = [protein_task_map[protein][cond] for cond in conds]
        
        # Verify we have both conditions
        assert len(cond_data) == 2, f"Protein {protein} should have exactly 2 conditions"
        
        # Create dataframe to verify that models were correctly joined
        join_control_df.loc[idx, 'protein_model_join_a'] = cond_data[0][2]['protein']
        join_control_df.loc[idx, 'protein_model_join_b'] = cond_data[1][2]['protein']
        join_control_df.loc[idx, 'condition_model_join_a'] = cond_data[0][2]['condition']
        join_control_df.loc[idx, 'condition_model_join_b'] = cond_data[1][2]['condition']
        
        # Extract training data for both conditions
        train_inputs = []
        train_targets = []
        state_dicts = []
        
        for chunk_idx, task_idx, meta in cond_data:
            # Get the batched model for this chunk
            model, likelihood, _ = full_model_list[chunk_idx]
            
            # Extract training data for this specific task
            train_x = model.train_inputs[0] 
            train_y = model.train_targets[task_idx]
            
            train_inputs.append(train_x)
            train_targets.append(train_y)
            
            # Extract state dict for this task from the batched model
            task_state_dict = {}
            full_state_dict = model.state_dict()
            
            for key, value in full_state_dict.items():
                if isinstance(value, torch.Tensor) and value.dim() > 0 and value.shape[0] == model.num_tasks:
                    task_state_dict[key] = value[task_idx]
                else:
                    task_state_dict[key] = value
            
            state_dicts.append(task_state_dict)

        # Track and store the split point for sampling of Null dataset - making sure that data gets correctly sampled 
        n_points_first_condition = train_inputs[0].size(0)
        
        # Combine training data from both conditions (stays on device)
        combined_train_input = torch.cat(train_inputs, dim=0)
        combined_train_target = torch.cat(train_targets, dim=0)
        
        # Verify that the combined data sizes are as expected
        assert combined_train_input.size(0) == train_inputs[0].size(0) + train_inputs[1].size(0), \
            'Size of training input sizes is incorrect - Joining models stopped!'
        assert combined_train_target.size(0) == train_targets[0].size(0) + train_targets[1].size(0), \
            'Size of training targets are incorrect - Joining models stopped!'
        
        # Create a new model and likelihood for the joint model
        combined_likelihood = gpytorch.likelihoods.GaussianLikelihood().to(device)
        combined_model = ExactGPModel(
            combined_train_input, 
            combined_train_target, 
            combined_likelihood,
            scale_kernel=scale_kernel,
            lengthscale_prior=parameters.get('lengthscale_prior'), 
            lengthscale_minconstraint=parameters.get('lengthscale_minconstraint'), 
            mean=parameters.get('MeanFunction'), 
            lengthscale_mult=parameters.get('lengthscale_mult', 1.0)
        ).to(device)
        
        # Combine state dictionaries and average all parameters
        combined_state_dict = {}
        for key in state_dicts[0].keys():
            # Average parameters from both conditions
            if isinstance(state_dicts[0][key], torch.Tensor):
                combined_state_dict[key] = (state_dicts[0][key] + state_dicts[1][key]) / 2
            else:
                combined_state_dict[key] = state_dicts[0][key]
        
        # Load dict for joint model
        combined_model.load_state_dict(combined_state_dict, strict=False)
        
        joint_models.append(combined_model)
        joint_likelihoods.append(combined_likelihood)
        combined_inputs.append(combined_train_input)
        combined_targets.append(combined_train_target)
        uniqueID_list.append(protein)
        
        # Update progress of joining
        if is_notebook():
            if null_dataset == False:
                live_update_message = f"(Models joined: {idx + 1}/{len(proteins_to_process)})"
                update_checklist(tasks, live_update={
                    "task": "2. Create a joint model and null dataset",
                    "subtask": f"Combine training input and targets for {parameters['control_condition']} & {parameters['perturbation']}",
                    "message": live_update_message
                })
            else:
                live_update_message = f"(Models joined: {idx + 1}/{len(proteins_to_process)})"
                update_checklist(tasks, live_update={
                    "task": "4. Build and fit null model",
                    "subtask": "Create a joint model",
                    "message": live_update_message
                })
        
        # Compute mll for joint model (all on same device now)
        mll = gpytorch.mlls.ExactMarginalLogLikelihood(combined_likelihood, combined_model)
        output = combined_model(*combined_model.train_inputs)
        mll_value = mll(output, combined_model.train_targets).item()
        mll_values.append(mll_value)
        
        # Sample from the joint model (sampling of true negatives)
        if null_dataset == False:
            # Get original training inputs for a protein across perturbation and control condition 
            test_x_joint_model = deepcopy(combined_model.train_inputs)
            test_x_joint_model = tuple(x.to(dtype) for x in test_x_joint_model)
            
            # Ensure the model parameters are of preset dtype
            for param in combined_model.parameters():
                param.data = param.data.to(dtype)
            
            for param in combined_likelihood.parameters():
                param.data = param.data.to(dtype)
            
            # Set the number of samples per protein ID
            samples_per_id = parameters['samples_per_id']
            
            for sample_idx in range(samples_per_id):
                with gpytorch.settings.prior_mode(True):
                    combined_model.eval()
                    combined_likelihood.eval()
                    prior_joint_model = combined_likelihood(combined_model(*test_x_joint_model))
                    prior_joint = prior_joint_model.sample().detach().cpu().numpy()  # Move to CPU for numpy
                
                # Get information about the sampled condition
                temp = combined_model.train_inputs
                temp_flat = temp[0]
                
                # Create condition assignments based on the known split point
                output_conds = []
                for j in range(len(temp_flat)):
                    if j < n_points_first_condition:
                        output_conds.append(conds[0])
                    else:
                        output_conds.append(conds[1])
                
                # Convert the list to a NumPy array
                output_conds_array = np.array(output_conds)
                
                # Make final dataframe for this sample
                prior_joint_df_sample = pd.DataFrame()
                prior_joint_df_sample.loc[:, 'uniqueID'] = np.repeat(f'{protein}_sample_{sample_idx + 1}', len(temp_flat))
                prior_joint_df_sample.loc[:, 'x'] = temp_flat.cpu().numpy()
                prior_joint_df_sample.loc[:, 'y'] = prior_joint
                prior_joint_df_sample.loc[:, 'condition'] = output_conds_array
                prior_joint_df_sample.loc[:, 'model'] = 'joint'
                prior_joint_df_sample.loc[:, 'comparison'] = f'{conds[0]} - {conds[1]}'
                prior_joint_df_sample.loc[:, 'data'] = 'joint model priors'
                
                all_samples.append(prior_joint_df_sample)
    
    # Convert the joint model and likelihood lists to IndependentModelList and LikelihoodList
    joint_model_list = gpytorch.models.IndependentModelList(*joint_models)
    joint_likelihood_list = gpytorch.likelihoods.LikelihoodList(*joint_likelihoods)
    
    # Update the checklist
    if is_notebook():
        if null_dataset == False:
            tasks[1][2][0] = (f"Combine training input and targets for {parameters['control_condition']} & {parameters['perturbation']}", True)
            update_checklist(tasks)
            tasks[1][2][1] = (f"Compute mll for joint model", True)
            update_checklist(tasks)
            tasks[1][2][2] = (f"Sample from the joint model (sampling of true negatives, null dataset)", True)
            update_checklist(tasks)
    
     
    # Save mll of joined model in dataframe
    mll_values_joint_model_df = pd.DataFrame({
        'protein': uniqueID_list,
        'condition': 'joint',
        'mll': mll_values
    })
    
    # Save parameters for joint model
    list_state_dict_joint = []
    for submodel in joint_model_list.models:
        list_state_dict_joint.append(submodel.state_dict())
    
    # Compute likelihood ratios for full vs joint model
    lr_values_full_vs_joint = compute_likelihood_ratio(
        mll_full_df=mll_values_full_model_df,
        mll_joint_df=mll_values_joint_model_df
    )
    
    # Create a dict with model dict and training information
    join_model_result_dict = deepcopy(result_dict)
    if null_dataset == False:
        join_model_result_dict.update({
            "joint_model_list": joint_model_list,
            "joint_likelihood_list": joint_likelihood_list,
            "joint_state_dict_list": list_state_dict_joint,
            "joint_mll_values": mll_values_joint_model_df,
            "lr_values_full_vs_joint": lr_values_full_vs_joint,
            "model_join_control_df": join_control_df,
            "sampled_prior_df": pd.concat(all_samples) if all_samples else pd.DataFrame()
        })
        # Mark the generation of a joint model and null dataset as complete
        if is_notebook():
            tasks[1][2][3] = (f"Save joint model, and null dataset", True)
            update_checklist(tasks)
            tasks[1] = ("2. Create a joint model and null dataset", True, tasks[1][2])
            update_checklist(tasks)
        
        # Save updates to summery file 
        with open(f'{output_path}/{pert}_Thermal_Tracks_summary.txt', 'a') as f:
            f.write("\n" + "="*70 + "\n")
            f.write("JOINT MODEL SUMMARY\n")
            f.write("="*70 + "\n")
            f.write(f"  Proteins processed: {len(proteins_to_process)}\n")
            f.write(f"  Joint models created: {len(joint_models)}\n")
            if null_dataset == False:
                f.write(f"  Samples generated: {len(all_samples)}\n")
            f.write(f"  Device used: {device}\n")
            f.write("="*70 + "\n\n")
    else:
        join_model_result_dict.update({
            "joint_model_list": joint_model_list,
            "joint_likelihood_list": joint_likelihood_list,
            "joint_state_dict_list": list_state_dict_joint,
            "joint_mll_values": mll_values_joint_model_df,
            "lr_values_full_vs_joint": lr_values_full_vs_joint,
            "model_join_control_df": join_control_df,
        })
        
        # Mark the generation of a joint model as complete
        if is_notebook():
            tasks[3][2][1] = ("Create a joint model", True)
            update_checklist(tasks)
    
    return join_model_result_dict
