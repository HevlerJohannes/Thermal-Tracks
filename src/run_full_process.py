#! /usr/bin/env python3
import pandas as pd
import random

# Import util functions
from utils import *
# Import gp modules
from fit_model import train_model
#from join_model import define_joint_model
from join_model_deep_sampling import define_joint_model
from predict_model import predict_and_evaluate
from compute_lr_test_statistics import compute_likelihood_statistics
from prepare_results import prepare_gp_results

####################################################################
#                    Run HeatBeats pipeline                        #
####################################################################

def full_gp_process(tpptr_df_input, parameters):

    # Print checklist
    tasks = [
        ("1. Build and fit full model", False),
        ("2. Creating a joint model and null dataset", False),
        ("3. Evaluate and predict models", False),
        ("4. Build and fit null model", False),
        ("5. Compute likelihood ratio test statistics", False),
        ("6. Combine and create result files", False)]
            
    # Update the checklist
    update_checklist(tasks)

    # Check that data input is correctly formatted and that parameters are correct
    live_update_message = f"(Input data validation check...)"
    update_checklist(tasks, live_update={"task": "1. Build and fit full model", "message": live_update_message})

    # Check whether Input DataFrame is correctly formatted
    tpptr_df = tpptr_df_input.copy()
    required_columns = ['condition', 'uniqueID', 'x', 'y']
    missing_columns = [col for col in required_columns if col not in tpptr_df.columns]
    if missing_columns:
         example_df = pd.DataFrame({col: [] for col in required_columns})
         raise ValueError(f"Missing columns: {missing_columns}. Example of required DataFrame format:\n{example_df}")
    else:
        live_update_message = f"(Required columns are present.)"
        update_checklist(tasks, live_update={"task": "1. Build and fit full model", "message": live_update_message})

    # Check that parameter file is correct
    required_keys = [
    "result_dir", "subset_test", "lengthscale_prior", "lengthscale_minconstraint", 
    "control_condition", "perturbation", "training_iterations", "learningRate", 
    "amsgrad", "n_predictions", "create_plots", "exclude_poor_fits", "samples_per_id"]
    
    missing_keys = [key for key in required_keys if key not in parameters]
    if missing_keys:
            raise ValueError(f"Missing keys: {missing_keys}")
    else:

        # Ensure that there are at least two unique conditions to check
        unique_conditions = tpptr_df['condition'].unique()

        assert len(unique_conditions) >= 2, "Input DataFrame must have at least two unique conditions."

        # Check whether the first two unique conditions are in the control or perturbation conditions
        for condition in unique_conditions[:2]:
            assert condition in parameters['control_condition'] or condition in parameters['perturbation'], \
                f"Condition '{condition}' Specified conditions: Input dataframe does not agree with the parameter file."
            
        live_update_message = f"(Parameter file is ok!)"
        update_checklist(tasks, live_update={"task": "1. Build and fit full model", "message": live_update_message})    

    # Update checklist
    live_update_message = f"(Input data validation complete.)"
    update_checklist(tasks, live_update={"task": "1. Build and fit full model", "message": live_update_message})

    # Check whether pipline should be run on a suubset of data first
    if parameters['subset_test'] == True:
        # Generate a list with 200 proteins at random
        random.seed(42)
        unique_ids = tpptr_df['uniqueID'].unique().tolist()
        random_ids = random.sample(unique_ids, 200)

        # Produce a reduced TPP-TR dataset with the 200 proteins
        tpptr_df = tpptr_df[tpptr_df['uniqueID'].isin(random_ids)].copy()
    else:
        # Skip subset creation, proceed as is
        pass
    
    # Build and train full model
    train_full_dict = train_model(tpptr_df, parameters, null_dataset = False)

    # Create a joint model (combination of two conditions) and sample from it to generate a null dataset
    joint_full_dict =  define_joint_model(train_full_dict, parameters, null_dataset=False)

    # Evaluate and predict
    predict_full_dict = predict_and_evaluate(joint_full_dict, parameters, null_dataset = False)

    # Build and train null model
    sampled_df = predict_full_dict['sampled_prior_df'].copy()
    train_null_dict = train_model(sampled_df, parameters, null_dataset = True)

    # Create a joint model (combination of two conditions) from null model - no sampling
    joint_null_dict =  define_joint_model(train_null_dict, parameters, null_dataset=True)

    # Evaluate and predict null model
    predict_null_dict = predict_and_evaluate(joint_null_dict, parameters, null_dataset = True)

    # Compute LR test statistics
    lr_test_stats_dict = compute_likelihood_statistics(predict_full_dict, predict_null_dict, parameters)

    # Combine and save results
    prepare_gp_results(lr_test_stats_dict, predict_null_dict, parameters)
    
     

 