from fit_model import *
from build_model import *
from join_model import define_joint_model_batched
from predict_model import *
from compute_lr_test_statistics import compute_likelihood_statistics
from prepare_results import prepare_gp_results_batched
import random

####################################################################
#                    Run Thermal-Tracks pipeline                   #
####################################################################

def run_ThermalTracks(tpptr_df_input, parameters):
    """
    Wrapper function to execute the Thermal-Tracks pipline
    
    Args:
        tpptr_df_input: Input dataframe with columns ['condition', 'uniqueID', 'x', 'y']
        parameters: Dictionary with required training parameters
    """
    import time
    
    start = time.time()

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
    live_update_message = "(Input data validation check...)"
    update_checklist(tasks, live_update={"task": "1. Build and fit full model", "message": live_update_message})

    # Check whether Input DataFrame is correctly formatted
    tpptr_df = tpptr_df_input.copy()
    required_columns = ['condition', 'uniqueID', 'x', 'y']
    missing_columns = [col for col in required_columns if col not in tpptr_df.columns]
    if missing_columns:
        example_df = pd.DataFrame({col: [] for col in required_columns})
        raise ValueError(f"Missing columns: {missing_columns}. Example of required DataFrame format:\n{example_df}")
    else:
        live_update_message = "(Required columns are present.)"
        update_checklist(tasks, live_update={"task": "1. Build and fit full model", "message": live_update_message})

    # Check that parameter file is correct
    required_keys = [
    # Directory and output settings
    "result_dir", 
    "perturbation",
    
    # Data and experimental settings
    "subset_test",
    "control_condition",
    "samples_per_id",
    
    # Device and resource management
    "device",
    "gpu_memory_fraction",
    "cpu_memory_fraction",
    "cpu_thread_strategy",
    "max_batch_size_gpu",
    "max_batch_size_cpu",
    "min_batch_size",
    
    # Model architecture
    "lengthscale_prior",
    "ScaledKernel",
    "MeanFunction",
    "lengthscale_minconstraint",
    "lengthscale_mult",
    "EarlyStop",
    
    # Training parameters
    "training_iterations",
    "learningRate",
    "amsgrad",
    
    # Prediction and output
    "n_predictions",
    "create_plots",
    "exclude_poor_fits"
]
    
    missing_keys = [key for key in required_keys if key not in parameters]
    if missing_keys:
        raise ValueError(f"Missing keys in parameters: {missing_keys}")
    
    # Ensure that there are at least two unique conditions to check
    unique_conditions = tpptr_df['condition'].unique()
    
    if len(unique_conditions) < 2:
        raise ValueError("Input DataFrame must have at least two unique conditions.")
    
    # Check whether the conditions match the parameter file
    for condition in unique_conditions[:2]:
        if condition not in parameters['control_condition'] and condition not in parameters['perturbation']:
            raise ValueError(
                f"Condition '{condition}' not found in parameters. "
                f"Specified conditions do not agree with parameter file."
            )
    
    live_update_message = "(Parameter file is ok!)"
    update_checklist(tasks, live_update={"task": "1. Build and fit full model", "message": live_update_message})    

    # Update checklist
    live_update_message = "(Input data validation complete.)"
    update_checklist(tasks, live_update={"task": "1. Build and fit full model", "message": live_update_message})
    

    # Check whether pipeline should be run on a subset of data first
    if parameters.get('subset_test', False):
        live_update_message = "(SUBSET TEST MODE: Using 200 random proteins.)"
        update_checklist(tasks, live_update={"task": "1. Build and fit full model", "message": live_update_message})
        # Generate a list with 200 proteins at random
        random.seed(42)
        unique_ids = tpptr_df['uniqueID'].unique().tolist()
        
        # Limit to 200 if more are available
        n_subset = min(200, len(unique_ids))
        random_ids = random.sample(unique_ids, n_subset)

        # Produce a reduced TPP-TR dataset with the subset proteins
        tpptr_df = tpptr_df[tpptr_df['uniqueID'].isin(random_ids)].copy()
        print(f"  Selected {n_subset} proteins for testing\n")
    
    # Build and train full model
    train_full_dict = train_model_batched_optimized(tpptr_df, parameters, null_dataset=False)

    # Create a joint model (combination of two conditions) and sample from it to generate a null dataset
    joint_full_dict =  define_joint_model_batched(train_full_dict, parameters, null_dataset=False)

    # Evaluate and predict
    predict_full_dict = predict_and_evaluate_batched(joint_full_dict, parameters, null_dataset = False)

    # Build and train null model
    sampled_df = predict_full_dict['sampled_prior_df'].copy()
    train_null_dict = train_model_batched_optimized(sampled_df, parameters, null_dataset = True)

    # Create a joint model (combination of two conditions) from null model - no sampling
    joint_null_dict =  define_joint_model_batched(train_null_dict, parameters, null_dataset=True)

    # Evaluate and predict null model
    predict_null_dict = predict_and_evaluate_batched(joint_null_dict, parameters, null_dataset = True)

    # Compute LR test statistics
    lr_test_stats_dict = compute_likelihood_statistics(predict_full_dict, predict_null_dict, parameters)

    # Combine and save results
    prepare_gp_results_batched(lr_test_stats_dict, predict_null_dict, parameters)
    
    elapsed_time = time.time() - start

    # Save Task Summary
    output_path = str(parameters['result_dir'])
    pert = parameters['perturbation']
    with open(f'{output_path}/{pert}_Thermal_Tracks_summary.txt', 'a') as f:
        f.write("\n" + "="*70 + "\n")
        f.write("Thermal Tracks Task Complete\n")
        f.write("="*70 + "\n")
        f.write(f"  Total elapsed time: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)\n")
        f.write(f"  Proteins trained: {len(train_full_dict['proteins'])}\n")
        f.write(f"  Conditions: {len(train_full_dict['conditions'])}\n")
        f.write(f"  Total tasks: {train_full_dict['n_total_tasks']}\n")
        f.write(f"  Batch size used: {train_full_dict['batch_size']}\n")
        f.write(f"  Device: {train_full_dict['device']}\n")
        f.write("="*70 + "\n")