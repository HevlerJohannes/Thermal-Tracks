#!/usr/bin/env python3
"""
Thermal Tracks - Command-line interface for batched GP training pipeline

Usage:
    python main.py -i input_data.csv -p parameters.json
    python main.py --input tpp_ecoli_mgcl2.csv --params params.json
    python main.py -i data.csv -p params.json -o ./results
"""

import argparse
import json
import sys
import os
import pandas as pd
import random
import time
from pathlib import Path

# Get absolute path of this script
script_dir = os.path.dirname(os.path.abspath(__file__))

# If we're inside 'src', use that; if not, look for it one level up
if os.path.basename(script_dir) == "src":
    src_dir = script_dir
else:
    src_dir = os.path.join(script_dir, "src")

# Add src_dir to Python path if not already included
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)

# Import util functions
from utils import *
# Import gp modules
from fit_model import *
from build_model import *
from join_model import define_joint_model_batched
from predict_model import *
from compute_lr_test_statistics import compute_likelihood_statistics
from prepare_results import prepare_gp_results_batched

# Ignore all warnings
import warnings
warnings.filterwarnings("ignore")

####################################################################
#                    Run ThermalTracks pipeline                    #
####################################################################

def run_ThermalTracks(tpptr_df_input, parameters):
    """
    Wrapper function to test batched GP training speed with input validation

    Args:
        tpptr_df_input: Input dataframe with columns ['condition', 'uniqueID', 'x', 'y']
        parameters: Dictionary with required training parameters

    Returns:
        train_full_dict: Training results dictionary
        elapsed_time: Total elapsed time in seconds
    """
    start = time.time()

    # Check that data input is correctly formatted and that parameters are correct
    print("Validating input data...")

    # Check whether Input DataFrame is correctly formatted
    tpptr_df = tpptr_df_input.copy()
    required_columns = ['condition', 'uniqueID', 'x', 'y']
    missing_columns = [col for col in required_columns if col not in tpptr_df.columns]
    if missing_columns:
        example_df = pd.DataFrame({col: [] for col in required_columns})
        raise ValueError(f"Missing columns: {missing_columns}. Example of required DataFrame format:\n{example_df}")
    else:
        print("Required columns are present.")

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

    print("Parameter file is ok!")
    print("Input data validation complete.")

    # Check whether pipeline should be run on a subset of data first
    if parameters.get('subset_test', False):
        print("SUBSET TEST MODE: Using 200 random proteins")
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
    print("\n" + "="*70)
    print("Step 1: Building and training full model")
    print("="*70)
    train_full_dict = train_model_batched_optimized(tpptr_df, parameters, null_dataset=False)

    # Create a joint model (combination of two conditions) and sample from it to generate a null dataset
    print("\n" + "="*70)
    print("Step 2: Creating joint model and null dataset")
    print("="*70)
    joint_full_dict = define_joint_model_batched(train_full_dict, parameters, null_dataset=False)

    # Evaluate and predict
    print("\n" + "="*70)
    print("Step 3: Evaluating and predicting models")
    print("="*70)
    predict_full_dict = predict_and_evaluate_batched(joint_full_dict, parameters, null_dataset=False)

    # Build and train null model
    print("\n" + "="*70)
    print("Step 4: Building and training null model")
    print("="*70)
    sampled_df = predict_full_dict['sampled_prior_df'].copy()
    train_null_dict = train_model_batched_optimized(sampled_df, parameters, null_dataset=True)

    # Create a joint model (combination of two conditions) from null model - no sampling
    print("\n" + "="*70)
    print("Step 5: Creating joint model from null dataset")
    print("="*70)
    joint_null_dict = define_joint_model_batched(train_null_dict, parameters, null_dataset=True)

    # Evaluate and predict null model
    print("\n" + "="*70)
    print("Step 6: Evaluating and predicting null model")
    print("="*70)
    predict_null_dict = predict_and_evaluate_batched(joint_null_dict, parameters, null_dataset=True)

    # Compute LR test statistics
    print("\n" + "="*70)
    print("Step 7: Computing likelihood ratio test statistics")
    print("="*70)
    lr_test_stats_dict = compute_likelihood_statistics(predict_full_dict, predict_null_dict, parameters)

    # Combine and save results
    print("\n" + "="*70)
    print("Step 8: Preparing and saving results")
    print("="*70)
    prepare_gp_results_batched(lr_test_stats_dict, predict_null_dict, parameters)

    elapsed_time = time.time() - start

    # Save Task Summary
    output_path = str(parameters['result_dir'])
    pert = parameters['perturbation']
    with open(f'{output_path}/{pert}_Thermal_Tracks_summary.txt', 'a') as f:
        f.write("\n" + "="*70 + "\n")
        f.write("Thermal Tracks Task Complete\n")
        f.write("="*70 + "\n")
        f.write(f"Total elapsed time: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)\n")
        f.write(f"Proteins trained: {len(train_full_dict['proteins'])}\n")
        f.write(f"Conditions: {len(train_full_dict['conditions'])}\n")
        f.write(f"Total tasks: {train_full_dict['n_total_tasks']}\n")
        f.write(f"Batch size used: {train_full_dict['batch_size']}\n")
        f.write(f"Device: {train_full_dict['device']}\n")
        f.write("="*70 + "\n")

    print("\n" + "="*70)
    print("THERMAL TRACKS TASK COMPLETE")
    print("="*70)
    print(f"Total elapsed time: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
    print(f"Results saved to: {output_path}")
    print("="*70)

    return elapsed_time


def load_input_data(filepath):
    """
    Load input CSV file

    Args:
        filepath: Path to the input CSV file

    Returns:
        DataFrame with the input data
    """
    try:
        df = pd.read_csv(filepath)
        print(f"Loaded input data from: {filepath}")
        print(f"Shape: {df.shape}")
        print(f"Columns: {list(df.columns)}")
        return df
    except FileNotFoundError:
        raise FileNotFoundError(f"Input file not found: {filepath}")
    except pd.errors.EmptyDataError:
        raise ValueError(f"Input file is empty: {filepath}")
    except Exception as e:
        raise Exception(f"Error loading input file: {e}")


def load_parameters(filepath):
    """
    Load parameter JSON file

    Args:
        filepath: Path to the parameter JSON file

    Returns:
        Dictionary with parameters
    """
    try:
        with open(filepath, 'r') as f:
            params = json.load(f)
        print(f"Loaded parameters from: {filepath}")
        return params
    except FileNotFoundError:
        raise FileNotFoundError(f"Parameter file not found: {filepath}")
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON in parameter file: {e}")
    except Exception as e:
        raise Exception(f"Error loading parameter file: {e}")


def validate_output_dir(output_dir):
    """
    Validate and create output directory if it doesn't exist

    Args:
        output_dir: Path to the output directory

    Returns:
        Path object for the output directory
    """
    output_path = Path(output_dir)
    if not output_path.exists():
        try:
            output_path.mkdir(parents=True, exist_ok=True)
            print(f"Created output directory: {output_dir}")
        except Exception as e:
            raise Exception(f"Could not create output directory: {e}")
    else:
        print(f"Using existing output directory: {output_dir}")
    return output_path


def main():
    """Main function to handle command-line execution"""

    banner = """
        ╔════════════════════════════════════════════════════════════════╗
        ║                                                                ║
        ║  ████████╗██╗  ██╗███████╗██████╗ ███╗   ███╗ █████╗ ██╗       ║
        ║  ╚══██╔══╝██║  ██║██╔════╝██╔══██╗████╗ ████║██╔══██╗██║       ║
        ║     ██║   ███████║█████╗  ██████╔╝██╔████╔██║███████║██║       ║
        ║     ██║   ██╔══██║██╔══╝  ██╔══██╗██║╚██╔╝██║██╔══██║██║       ║
        ║     ██║   ██║  ██║███████╗██║  ██║██║ ╚═╝ ██║██║  ██║███████╗  ║
        ║     ╚═╝   ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚═╝     ╚═╝╚═╝  ╚═╝╚══════╝  ║
        ║                                                                ║
        ║  ████████╗██████╗  █████╗  ██████╗██╗  ██╗███████╗             ║
        ║  ╚══██╔══╝██╔══██╗██╔══██╗██╔════╝██║ ██╔╝██╔════╝             ║
        ║     ██║   ██████╔╝███████║██║     █████╔╝ ███████╗             ║
        ║     ██║   ██╔══██╗██╔══██║██║     ██╔═██╗ ╚════██║             ║
        ║     ██║   ██║  ██║██║  ██║╚██████╗██║  ██╗███████║             ║
        ║     ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝╚══════╝             ║
        ║                                                                ║
        ╠════════════════════════════════════════════════════════════════╣
        ║  Copyright (c) 2025 Johannes F. Hevler                         ║
        ╚════════════════════════════════════════════════════════════════╝
        """
    print(banner)
    # Create argument parser
    parser = argparse.ArgumentParser(
        description='Thermal Tracks - Batched GP training pipeline for thermal proteome profiling data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""

        Example usage:
  python main.py -i tpp_ecoli_mgcl2.csv -p params.json
  python main.py --input data.csv --params parameters.json --output ./results
  python main.py -i input.csv -p config.json -v

Input CSV format:
  Required columns: condition, uniqueID, x, y
  - condition: experimental condition label
  - uniqueID: unique identifier for each protein
  - x: temperature values
  - y: intensity values
        """
    )

    # Add arguments
    parser.add_argument(
        '-i', '--input',
        type=str,
        required=True,
        help='Path to input CSV file with TPP data'
    )

    parser.add_argument(
        '-p', '--params',
        type=str,
        required=True,
        help='Path to parameter JSON file'
    )

    parser.add_argument(
        '-o', '--output',
        type=str,
        default=None,
        help='Output directory path (overrides parameter file setting)'
    )

    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output'
    )

    parser.add_argument(
        '--version',
        action='version',
        version='Thermal-Tracks v1.0.0'
    )

    # Parse arguments
    args = parser.parse_args()

    # Print header
    print("\n" + "="*70)
    print("THERMAL TRACKS")
    print("="*70)

    try:
        # Load input data
        print("\n1. Loading input data...")
        input_df = load_input_data(args.input)

        # Load parameters
        print("\n2. Loading parameters...")
        parameters = load_parameters(args.params)

        # Replace Boolean
        for key in parameters.keys():
            if str(parameters[key]).lower() == "true":
                parameters[key] = True
            elif str(parameters[key]).lower() == "false":
                parameters[key] = False

        # Override output directory if specified
        if args.output:
            print(f"\n3. Overriding output directory to: {args.output}")
            parameters['result_dir'] = args.output

        # Validate and create output directory
        print(f"\n3. Setting up output directory...")
        output_path = validate_output_dir(parameters['result_dir'])
        parameters['result_dir'] = str(output_path)

        # Properly configure parameters
        ## Set dtype
        if parameters['dtype'] == 64:
            parameters['dtype'] = torch.float64
        elif parameters['dtype'] == 32:
            parameters['dtype'] = torch.float32
        else:
            parameters['dtype'] = torch.float64

        ## Training parameters
        if parameters['MeanFunction'] == 'zero':
            parameters['MeanFunction'] = gpytorch.means.ZeroMean()
        else:
            parameters['MeanFunction'] = gpytorch.means.ConstantMean()

        if parameters['lengthscale_prior'] is not None:
           lengthscale_prior_value = parameters['lengthscale_prior']
           parameters['lengthscale_prior'] = gpytorch.priors.GammaPrior(lengthscale_prior_value, 1)
        else:
            parameters['lengthscale_prior'] = None        

        # Print summary
        print("\n" + "="*70)
        print("CONFIGURATION SUMMARY")
        print("="*70)
        print(f"Input file: {args.input}")
        print(f"Parameter file: {args.params}")
        print(f"Output directory: {parameters['result_dir']}")
        print(f"Perturbation: {parameters.get('perturbation', 'Not specified')}")
        print(f"Control condition: {parameters.get('control_condition', 'Not specified')}")
        print(f"Device: {parameters.get('device', 'auto')}")
        print(f"Training iterations: {parameters.get('training_iterations', 'Not specified')}")
        print(f"Create plots: {parameters.get('create_plots', False)}")
        print("="*70)

        # Run the pipeline
        print("\nStarting Thermal Tracks pipeline...")
        elapsed_time = run_ThermalTracks(input_df, parameters)

        # Success message
        print("\n" + "="*70)
        print("PIPELINE COMPLETED SUCCESSFULLY")
        print("="*70)
        print(f"Total runtime: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
        print(f"Output saved to: {parameters['result_dir']}")
        print("="*70 + "\n")

        return 0

    except FileNotFoundError as e:
        print(f"\n❌ ERROR: {e}", file=sys.stderr)
        return 1
    except ValueError as e:
        print(f"\n❌ ERROR: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"\n❌ UNEXPECTED ERROR: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())