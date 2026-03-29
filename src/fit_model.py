import torch
import gpytorch
import pandas as pd
import matplotlib
import csv
import os
import numpy as np
from gpytorch.mlls import ExactMarginalLogLikelihood
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
import time

# Import util functions
from utils import update_checklist, is_notebook

if is_notebook():
    # Notebook settings
    import matplotlib.pyplot as plt
else:
    # Terminal settings
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

# Import GP model building functions
from build_model import (
    calculate_optimal_batch_size,
    building_exactgp_model_batched_optimized,
)


def compute_batch_mll_parallel(model_data):
    """
    Compute per-task marginal log-likelihood for one batched model.

    Designed for use with ThreadPoolExecutor. Handles both scalar MLL
    (single-task batches where the MLL is divided equally) and batched MLL
    tensors (multi-task batches with per-task values).

    Args:
        model_data: Tuple of (model, likelihood, metadata_list) as stored
            in model_list.

    Returns:
        List of dicts with keys 'protein', 'condition', 'mll'.
    """
    model, likelihood, meta_chunk = model_data
    train_x = model.train_inputs[0]
    train_y = model.train_targets

    output = model(train_x)
    mll = ExactMarginalLogLikelihood(likelihood, model)

    mll_results = []

    if model.num_tasks >= 1:
        mll_batch = mll(output, train_y)

        if mll_batch.ndim == 0:  # Scalar
            per_task_mll = mll_batch.item() / model.num_tasks
            for meta in meta_chunk:
                mll_results.append(
                    {
                        "protein": meta["protein"],
                        "condition": meta["condition"],
                        "mll": per_task_mll,
                    }
                )
        else:  # Batched [num_tasks]
            for task_idx, meta in enumerate(meta_chunk):
                mll_results.append(
                    {
                        "protein": meta["protein"],
                        "condition": meta["condition"],
                        "mll": mll_batch[task_idx].item(),
                    }
                )
    else:
        meta = meta_chunk[0]
        mll_batch = mll(output, train_y)
        mll_results.append(
            {
                "protein": meta["protein"],
                "condition": meta["condition"],
                "mll": mll_batch.item(),
            }
        )

    return mll_results


def train_model_batched_optimized(tpptr_df_input, parameters, null_dataset=False):
    """
    Train all batched GP models with a single Adam optimizer.

    Handles model building, iterative training with optional early stopping,
    per-task MLL computation, and result serialization. When the 'monotonicity_penalty'
    parameter is set (> 0), an additional regularization term penalizes non-decreasing
    posterior means, encouraging the expected melting-curve shape (high at low temp,
    low at high temp).

    Args:
        tpptr_df_input: Input dataframe with columns 'uniqueID', 'condition', 'x', 'y'.
        parameters: Full parameter dictionary. Key training-related entries:
            training_iterations, learningRate, amsgrad, EarlyStop,
            monotonicity_penalty (float, default 0.0).
        null_dataset: If True, trains on the sampled null data instead of observed data.

    Returns:
        train_results_dict: Dictionary containing trained model_list, state dicts,
            MLL dataframe, metadata, and run statistics.
    """

    # Get parallel processing parameters
    use_parallel = parameters.get("use_parallel_processing", True)
    n_workers = parameters.get("n_parallel_workers", 4)

    # Setup tasks for progress tracking (same as before)
    if is_notebook():
        if null_dataset == False:
            tasks = [
                (
                    "1. Build and fit full model",
                    False,
                    [
                        ("Prepare input for GP process", False),
                        ("Calculate optimal batch size", False),
                        ("Build models", False),
                        ("Initialize models", False),
                        ("Train models", False),
                        ("Compute mll for full models", False),
                        ("Save models and training results", False),
                    ],
                ),
                ("2. Creating a joint model and null dataset", False),
                ("3. Evaluate and predict models", False),
                ("4. Build and fit null model", False),
                ("5. Compute likelihood ratio test statistics", False),
                ("6. Combine and create result files", False),
            ]
        else:
            tasks = [
                ("1. Build and fit full model", True),
                ("2. Create a joint model and null dataset", True),
                ("3. Evaluate and predict full and joint models", True),
                (
                    "4. Build and fit null model",
                    False,
                    [
                        ("Train null model", False),
                        ("Create a joint model", False),
                        ("Evaluate and predict null model", False),
                    ],
                ),
                ("5. Compute likelihood ratio test statistics", False),
                ("6. Combine and create result files", False),
            ]

        update_checklist(tasks)

    ############################## PREPARATION ##############################
    # Prepare input
    tpptr_df = tpptr_df_input.copy()
    pert = parameters["perturbation"]
    control = parameters["control_condition"]
    output_path = str(parameters["result_dir"])
    proteins2test = tpptr_df["uniqueID"].unique()
    conds = tpptr_df["condition"].unique()
    scale_kernel = parameters["ScaledKernel"]
    mean = parameters["MeanFunction"]
    dtype = parameters["dtype"]

    if is_notebook():
        # Update checklist
        if null_dataset == False:
            tasks[0][2][0] = ("Prepare input for GP process", True)
            update_checklist(tasks)

    # Training parameters
    lengthscale_prior = parameters["lengthscale_prior"]
    lengthscale_minconstraint = parameters["lengthscale_minconstraint"]
    lengthscale_mult = parameters["lengthscale_mult"]
    n_iterations = parameters["training_iterations"]

    # Calculate number of datapoints (for batch size calculation)
    sample_protein = proteins2test[0]
    sample_cond = tpptr_df[tpptr_df["uniqueID"] == sample_protein][
        "condition"
    ].unique()[0]
    sample_df = tpptr_df[
        (tpptr_df["uniqueID"] == sample_protein)
        & (tpptr_df["condition"] == sample_cond)
    ]
    n_datapoints = len(sample_df)

    ############################## CALCULATE OPTIMAL BATCH SIZE ##############################
    # Determine device and batch size
    device_preference = parameters.get("device", "auto")
    batch_size, use_gpu, torch_device = calculate_optimal_batch_size(
        parameters, n_datapoints, device=device_preference, null_dataset=null_dataset
    )

    device = torch.device(torch_device)

    # Update checklist
    if is_notebook():
        if null_dataset == False:
            tasks[0][2][1] = ("Calculate optimal batch size", True)
            update_checklist(tasks)

    ############################## BUILD MODELS ##############################
    start_time = time.time()

    # Use the optimized model building function
    model_list, n_cond, n_prot, n_models = building_exactgp_model_batched_optimized(
        tpptr_df,
        parameters,
        proteins2test,
        conds,
        scale_kernel,
        lengthscale_prior,
        lengthscale_minconstraint,
        lengthscale_mult,
        mean,
        batch_size=batch_size,
        device=device,
        null_dataset=null_dataset,
        use_parallel=use_parallel,
        n_workers=n_workers,
        dtype=dtype,
    )

    build_time = time.time() - start_time

    # Update checklist
    if is_notebook():
        if null_dataset == False:
            tasks[0][2][2] = ("Build models", True)
            update_checklist(tasks)

    # Flatten metadata list for easier access
    metadata_list = []
    for model, likelihood, meta_chunk in model_list:
        metadata_list.extend(meta_chunk)

    ############################## TRAIN MODELS ##############################
    # Set all models to training mode
    for model, likelihood, _ in model_list:
        model.train()
        likelihood.train()

    # Update checklist
    if is_notebook():
        if null_dataset == False:
            tasks[0][2][3] = ("Initialize models", True)
            tasks[0][2][4] = ("Train models", False)
            update_checklist(tasks)

    train_start_time = time.time()

    # Collect all parameters
    all_params = []
    for model, likelihood, _ in model_list:
        all_params.extend(list(model.parameters()))
        all_params.extend(list(likelihood.parameters()))

    # Optimizer
    optimizer = torch.optim.Adam(
        all_params, lr=parameters["learningRate"], amsgrad=parameters["amsgrad"]
    )

    # Pre-compute MLL objects for each model
    mll_objects = []
    for model, likelihood, _ in model_list:
        mll_objects.append(ExactMarginalLogLikelihood(likelihood, model))

    # Training loop
    training_iterations = n_iterations
    LossValues = []
    training_info = []

    for i in range(training_iterations):
        optimizer.zero_grad()

        # Vectorized loss computation
        losses = []
        task_counts = []

        # Compute losses for all models
        for idx, (model, likelihood, meta_chunk) in enumerate(model_list):
            train_x = model.train_inputs[0]
            train_y = model.train_targets

            # Forward pass
            output = model(train_x)

            # Use pre-computed MLL object
            mll_values = mll_objects[idx](output, train_y)

            # Collect loss and task count
            losses.append(-mll_values.sum())
            task_counts.append(model.num_tasks)

        # Stack and normalize losses efficiently
        total_loss = torch.stack(losses).sum()
        total_tasks = sum(task_counts)
        normalized_loss = total_loss / total_tasks

        # Monotonicity penalty: curves must only decrease as temperature increases
        mono_weight = parameters.get("monotonicity_penalty", 0.0)
        if mono_weight > 0:
            mono_penalty = torch.tensor(0.0, device=device, dtype=dtype)
            for model, likelihood, meta_chunk in model_list:
                # Switch to eval to get posterior (fitted curve)
                model.eval()
                likelihood.eval()
                train_x = model.train_inputs[0]

                # Evaluate at dense grid (catches waviness between data points)
                dense_x = torch.linspace(
                    train_x.min().item(),
                    train_x.max().item(),
                    100,
                    device=device,
                    dtype=dtype,
                )
                with gpytorch.settings.fast_pred_var():
                    posterior = model(dense_x)
                pred_mean = posterior.mean

                # Switch back to training mode
                model.train()
                likelihood.train()

                if pred_mean.dim() == 1:
                    pred_mean = pred_mean.unsqueeze(0)

                # 1) Monotonic decrease: LINEAR penalty (constant gradient signal)
                #    even tiny increases get a strong push back
                diffs = pred_mean[:, 1:] - pred_mean[:, :-1]
                mono_penalty = mono_penalty + torch.relu(diffs).sum()

                # 2) High at low temp: first point should be near 1.0
                mono_penalty = mono_penalty + torch.relu(0.8 - pred_mean[:, 0]).sum()

                # 3) Low at high temp: last point should be near 0.0
                mono_penalty = mono_penalty + torch.relu(pred_mean[:, -1] - 0.2).sum()

            normalized_loss = normalized_loss + mono_weight * mono_penalty / total_tasks

        # Backward pass
        normalized_loss.backward()

        optimizer.step()

        LossValues.append(normalized_loss.item())

        # Live update for checklist (reduce overhead not every epoch is reported)
        if is_notebook():
            if (
                i % max(1, training_iterations // 20) == 0
                or i == training_iterations - 1
            ):
                if null_dataset == False:
                    live_update_message = f"(Training iteration {i + 1}/{n_iterations} - Loss: {normalized_loss.item():.6f})"
                    plot_data = {
                        "task": "1. Build and fit full model",
                        "subtask": "Train models",
                        "x": list(range(1, len(LossValues) + 1)),
                        "y": LossValues,
                        "label": "Loss",
                        "xlabel": "Iteration",
                        "ylabel": "Loss",
                        "title": "Training Loss Over Iterations",
                    }
                    update_checklist(
                        tasks,
                        live_update={
                            "task": "1. Build and fit full model",
                            "subtask": "Train models",
                            "message": live_update_message,
                        },
                        plot_data=plot_data,
                    )
                else:
                    live_update_message = f"(Training iteration {i + 1}/{n_iterations} - Loss: {normalized_loss.item():.6f})"
                    update_checklist(
                        tasks,
                        live_update={
                            "task": "4. Build and fit null model",
                            "subtask": "Train null model",
                            "message": live_update_message,
                        },
                    )

        # Early stopping check
        value = parameters.get("EarlyStop", 0.0001)
        if is_notebook():
            if i > 1 and abs(LossValues[i] - LossValues[i - 1]) <= value:
                live_update_message = f"(Early stopping at iteration {i + 1}"

                # Create fresh plot_data for early stopping
                early_stop_plot_data = {
                    "task": "1. Build and fit full model",
                    "subtask": "Train models",
                    "x": list(range(1, len(LossValues) + 1)),
                    "y": LossValues,
                    "label": "Loss",
                    "xlabel": "Iteration",
                    "ylabel": "Loss",
                    "title": "Training Loss Over Iterations",
                }
                update_checklist(
                    tasks,
                    live_update={
                        "task": "1. Build and fit full model",
                        "subtask": "Train models",
                        "message": live_update_message,
                    },
                    plot_data=early_stop_plot_data,
                )
                break
            else:
                if i > 1 and abs(LossValues[i] - LossValues[i - 1]) <= value:
                    print(f"(Early stopping at iteration {i + 1}")
                    break

        # Collect training info
        for model_idx, (model, likelihood, meta_chunk) in enumerate(model_list):
            num_tasks_in_model = model.num_tasks

            for task_idx in range(num_tasks_in_model):
                # Extract parameters based on kernel type
                if isinstance(model.covar_module, gpytorch.kernels.ScaleKernel):
                    if num_tasks_in_model >= 1:
                        lengthscale = model.covar_module.base_kernel.lengthscale[
                            task_idx
                        ].item()
                        outputscale = model.covar_module.outputscale[task_idx].item()
                    else:
                        lengthscale = model.covar_module.base_kernel.lengthscale.item()
                        outputscale = model.covar_module.outputscale.item()
                else:
                    if num_tasks_in_model >= 1:
                        lengthscale = model.covar_module.lengthscale[task_idx].item()
                    else:
                        lengthscale = model.covar_module.lengthscale.item()
                    outputscale = 1.0

                # Get noise
                try:
                    if num_tasks_in_model >= 1:
                        noise = likelihood.noise[task_idx].item()
                    else:
                        noise = likelihood.noise.item()
                except:
                    noise = float("nan")

                # Get metadata for this task
                if task_idx < len(meta_chunk):
                    meta = meta_chunk[task_idx]
                    training_info.append(
                        [
                            i + 1,
                            normalized_loss.item(),
                            meta["global_index"] + 1,
                            lengthscale,
                            outputscale,
                            noise,
                            meta["protein"],
                            meta["condition"],
                        ]
                    )

    train_time = time.time() - train_start_time

    # Update checklist
    if is_notebook():
        live_update_message = (
            f"Training complete after {len(LossValues)} iterations ({train_time:.1f}s)"
        )
        update_checklist(
            tasks,
            live_update={
                "task": "1. Build and fit full model",
                "subtask": "Train models",
                "message": live_update_message,
            },
        )
        if null_dataset == False:
            tasks[0][2][4] = ("Train models", True)
            update_checklist(tasks)

    ############################## COMPUTE MLL ##############################
    mll_start_time = time.time()

    if use_parallel and len(model_list) > 1:
        # Parallel MLL computation
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            mll_results_nested = list(
                executor.map(compute_batch_mll_parallel, model_list)
            )

        # Flatten results
        mll_values = []
        for batch_results in mll_results_nested:
            mll_values.extend(batch_results)
    else:
        # Sequential MLL computation (original code)
        mll_values = []
        for model, likelihood, meta_chunk in model_list:
            train_x = model.train_inputs[0]
            train_y = model.train_targets

            output = model(train_x)
            mll = ExactMarginalLogLikelihood(likelihood, model)

            if model.num_tasks >= 1:
                mll_batch = mll(output, train_y)

                if mll_batch.ndim == 0:  # Scalar
                    per_task_mll = mll_batch.item() / model.num_tasks
                    for meta in meta_chunk:
                        mll_values.append(
                            {
                                "protein": meta["protein"],
                                "condition": meta["condition"],
                                "mll": per_task_mll,
                            }
                        )
                else:  # Batched [num_tasks]
                    for task_idx, meta in enumerate(meta_chunk):
                        mll_values.append(
                            {
                                "protein": meta["protein"],
                                "condition": meta["condition"],
                                "mll": mll_batch[task_idx].item(),
                            }
                        )
            else:
                meta = meta_chunk[0]
                mll_batch = mll(output, train_y)
                mll_values.append(
                    {
                        "protein": meta["protein"],
                        "condition": meta["condition"],
                        "mll": mll_batch.item(),
                    }
                )

    mll_values_df = pd.DataFrame(mll_values)
    mll_time = time.time() - mll_start_time

    # Update checklist
    if is_notebook():
        if null_dataset == False:
            tasks[0][2][5] = ("Compute mll for full models", True)
            update_checklist(tasks)

    ############################## SAVE RESULTS ##############################

    # Create training info dataframe
    if training_info:
        training_info_df = pd.DataFrame(
            training_info,
            columns=[
                "iteration",
                "loss",
                "model",
                "lengthscale",
                "outputscale",
                "noise",
                "protein",
                "condition",
            ],
        )
    else:  # create an empty dataframe
        training_info_df = pd.DataFrame(
            columns=[
                "iteration",
                "loss",
                "model",
                "lengthscale",
                "outputscale",
                "noise",
                "protein",
                "condition",
            ]
        )

    # Generate convergence plot
    plt.figure(figsize=(10, 6))
    plt.plot(range(len(LossValues)), LossValues, linewidth=2)
    plt.xlabel("Number of iterations", fontsize=12)
    plt.ylabel("Loss", fontsize=12)
    title_suffix = "Full Model" if not null_dataset else "Null Model"
    plt.title(f"{title_suffix} training - {pert}", fontsize=14)
    plt.grid(True, alpha=0.3)

    if null_dataset == False:
        plt.savefig(
            f"{output_path}/loss_full_model_gp_{control}_{pert}.pdf",
            bbox_inches="tight",
            dpi=300,
        )
    else:
        plt.savefig(
            f"{output_path}/loss_null_model_gp_{control}_{pert}.pdf",
            bbox_inches="tight",
            dpi=300,
        )
    plt.close()

    # Save loss values
    loss_filename = "loss_full_model.csv" if not null_dataset else "loss_null_model.csv"
    with open(f"{output_path}/{pert}_{loss_filename}", mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Iteration", "Loss"])
        for idx, loss_value in enumerate(LossValues):
            writer.writerow([idx + 1, loss_value])

    # Save training info
    if null_dataset == False:
        output_file_path = f"{output_path}/{pert}_training_info_df.csv"
        training_info_df.to_csv(output_file_path, index=False)

    # Save model states
    list_state_dict_full = []
    for model, likelihood, _ in model_list:
        list_state_dict_full.append(model.state_dict())

    # Create results dictionary
    train_results_dict = {
        "full_model_list": model_list,  # List of (model, likelihood, metadata) tuples
        "full_state_dict_list": list_state_dict_full,
        "full_fit_parameters_df": training_info_df,
        "full_mll_values": mll_values_df,
        "modeled_proteins_list": proteins2test,
        "full_model_metadata": metadata_list,
        "batch_size": batch_size,
        "n_batched_models": n_models,
        "n_total_tasks": len(metadata_list),
        "proteins": proteins2test,
        "n_proteins": n_prot,
        "conditions": conds,
        "n_conditions": n_cond,
        "exactgp_input": tpptr_df,
        "output_path": output_path,
        "device": str(device),
    }

    # Update checklist
    if is_notebook():
        if null_dataset == False:
            tasks[0][2][6] = ("Save models and training results", True)
            tasks[0] = ("1. Build and fit full model", True, tasks[0][2])
            update_checklist(tasks)

    if null_dataset == False:
        # Save updates to summery file
        with open(f"{output_path}/{pert}_Thermal_Tracks_summary.txt", "a") as f:
            f.write("\n" + "=" * 70 + "\n")
            f.write("TRAINING SUMMARY\n")
            f.write("=" * 70 + "\n")
            f.write(f"Total tasks trained: {len(metadata_list)}\n")
            f.write(f"Batched into: {n_models} models\n")
            f.write(f"Batch size: {batch_size}\n")
            f.write(f"Device: {device}\n")
            f.write(f"Parallel processing: {use_parallel}\n")
            if use_parallel:
                f.write(f"Number of workers: {n_workers}\n")
            f.write(f"Model building time: {build_time:.2f} seconds\n")
            f.write(f"Training time: {train_time:.2f} seconds\n")
            f.write(f"MLL computation time: {mll_time:.2f} seconds\n")
            f.write(f"Total time: {build_time + train_time + mll_time:.2f} seconds\n")
            if LossValues:
                f.write(f"Final loss: {LossValues[-1]:.6f}\n")
            f.write(f"Iterations: {len(LossValues)}\n")
            f.write("=" * 70 + "\n")

    else:
        if is_notebook():
            tasks[3][2][0] = ("Train null model", True)
            update_checklist(tasks)
        else:
            pass

    return train_results_dict
