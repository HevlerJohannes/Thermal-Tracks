import torch
import gpytorch
import numpy as np
import pandas as pd
import matplotlib
import psutil
from tqdm import tqdm
import os
from concurrent.futures import ThreadPoolExecutor
import pickle
from functools import partial
from typing import List, Tuple, Dict, Optional
import warnings

warnings.filterwarnings("ignore")

from utils import is_notebook

if is_notebook():
    # Notebook settings
    import matplotlib.pyplot as plt
else:
    # Terminal settings
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt


# Sigmoid decreasing mean function for melting curve prior
class SigmoidDecreasingMean(gpytorch.means.Mean):
    """
    Sigmoid decreasing mean function: m(x) = 1 / (1 + exp(rate * (x - midpoint))).

    Models a typical melting curve prior where intensity is ~1 at low temperatures
    and decays to ~0 at high temperatures. Both midpoint and rate are learned
    during training (rate is constrained positive via softplus).

    Args:
        batch_shape: Batch dimensions for independent GP tasks.
    """

    def __init__(self, batch_shape=torch.Size([])):
        super().__init__()
        self.register_parameter(
            "raw_midpoint", torch.nn.Parameter(torch.full((*batch_shape, 1), 55.0))
        )
        self.register_parameter(
            "raw_rate", torch.nn.Parameter(torch.full((*batch_shape, 1), 0.15))
        )

    def forward(self, x):
        midpoint = self.raw_midpoint  # [*batch, 1]
        rate = torch.nn.functional.softplus(self.raw_rate)  # ensure positive
        # x may be [n, d] — squeeze feature dim to get [n]
        if x.dim() >= 2:
            x = x.squeeze(-1)
        # For batched params [batch, 1], expand x to [1, n] for broadcasting
        if midpoint.dim() > x.dim():
            x = x.unsqueeze(0)
        return 1.0 / (1.0 + torch.exp(rate * (x - midpoint)))


def create_mean_function(mean_spec, batch_shape=torch.Size([])):
    """
    Factory for GP mean functions.

    Args:
        mean_spec: One of "sigmoid", "constant", "zero" (str), None (defaults to
            ConstantMean), or a pre-instantiated gpytorch.means.Mean object.
        batch_shape: Batch dimensions passed to the mean constructor.

    Returns:
        A gpytorch.means.Mean instance.
    """
    if isinstance(mean_spec, str):
        if mean_spec == "sigmoid":
            return SigmoidDecreasingMean(batch_shape=batch_shape)
        elif mean_spec == "constant":
            return gpytorch.means.ConstantMean(batch_shape=batch_shape)
        elif mean_spec == "zero":
            return gpytorch.means.ZeroMean(batch_shape=batch_shape)
        else:
            raise ValueError(f"Unknown mean function: {mean_spec}")
    elif mean_spec is None:
        return gpytorch.means.ConstantMean(batch_shape=batch_shape)
    else:
        return mean_spec


# Batch-independent GP model with batching support for INDEPENDENT tasks
class BatchIndependentMultitaskGPModel(gpytorch.models.ExactGP):
    """
    Batch-independent GP model that stacks multiple protein-condition pairs into
    a single model via batch_shape for vectorized training and prediction.

    Args:
        train_x: Training inputs, shape [n_points] (broadcast across batch).
        train_y: Training targets, shape [num_tasks, n_points].
        likelihood: GaussianLikelihood with matching batch_shape.
        scale_kernel: If True, wraps RBFKernel in a ScaleKernel.
        mean: Mean function specification (str, None, or gpytorch Mean object).
        num_tasks: Number of independent GP tasks in the batch.
        lengthscale_prior: Prior on the RBF lengthscale parameter.
        lengthscale_constraint: Pre-computed constraint on the lengthscale.
        lengthscale_mult: Multiplier applied to the lengthscale constraint.
        outputscale_prior: Prior on the ScaleKernel outputscale (only used
            when scale_kernel=True).
    """

    def __init__(
        self,
        train_x,
        train_y,
        likelihood,
        scale_kernel,
        mean=None,
        num_tasks=1,
        lengthscale_prior=None,
        lengthscale_constraint=None,
        lengthscale_mult=1.0,
        outputscale_prior=None,
    ):
        super(BatchIndependentMultitaskGPModel, self).__init__(
            train_x, train_y, likelihood
        )

        # Set batch_shape based on num_tasks
        batch_shape = torch.Size([num_tasks])

        # Mean module with batch support
        self.mean_module = create_mean_function(mean, batch_shape=batch_shape)

        # Covariance module with batch support
        if scale_kernel:
            self.covar_module = gpytorch.kernels.ScaleKernel(
                gpytorch.kernels.RBFKernel(
                    lengthscale_prior=lengthscale_prior,
                    lengthscale_constraint=lengthscale_constraint,
                    batch_shape=batch_shape,
                ),
                batch_shape=batch_shape,
                outputscale_prior=outputscale_prior,
            )
        else:
            self.covar_module = gpytorch.kernels.RBFKernel(
                lengthscale_prior=lengthscale_prior,
                lengthscale_constraint=lengthscale_constraint,
                batch_shape=batch_shape,
            )

        self.num_tasks = num_tasks

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)


def calculate_optimal_batch_size(
    parameters, n_datapoints, device="auto", null_dataset=[True, False]
):
    """
    Calculate optimal batch size based on available resources

    Args:
        parameters: Dictionary containing:
            - 'cpu_memory_fraction': Fraction of RAM to use (default: 0.5)
            - 'gpu_memory_fraction': Fraction of GPU to use (default: 0.8)
            - 'cpu_thread_strategy': Threading strategy (default: 'all_physical')
                Options: 'all_physical', 'all_logical', 'half_physical',
                        'leave_one', 'leave_two', or integer
            - 'max_batch_size_gpu': Maximum batch size for GPU (default: 64)
            - 'max_batch_size_cpu': Maximum batch size for CPU (default: 16)
            - 'min_batch_size': Minimum batch size (default: 1)
        n_datapoints: Number of data points per model
        device: 'auto', 'cpu', or 'cuda'

    Returns:
        batch_size: Optimal number of tasks per batch
        use_gpu: Whether to use GPU
        torch_device: PyTorch device object
    """
    # Get configuration parameters with sensible defaults
    cpu_memory_fraction = parameters.get("cpu_memory_fraction", 0.5)
    gpu_memory_fraction = parameters.get("gpu_memory_fraction", 0.8)
    thread_strategy = parameters.get("cpu_thread_strategy", "all_physical")
    max_batch_gpu = parameters.get("max_batch_size_gpu", 64)  # Default: 64
    max_batch_cpu = parameters.get("max_batch_size_cpu", 16)  # Default: 16
    min_batch = parameters.get("min_batch_size", 2)

    # Memory per task (in GB): n² × 8 bytes × 5 (overhead factor)
    gb_per_task = (n_datapoints**2 * 8 * 5) / 1e9

    # Determine device
    if device == "auto":
        use_gpu = torch.cuda.is_available()
        device_type = "cuda" if use_gpu else "cpu"
    elif device == "cuda":
        use_gpu = torch.cuda.is_available()
        device_type = "cuda" if use_gpu else "cpu"
        if not use_gpu:
            print("CUDA requested but not available, falling back to CPU")
    elif device == "cpu":
        use_gpu = False
        device_type = "cpu"
        print("CPU mode explicitly selected")
    else:
        print(f"Unknown device '{device}', defaulting to auto")
        use_gpu = torch.cuda.is_available()
        device_type = "cuda" if use_gpu else "cpu"

    # Set torch device
    torch_device = torch.device(device_type)

    # Calculate available memory and configure threads
    if use_gpu:
        gpu_memory_gb = torch.cuda.get_device_properties(0).total_memory / 1e9
        available_memory = gpu_memory_gb * gpu_memory_fraction
        device_name = torch.cuda.get_device_properties(0).name
        print(f"Using GPU: {device_name}")
        print(f"Total GPU memory: {gpu_memory_gb*1000:.1f} GB")
        print(f"Using {gpu_memory_fraction*100:.0f}%: {available_memory:.100f} GB")
        n_threads = None  # Not applicable for GPU
    else:
        # RAM calculation
        total_ram_gb = psutil.virtual_memory().total / 1e9
        available_memory = total_ram_gb * cpu_memory_fraction
        device_name = "CPU"

        # CPU thread configuration
        physical_cores = psutil.cpu_count(logical=False)
        logical_cores = psutil.cpu_count(logical=True)

        if thread_strategy == "all_physical":
            n_threads = physical_cores
        elif thread_strategy == "all_logical":
            n_threads = logical_cores
        elif thread_strategy == "half_physical":
            n_threads = max(1, physical_cores // 2)
        elif thread_strategy == "leave_one":
            n_threads = max(1, physical_cores - 1)
        elif thread_strategy == "leave_two":
            n_threads = max(1, physical_cores - 2)
        elif isinstance(thread_strategy, int):
            n_threads = max(1, min(thread_strategy, logical_cores))
        else:
            print(
                f"Unknown thread strategy '{thread_strategy}', using all physical cores"
            )
            n_threads = physical_cores

        # Set thread limits
        torch.set_num_threads(n_threads)
        os.environ["OMP_NUM_THREADS"] = str(n_threads)
        os.environ["MKL_NUM_THREADS"] = str(n_threads)
        os.environ["OPENBLAS_NUM_THREADS"] = str(n_threads)

        print(f"Using CPU with {n_threads} threads")
        print(f"Physical cores: {physical_cores}, Logical cores: {logical_cores}")
        print(f"Thread strategy: {thread_strategy}")
        print(f"Total RAM: {total_ram_gb:.1f} GB")
        print(f"Using {cpu_memory_fraction*100:.0f}%: {available_memory:.10f} GB")
        print(f"Reserved for system: {total_ram_gb - available_memory:.10f} GB")

    # Calculate batch size based on memory
    max_batch_size = int(available_memory / gb_per_task)

    # Round down to power of 2 for efficiency
    batch_size_power2 = 2 ** int(np.log2(max_batch_size)) if max_batch_size > 0 else 1

    # Apply user-defined limits (with sensible defaults)
    if use_gpu:
        if max_batch_gpu is not None:
            batch_size = min(batch_size_power2, max_batch_gpu)
            limit_applied = max_batch_gpu if batch_size_power2 > max_batch_gpu else None
        else:
            batch_size = batch_size_power2
            limit_applied = None
    else:
        if max_batch_cpu is not None:
            batch_size = min(batch_size_power2, max_batch_cpu)
            limit_applied = max_batch_cpu if batch_size_power2 > max_batch_cpu else None
        else:
            batch_size = batch_size_power2
            limit_applied = None

    # Apply minimum
    batch_size = max(min_batch, batch_size)

    # Get Model information
    pert = parameters.get("perturbation", "unknown")
    output_path = str(parameters.get("result_dir", "./"))

    # Write summary to file
    if null_dataset == False:
        with open(f"{output_path}/{pert}_Thermal_Tracks_summary.txt", "w") as f:
            f.write("\n" + "=" * 70 + "\n")
            f.write("PARAMETERS\n")
            f.write("=" * 70 + "\n")
            for key, value in parameters.items():
                f.write(f"  {key}: {value}\n")
            f.write("\n" + "=" * 70 + "\n")
            f.write("RESOURCE CALCULATION\n")
            f.write("=" * 70 + "\n")
            f.write(f"  Device: {device_name}\n")
            f.write(f"  Data points per model: {n_datapoints}\n")
            f.write(f"  Memory per task: {gb_per_task:.4f} GB\n")

            if use_gpu:
                f.write(f"Total GPU memory: {gpu_memory_gb:.2f} GB\n")
                f.write(f"GPU memory fraction: {gpu_memory_fraction*1000:.0f}%\n")
                f.write(f"Available memory: {available_memory:.2f} GB\n")
                f.write(
                    f"Max batch size limit: {max_batch_gpu if max_batch_gpu else 'None'}\n"
                )
            else:
                f.write(f"Total RAM: {total_ram_gb:.2f} GB\n")
                f.write(f"CPU memory fraction: {cpu_memory_fraction*1000:.0f}%\n")
                f.write(f"Available memory: {available_memory:.2f} GB\n")
                f.write(
                    f"Reserved for system: {total_ram_gb - available_memory:.2f} GB\n"
                )
                f.write(f"Physical cores: {physical_cores}\n")
                f.write(f"Logical cores: {logical_cores}\n")
                f.write(f"Thread strategy: {thread_strategy}\n")
                f.write(f"CPU threads: {n_threads}\n")
                f.write(
                    f"Max batch size limit: {max_batch_cpu if max_batch_cpu else 'None'}\n"
                )

            f.write(f"Min batch size: {min_batch}\n")
            f.write(f"Memory-based batch size: {max_batch_size}\n")
            f.write(f"Power-of-2 batch size: {batch_size_power2}\n")
            f.write(f"Final batch size: {batch_size}\n")
            f.write(f"Memory per batch: {gb_per_task * batch_size:.2f} GB\n")
            f.write("=" * 70 + "\n\n")

    return batch_size, use_gpu, torch_device


def precompute_lengthscale_constraint(
    train_x_values, lengthscale_minconstraint, lengthscale_mult
):
    """
    Pre-compute a GreaterThan lengthscale constraint from the spacing of temperature points.

    The constraint prevents the GP from fitting noise at scales smaller than the
    data resolution. Computed once and shared across all models in a batch.

    Args:
        train_x_values: Sorted temperature tensor for one signature group.
        lengthscale_minconstraint: Statistic to use on consecutive differences
            ("min", "mean", "median", "max") or None to skip.
        lengthscale_mult: Multiplier applied to the chosen statistic.

    Returns:
        gpytorch.constraints.GreaterThan or None.
    """
    if lengthscale_minconstraint is None:
        return None

    if len(train_x_values) <= 1:
        return None

    DistVec = train_x_values[1:] - train_x_values[:-1]

    if lengthscale_minconstraint == "min":
        Constt = lengthscale_mult * torch.min(DistVec)
    elif lengthscale_minconstraint == "mean":
        Constt = lengthscale_mult * torch.mean(DistVec)
    elif lengthscale_minconstraint == "median":
        Constt = lengthscale_mult * torch.median(DistVec)
    elif lengthscale_minconstraint == "max":
        Constt = lengthscale_mult * torch.max(DistVec)
    else:
        return None

    return gpytorch.constraints.GreaterThan(Constt)


def create_single_batch_model(batch_data, params):
    """
    Create a single BatchIndependentMultitaskGPModel for one chunk of protein-condition pairs.

    Called by building_exactgp_model_batched_optimized, either sequentially or via
    ThreadPoolExecutor. All configuration is packed into batch_data to support
    parallel dispatch.

    Args:
        batch_data: Tuple of (chunk_idx, chunk_pairs, train_x_template, device,
            scale_kernel, lengthscale_prior, lengthscale_constraint,
            lengthscale_mult, start_idx, dtype, mean_spec, noise_prior,
            noise_constraint, outputscale_prior).
        params: Unused (kept for interface compatibility with ThreadPoolExecutor).

    Returns:
        Tuple of (model, likelihood, metadata_list) for this chunk.
    """
    (
        chunk_idx,
        chunk_pairs,
        train_x_template,
        device,
        scale_kernel,
        lengthscale_prior,
        lengthscale_constraint,
        lengthscale_mult,
        start_idx,
        dtype,
        mean_spec,
        noise_prior,
        noise_constraint,
        outputscale_prior,
    ) = batch_data

    current_batch_size = len(chunk_pairs)

    # Use pre-allocated template x values with correct dtype
    train_x = train_x_template.clone().to(dtype=dtype)

    # Stack y values efficiently with correct dtype
    train_y = torch.stack(
        [pair["intens"].to(dtype=dtype) for pair in chunk_pairs], dim=0
    )

    # Create likelihood with noise prior for outlier robustness
    likelihood = gpytorch.likelihoods.GaussianLikelihood(
        batch_shape=torch.Size([current_batch_size]),
        noise_prior=noise_prior,
        noise_constraint=noise_constraint,
    ).to(device)

    # Create model with pre-computed constraint
    model = BatchIndependentMultitaskGPModel(
        train_x,
        train_y,
        likelihood,
        scale_kernel=scale_kernel,
        mean=mean_spec,
        num_tasks=current_batch_size,
        lengthscale_prior=lengthscale_prior,
        lengthscale_constraint=lengthscale_constraint,
        lengthscale_mult=lengthscale_mult,
        outputscale_prior=outputscale_prior,
    ).to(device)

    # Create metadata
    metadata = []
    for local_idx, pair in enumerate(chunk_pairs):
        global_idx = start_idx + local_idx
        metadata.append(
            {
                "protein": pair["protein"],
                "condition": pair["condition"],
                "chunk_index": chunk_idx,
                "task_index": local_idx,
                "global_index": global_idx,
            }
        )

    return (model, likelihood, metadata)


def building_exactgp_model_batched_optimized(
    tpptr_df,
    parameters,
    proteins2test,
    conds,
    scale_kernel,
    lengthscale_prior,
    lengthscale_minconstraint,
    lengthscale_mult,
    mean,
    batch_size=8,
    device=None,
    null_dataset=[True, False],
    use_parallel=True,
    n_workers=4,
    dtype=torch.float64,
):
    """
    Build batched GP models for all protein-condition pairs.

    Pairs are grouped by temperature signature (identical x-values) because
    batching requires uniform tensor dimensions. Each group is split into chunks
    of at most batch_size, and a BatchIndependentMultitaskGPModel is created per
    chunk. Model creation can optionally be parallelized across chunks.

    Args:
        tpptr_df: Input dataframe with columns 'uniqueID', 'condition', 'x', 'y'.
        parameters: Full parameter dictionary (used to extract MeanFunction,
            noise_prior, noise_constraint, outputscale_prior, and output paths).
        proteins2test: Array of protein IDs to include.
        conds: Array of condition labels.
        scale_kernel: Whether to wrap RBFKernel in ScaleKernel.
        lengthscale_prior: Prior on lengthscale (or None).
        lengthscale_minconstraint: Statistic for lengthscale lower bound.
        lengthscale_mult: Multiplier on the lengthscale constraint.
        mean: Mean function specification (forwarded to create_mean_function).
        batch_size: Maximum number of tasks per batched model.
        device: torch.device for model placement.
        null_dataset: Whether this is building null (sampled) models.
        use_parallel: Enable parallel chunk creation via ThreadPoolExecutor.
        n_workers: Number of worker threads for parallel creation.
        dtype: torch.float32 or torch.float64.

    Returns:
        model_list: List of (model, likelihood, metadata) tuples.
        n_cond: Number of conditions.
        n_prot: Number of proteins.
        n_models: Number of batched models created.
    """
    if device is None:
        device = torch.device("cpu")

    pert = parameters.get("perturbation", "unknown")
    output_path = str(parameters.get("result_dir", "./"))

    # Write header
    mode_text = "NULL GP MODELS" if null_dataset else "FULL GP MODELS"
    with open(f"{output_path}/{pert}_Thermal_Tracks_summary.txt", "a") as f:
        f.write("\n" + "=" * 70 + "\n")
        f.write(f"BUILDING {mode_text}\n")
        f.write("=" * 70 + "\n")

    # ========== Vectorize data preparation ==========
    # Determine numpy dtype based on specified torch dtype
    if dtype == torch.float32:
        np_dtype = np.float32
    elif dtype == torch.float64:
        np_dtype = np.float64
    else:
        np_dtype = np.float32

    # Convert to specified dtype once
    tpptr_df = tpptr_df.copy()
    tpptr_df["x"] = tpptr_df["x"].astype(np_dtype)
    tpptr_df["y"] = tpptr_df["y"].astype(np_dtype)

    # Group by protein and condition efficiently
    grouped = tpptr_df[tpptr_df["uniqueID"].isin(proteins2test)].groupby(
        ["uniqueID", "condition"]
    )

    # Pre-allocate list for all pairs
    all_pairs = []

    # Process groups in vectorized manner
    for (prot, cond), group in grouped:
        if cond in conds:
            # Convert to tensors with specified dtype and move to device
            temp_tensor = torch.from_numpy(group["x"].values).to(
                dtype=dtype, device=device
            )
            intens_tensor = torch.from_numpy(group["y"].values).to(
                dtype=dtype, device=device
            )

            all_pairs.append(
                {
                    "protein": prot,
                    "condition": cond,
                    "temp": temp_tensor,
                    "intens": intens_tensor,
                }
            )

    n_total_tasks = len(all_pairs)

    if not all_pairs:
        print("No data pairs found!")
        return [], 0, 0, 0

    # ========== Group pairs by temperature signature ==========
    # Proteins may have different numbers of observed temperatures.
    # Batching requires uniform tensor sizes, so we group pairs by their
    # exact set of temperature points — only same-signature pairs are batched.
    temp_groups = {}
    for pair in all_pairs:
        temp_key = tuple(pair["temp"].cpu().tolist())
        if temp_key not in temp_groups:
            temp_groups[temp_key] = []
        temp_groups[temp_key].append(pair)

    n_temp_groups = len(temp_groups)
    group_summary = {
        len(k): sum(1 for tk in temp_groups if len(tk) == len(k)) for k in temp_groups
    }

    # Write stats
    with open(f"{output_path}/{pert}_Thermal_Tracks_summary.txt", "a") as f:
        f.write(f"Total protein-condition pairs: {n_total_tasks}\n")
        f.write(f"Unique temperature signatures: {n_temp_groups}\n")
        for temp_key, grp in temp_groups.items():
            f.write(f"  {len(temp_key)} temps ({len(grp)} pairs)\n")
        f.write(f"Batch size: {batch_size}\n")
        f.write(f"Device: {device}\n")
        f.write(f"Parallel processing: {use_parallel}\n")
        if use_parallel:
            f.write(f"Number of workers: {n_workers}\n")

    # ========== Prepare chunks for parallel processing ==========
    mean_spec = parameters.get("MeanFunction", "constant")
    noise_prior = parameters.get("noise_prior", None)
    noise_constraint = parameters.get("noise_constraint", None)
    outputscale_prior = parameters.get("outputscale_prior", None)

    chunks_data = []
    chunk_idx = 0
    global_idx = 0

    for temp_key, group_pairs in temp_groups.items():
        # Each temperature group gets its own x-template and lengthscale constraint
        group_x_template = group_pairs[0]["temp"]
        group_ls_constraint = precompute_lengthscale_constraint(
            group_x_template, lengthscale_minconstraint, lengthscale_mult
        )

        n_group_chunks = int(np.ceil(len(group_pairs) / batch_size))
        for i in range(n_group_chunks):
            start = i * batch_size
            end = min(start + batch_size, len(group_pairs))
            chunk_pairs = group_pairs[start:end]

            chunks_data.append(
                (
                    chunk_idx,
                    chunk_pairs,
                    group_x_template,
                    device,
                    scale_kernel,
                    lengthscale_prior,
                    group_ls_constraint,
                    lengthscale_mult,
                    global_idx,
                    dtype,
                    mean_spec,
                    noise_prior,
                    noise_constraint,
                    outputscale_prior,
                )
            )

            chunk_idx += 1
            global_idx += len(chunk_pairs)

    n_chunks = len(chunks_data)

    # ========== Parallel or sequential model creation ==========
    if use_parallel and n_chunks > 1:

        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            model_list = list(
                tqdm(
                    executor.map(
                        lambda x: create_single_batch_model(x, None), chunks_data
                    ),
                    total=n_chunks,
                    desc="Creating batched models",
                )
            )
    else:
        model_list = []
        for chunk_data in tqdm(chunks_data, desc="Creating batched models"):
            model_list.append(create_single_batch_model(chunk_data, None))

    n_models = len(model_list)
    n_cond = len(conds)
    n_prot = len(proteins2test)

    # Store dtype information in the first model for later reference
    if model_list:
        for model, likelihood, metadata in model_list:
            model.dtype_used = dtype  # Store dtype for prediction phase

    with open(f"{output_path}/{pert}_Thermal_Tracks_summary.txt", "a") as f:
        f.write(
            f"Created {n_models} batched models containing {n_total_tasks} total tasks\n"
        )
        f.write(f"Data type: {dtype}\n")
        f.write(f"=" * 70 + "\n")

    return model_list, n_cond, n_prot, n_models
