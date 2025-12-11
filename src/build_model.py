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
warnings.filterwarnings('ignore')

from utils import is_notebook

if is_notebook():
    # Notebook settings
    import matplotlib.pyplot as plt
else:
    # Terminal settings
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

# Batch-independent GP model with batching support for INDEPENDENT tasks
class BatchIndependentMultitaskGPModel(gpytorch.models.ExactGP):
    """
    Batch-independent GP model with batching support for INDEPENDENT tasks
    """
    def __init__(self, train_x, train_y, likelihood, scale_kernel, mean=None, 
                 num_tasks=1,
                 lengthscale_prior=None, 
                 lengthscale_constraint=None,  # Changed to accept pre-computed constraint
                 lengthscale_mult=1.0):
        super(BatchIndependentMultitaskGPModel, self).__init__(train_x, train_y, likelihood)
        
        # Set batch_shape based on num_tasks
        batch_shape = torch.Size([num_tasks]) if num_tasks > 1 else torch.Size([])
        
        # Mean module with batch support
        if num_tasks > 1:
            self.mean_module = gpytorch.means.ConstantMean(batch_shape=batch_shape)
        else:
            self.mean_module = mean if mean is not None else gpytorch.means.ConstantMean()
        
        # Covariance module with batch support
        if scale_kernel:
            if num_tasks > 1:
                self.covar_module = gpytorch.kernels.ScaleKernel(
                    gpytorch.kernels.RBFKernel(
                        lengthscale_prior=lengthscale_prior,
                        lengthscale_constraint=lengthscale_constraint,
                        batch_shape=batch_shape
                    ),
                    batch_shape=batch_shape
                )
            else:
                self.covar_module = gpytorch.kernels.ScaleKernel(
                    gpytorch.kernels.RBFKernel(
                        lengthscale_prior=lengthscale_prior,
                        lengthscale_constraint=lengthscale_constraint
                    )
                )
        else:
            if num_tasks > 1:
                self.covar_module = gpytorch.kernels.RBFKernel(
                    lengthscale_prior=lengthscale_prior,
                    lengthscale_constraint=lengthscale_constraint,
                    batch_shape=batch_shape
                )
            else:
                self.covar_module = gpytorch.kernels.RBFKernel(
                    lengthscale_prior=lengthscale_prior,
                    lengthscale_constraint=lengthscale_constraint
                )
        
        self.num_tasks = num_tasks
    
    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)
    

def calculate_optimal_batch_size(parameters, n_datapoints, device='auto', null_dataset=[True,False]):
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
    cpu_memory_fraction = parameters.get('cpu_memory_fraction', 0.5)
    gpu_memory_fraction = parameters.get('gpu_memory_fraction', 0.8)
    thread_strategy = parameters.get('cpu_thread_strategy', 'all_physical')
    max_batch_gpu = parameters.get('max_batch_size_gpu', 64)  # Default: 64
    max_batch_cpu = parameters.get('max_batch_size_cpu', 16)  # Default: 16
    min_batch = parameters.get('min_batch_size', 2)
    
    # Memory per task (in GB): n² × 8 bytes × 5 (overhead factor)
    gb_per_task = (n_datapoints**2 * 8 * 5) / 1e9
    
    # Determine device
    if device == 'auto':
        use_gpu = torch.cuda.is_available()
        device_type = 'cuda' if use_gpu else 'cpu'
    elif device == 'cuda':
        use_gpu = torch.cuda.is_available()
        device_type = 'cuda' if use_gpu else 'cpu'
        if not use_gpu:
            print("CUDA requested but not available, falling back to CPU")
    elif device == 'cpu':
        use_gpu = False
        device_type = 'cpu'
        print("CPU mode explicitly selected")
    else:
        print(f"Unknown device '{device}', defaulting to auto")
        use_gpu = torch.cuda.is_available()
        device_type = 'cuda' if use_gpu else 'cpu'
    
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
        
        if thread_strategy == 'all_physical':
            n_threads = physical_cores
        elif thread_strategy == 'all_logical':
            n_threads = logical_cores
        elif thread_strategy == 'half_physical':
            n_threads = max(1, physical_cores // 2)
        elif thread_strategy == 'leave_one':
            n_threads = max(1, physical_cores - 1)
        elif thread_strategy == 'leave_two':
            n_threads = max(1, physical_cores - 2)
        elif isinstance(thread_strategy, int):
            n_threads = max(1, min(thread_strategy, logical_cores))
        else:
            print(f"Unknown thread strategy '{thread_strategy}', using all physical cores")
            n_threads = physical_cores
        
        # Set thread limits
        torch.set_num_threads(n_threads)
        os.environ['OMP_NUM_THREADS'] = str(n_threads)
        os.environ['MKL_NUM_THREADS'] = str(n_threads)
        os.environ['OPENBLAS_NUM_THREADS'] = str(n_threads)
        
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
    pert = parameters.get('perturbation', 'unknown')
    output_path = str(parameters.get('result_dir', './'))
    
    # Write summary to file
    if null_dataset== False:
        with open(f'{output_path}/{pert}_Thermal_Tracks_summary.txt', 'w') as f:
            f.write("\n" + "="*70 + "\n")
            f.write("PARAMETERS\n")
            f.write("="*70 + "\n")
            for key, value in parameters.items():
                f.write(f"  {key}: {value}\n")
            f.write("\n" + "="*70 + "\n")
            f.write("RESOURCE CALCULATION\n")
            f.write("="*70 + "\n")
            f.write(f"  Device: {device_name}\n")
            f.write(f"  Data points per model: {n_datapoints}\n")
            f.write(f"  Memory per task: {gb_per_task:.4f} GB\n")
            
            if use_gpu:
                f.write(f"Total GPU memory: {gpu_memory_gb:.2f} GB\n")
                f.write(f"GPU memory fraction: {gpu_memory_fraction*1000:.0f}%\n")
                f.write(f"Available memory: {available_memory:.2f} GB\n")
                f.write(f"Max batch size limit: {max_batch_gpu if max_batch_gpu else 'None'}\n")
            else:
                f.write(f"Total RAM: {total_ram_gb:.2f} GB\n")
                f.write(f"CPU memory fraction: {cpu_memory_fraction*1000:.0f}%\n")
                f.write(f"Available memory: {available_memory:.2f} GB\n")
                f.write(f"Reserved for system: {total_ram_gb - available_memory:.2f} GB\n")
                f.write(f"Physical cores: {physical_cores}\n")
                f.write(f"Logical cores: {logical_cores}\n")
                f.write(f"Thread strategy: {thread_strategy}\n")
                f.write(f"CPU threads: {n_threads}\n")
                f.write(f"Max batch size limit: {max_batch_cpu if max_batch_cpu else 'None'}\n")
            
            f.write(f"Min batch size: {min_batch}\n")
            f.write(f"Memory-based batch size: {max_batch_size}\n")
            f.write(f"Power-of-2 batch size: {batch_size_power2}\n")
            f.write(f"Final batch size: {batch_size}\n")
            f.write(f"Memory per batch: {gb_per_task * batch_size:.2f} GB\n")
            f.write("="*70 + "\n\n")
    
    return batch_size, use_gpu, torch_device


def precompute_lengthscale_constraint(train_x_values, lengthscale_minconstraint, lengthscale_mult):
    """Pre-compute lengthscale constraint once for all models to further reduce computational time"""
    if lengthscale_minconstraint is None:
        return None
    
    if len(train_x_values) <= 1:
        return None
    
    DistVec = train_x_values[1:] - train_x_values[:-1]
    
    if lengthscale_minconstraint == 'min':
        Constt = lengthscale_mult * torch.min(DistVec)
    elif lengthscale_minconstraint == 'mean':
        Constt = lengthscale_mult * torch.mean(DistVec)
    elif lengthscale_minconstraint == 'median':
        Constt = lengthscale_mult * torch.median(DistVec)
    elif lengthscale_minconstraint == 'max':
        Constt = lengthscale_mult * torch.max(DistVec)
    else:
        return None
    
    return gpytorch.constraints.GreaterThan(Constt)


def create_single_batch_model(batch_data, params):
    """Create a single batch model for parallel processing"""
    chunk_idx, chunk_pairs, train_x_template, device, scale_kernel, lengthscale_prior, \
    lengthscale_constraint, lengthscale_mult, start_idx, dtype = batch_data
    
    current_batch_size = len(chunk_pairs)
    
    # Use pre-allocated template x values with correct dtype
    train_x = train_x_template.clone().to(dtype=dtype)
    
    # Stack y values efficiently with correct dtype
    train_y = torch.stack([pair['intens'].to(dtype=dtype) for pair in chunk_pairs], dim=0)
    
    # Create likelihood
    likelihood = gpytorch.likelihoods.GaussianLikelihood(
        batch_shape=torch.Size([current_batch_size])
    ).to(device)
    
    # Create model with pre-computed constraint
    model = BatchIndependentMultitaskGPModel(
        train_x, train_y, likelihood,
        scale_kernel=scale_kernel,
        mean=None,
        num_tasks=current_batch_size,
        lengthscale_prior=lengthscale_prior,
        lengthscale_constraint=lengthscale_constraint,
        lengthscale_mult=lengthscale_mult
    ).to(device)
    
    # Create metadata
    metadata = []
    for local_idx, pair in enumerate(chunk_pairs):
        global_idx = start_idx + local_idx
        metadata.append({
            'protein': pair['protein'],
            'condition': pair['condition'],
            'chunk_index': chunk_idx,
            'task_index': local_idx,
            'global_index': global_idx
        })
    
    return (model, likelihood, metadata)


def building_exactgp_model_batched_optimized(
    tpptr_df, parameters, proteins2test, conds, scale_kernel, 
    lengthscale_prior, lengthscale_minconstraint, 
    lengthscale_mult, mean, batch_size=8, device=None, 
    null_dataset=[True,False], use_parallel=True, n_workers=4, 
    dtype=torch.float64):
    """
    Optimized GP model building with parallel processing and pre-computation
    """
    if device is None:
        device = torch.device('cpu')
    
    pert = parameters.get('perturbation', 'unknown')
    output_path = str(parameters.get('result_dir', './'))
    
    # Write header
    mode_text = "NULL GP MODELS" if null_dataset else "FULL GP MODELS"
    with open(f'{output_path}/{pert}_Thermal_Tracks_summary.txt', 'a') as f:
        f.write("\n" + "="*70 + "\n")
        f.write(f"BUILDING {mode_text}\n")
        f.write("="*70 + "\n")
    
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
    tpptr_df['x'] = tpptr_df['x'].astype(np_dtype)
    tpptr_df['y'] = tpptr_df['y'].astype(np_dtype)
    
    # Group by protein and condition efficiently
    grouped = tpptr_df[tpptr_df['uniqueID'].isin(proteins2test)].groupby(['uniqueID', 'condition'])
    
    # Pre-allocate list for all pairs
    all_pairs = []
    
    # Process groups in vectorized manner
    for (prot, cond), group in grouped:
        if cond in conds:
            # Convert to tensors with specified dtype and move to device
            temp_tensor = torch.from_numpy(group['x'].values).to(dtype=dtype, device=device)
            intens_tensor = torch.from_numpy(group['y'].values).to(dtype=dtype, device=device)
            
            all_pairs.append({
                'protein': prot,
                'condition': cond,
                'temp': temp_tensor,
                'intens': intens_tensor
            })
    
    n_total_tasks = len(all_pairs)
    n_chunks = int(np.ceil(n_total_tasks / batch_size))
    
    # ========== Pre-compute common values ==========
    # Assuming all pairs have the same temperature points
    if all_pairs:
        train_x_template = all_pairs[0]['temp']
        
        # Pre-compute lengthscale constraint once
        lengthscale_constraint = precompute_lengthscale_constraint(
            train_x_template, 
            lengthscale_minconstraint, 
            lengthscale_mult
        )
    else:
        print("No data pairs found!")
        return [], 0, 0, 0
    
    # Write stats
    with open(f'{output_path}/{pert}_Thermal_Tracks_summary.txt', 'a') as f:
        f.write(f"Total protein-condition pairs: {n_total_tasks}\n")
        f.write(f"Batch size: {batch_size}\n")
        f.write(f"Number of batched models (chunks): {n_chunks}\n")
        f.write(f"Device: {device}\n")
        f.write(f"Parallel processing: {use_parallel}\n")
        if use_parallel:
            f.write(f"Number of workers: {n_workers}\n")
    
    # ========== Prepare chunks for parallel processing ==========
    chunks_data = []
    for chunk_idx in range(n_chunks):
        start_idx = chunk_idx * batch_size
        end_idx = min(start_idx + batch_size, n_total_tasks)
        chunk_pairs = all_pairs[start_idx:end_idx]
        
        chunks_data.append((
            chunk_idx, chunk_pairs, train_x_template, device, scale_kernel,
            lengthscale_prior, lengthscale_constraint, lengthscale_mult, start_idx, dtype
        ))
    
    # ========== Parallel or sequential model creation ==========
    if use_parallel and n_chunks > 1:

        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            model_list = list(tqdm(
                executor.map(lambda x: create_single_batch_model(x, None), chunks_data),
                total=n_chunks,
                desc="Creating batched models"
            ))
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
    
    with open(f'{output_path}/{pert}_Thermal_Tracks_summary.txt', 'a') as f:
        f.write(f"Created {n_models} batched models containing {n_total_tasks} total tasks\n")
        f.write(f"Data type: {dtype}\n")
        f.write(f"="*70 + "\n")

    return model_list, n_cond, n_prot, n_models
