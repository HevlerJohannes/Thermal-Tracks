#! /usr/bin/env python3
import numpy as np
import gpytorch
import torch
from matplotlib import pyplot as plt

####################################################################
#                      Multitask Exact GP                          #
####################################################################

# define a GP model using exact inference with constant mean function and squared exponential (RBF) kerclass ExactGPModel(gpytorch.models.ExactGP):
class ExactGPModel(gpytorch.models.ExactGP):
        def __init__(self, train_x, train_y, likelihood, mean = gpytorch.means.ZeroMean(), lengthscale_prior = None, lengthscale_minconstraint = ['min', 'mean', 'median', 'max', None], lengthscale_mult: float = 1.0):
            super(ExactGPModel, self).__init__(train_x, train_y, likelihood)

            self.mean_module = mean
            
            if not lengthscale_minconstraint is None:
                train_x_values = train_x.unique()
                DistVec = train_x_values[1:len(train_x_values)]-train_x_values[0:len(train_x_values)-1] # compute the distance between each pair of consecutive temperatures
                
                if lengthscale_minconstraint == 'min':
                    Constt = lengthscale_mult * min(DistVec) # we use the minium distance between temperatures as a lower limit for the lengthscale
                elif lengthscale_minconstraint == 'mean':
                    Constt = lengthscale_mult * torch.mean(DistVec) # we use the mean distance between temperatures as a lower limit for the lengthscale
                elif lengthscale_minconstraint == 'median':
                    Constt = lengthscale_mult * torch.median(DistVec) # we use the median distance between temperatures as a lower limit for the lengthscale
                elif lengthscale_minconstraint == 'max':
                    Constt = lengthscale_mult * max(DistVec) # we use the max distance between temperature as a lower limit for the lengthscale
                
                lengthscale_constraint = gpytorch.constraints.GreaterThan(Constt)
            else:
                lengthscale_constraint = None # there is no lower limit for the lengthscale
        
        
            self.covar_module = gpytorch.kernels.RBFKernel(lengthscale_prior = lengthscale_prior, lengthscale_constraint = lengthscale_constraint)
                
    
        def forward(self, x):
            mean_x = self.mean_module(x)
            covar_x = self.covar_module(x)
            return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)
    
# Model definition
def building_exactgp_model(tpptr_df, proteins2test, conds, lengthscale_prior, lengthscale_minconstraint, lengthscale_mult, mean = gpytorch.means.ZeroMean()):
    """
    Defines a GP model using exact inference with zero mean function and squared exponential (RBF) kernel class ExactGPModel(gpytorch.models.ExactGP)
    Args:
        - tpptr_df (pandas.core.frame.DataFrame):
            Expected is a pandas data frame with variables
            'uniqueID' for protein ids,
            'x' for profiled temperatures,
            'y' for scaled protein intensities,
            'condition' for performed comparison (e.g. vehicle vs drug)
            
            Note that the input data needs to be already filtered and formatted (e.g. proteins that should be analyzed need to be identified in both conditions).
.
        - conds (numpy.ndarray): array containing conditions that melting curves will be modeled for.
        - lengthscale_prior (str): Set this if you want to apply a prior to the lengthscale parameter. (Default: None). See: https://docs.gpytorch.ai/en/stable/priors.html for possible priors.
        - lengthscale_minconstraint (str): Set this if you want to apply a constraint to the lengthscale parameter. Possible options are: 'min', 'mean', 'median', 'max', None.
            
            This constraint will limit the lowest possible value for the lengthscale (implemented via gpytorch.constraints.GreaterThan) and was adopted from:
                LeSueur, Cecile, Magnus Rattray, and Mikhail Savitski. 2023. “Hierarchical Gaussian Process Models Explore the Dark Meltome of Thermal Proteome Profiling Experiments.” bioRxiv. https://doi.org/10.1101/2023.10.26.564129.
        
        - mean (mean function): Defines the mean function of the GP. (Default: gpytorch.means.ZeroMean()). See: https://docs.gpytorch.ai/en/stable/means.html for alternatives.
            
    Returns:
        gpytorch.models.model_list.IndependentModelList: ModelList of independent GP models. 
    For more details, see the GPyTorch documentation:
    https://docs.gpytorch.ai/en/stable/examples/03_Multitask_Exact_GPs/ModelList_GP_Regression.html
    """
    model_list = []
    l_list = []
    metadata_list = []  # List to store metadata
    for prot in proteins2test:
        # filter data
        df = tpptr_df[tpptr_df['uniqueID'] == prot]
        # define alternative models (one model for each condition)
        for cond in df['condition'].unique():
            dfcond = df[df['condition'] == cond]
            temp = torch.as_tensor(np.asarray(dfcond['x'])).double()
            intens = torch.as_tensor(np.asarray(dfcond['y'])).double()
            lik = gpytorch.likelihoods.GaussianLikelihood()
            l_list.append(lik)
            model = ExactGPModel(temp, intens, lik, mean, lengthscale_prior, lengthscale_minconstraint, lengthscale_mult)
            model_list.append(model)
            model_id = len(model_list)
            metadata_list.append((prot, cond, model_id))
    # dimensions
    n_models = len(model_list)
    n_cond = len(conds)
    n_prot = len(proteins2test)
    assert n_models == n_cond * n_prot, "Error in building the model list."
    
    # create ModelList Multioutput GP
    return model_list, l_list, metadata_list, n_cond, n_prot, n_models 
