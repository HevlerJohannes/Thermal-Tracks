# Thermal Tracks

## Description

**Thermal Tracks** is a robust statistical framework designed to improve the analysis of protein melting curves in TPP-TR experiments.  
While originally developed for temperature-dependent abundance profiles, the framework is applicable to any continuous variable (e.g., temperature, time, concentration) that modulates protein intensity across two conditions (e.g., control and perturbation).

<p align="center">
  <img width="500" alt="Image" src=("https://github.com/user-attachments/assets/2c58a1a5-17f4-4ad9-9236-ec89efa4a949") />
</p>

<p align="center"><small><em>Figure 1. Overview of the Thermal Tracks workflow. The framework fits scaled intensity data from TPP experiments (1) using an RBF kernel (prior distribution) and a Multi-output Gaussian Process (MOGP) model to compare control and treatment conditions. The full model, which includes fits (posterior distribution) for both control (blue) and perturbation (green) conditions, is fitted using Type II maximum likelihood estimation (MLE). Estimated parameters are then incorporated into joint models (black dotted lines), which combine control and treatment conditions. These joint models are used to quantify changes in melting behavior in control vs perturbation by approximating the null distribution of the Thermal Tracks test statistic Λ (2). The observed statistics Λcontrol vs. perturbation are compared to this null distribution approximation to compute empirical p-values. Further, an effect size is determined based on the predicted fits from the full model.</em></small></p>


## Abstract

Thermal Proteome Profiling [TPP](https://www.science.org/doi/10.1126/science.1255784?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed) is an advanced scientific technique that investigates protein behavior under different conditions. By subjecting proteins to heat stress and measuring their solubility at various temperatures using mass spectrometry, researchers can create detailed melting profiles for each protein. These profiles provide rich information about a protein's specific context and can reveal how different factors like small molecule interactions, nucleic acids, protein-protein interactions, and post-translational modifications influence protein stability. The method combines principles from cellular thermal shift assays and quantitative mass spectrometry, offering a powerful approach to understanding protein interactions and changes in living systems.
Until recently, the analysis of melting curves was constrained by limitations in statistical methods, primarily relying on sigmoidal melting curve models that assumed most observations came from a null distribution. This approach becomes problematic in scenarios with significant protein changes, such as observed in disease-related perturbations like alterations in the secretory pathway, as well as where unconventional protein melting curves (e.g. phase separation proteins) may be expected.

Thermal Tracks is a robust statistical framework designed to improve the analysis of protein melting curves, specifically addressing the challenges of detecting and characterizing system-wide physiological perturbations, such as altered protein glycosylation patterns. Inspired by recent workflows introduced by [Fang et al. (2021)](https://www.nature.com/articles/s42003-021-02306-8), [Pepelnjak et al. (2024)](https://www.nature.com/articles/s41589-024-01568-7) and [Le Sueur et al. (2024)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011632), Thermal Tracks adapts a Bayesian framework to model protein melting curves, offering a robust and flexible approach to analyze complex biological systems. The resulting workflow is fully automated, integrating data preprocessing, statistical analysis, downstream processing, and result visualization into a unified pipeline. Implemented through accessible Jupyter notebooks, it enables researchers to execute complete analyses on standard computing resources without requiring specialized infrastructure.

<p align="center">
  <img width="500" alt="Image" src="https://github.com/user-attachments/assets/626f7723-f462-42c9-908a-063ee1151311" />
</p>

<p align="center"><small><em>Figure 2. Predicted melting curves for membrane proteins based on a previously published TPP experiment (<a href="https://www.embopress.org/doi/full/10.15252/msb.20188242" target="_blank" rel="noopener noreferrer">1</a>) on Escherichia coli lysate with and without MgCl<sub>2</sub>. While NPARC (<a href="https://pubmed.ncbi.nlm.nih.gov/31582558/" target="_blank" rel="noopener noreferrer">2</a>) struggles with fitting the non-typical and complex melting curves of selected membrane proteins, Thermal Tracks enables the accurate fitting of the melting curves, revealing significant differences (*** - BH adjusted p-values < 0.001) in melting behaviors of OmpA. OmpC and LPP upon the addition of MgCl<sub>2</sub>.</em></small></p>

## Installation

### Prerequisites
- Python 3.12 or higher
- [Conda](https://docs.conda.io/en/latest/)
- Jupyter Notebook or JupyterLab
- Clone repository

### Create and activate environment
- conda env create -f environment.yml
- conda activate your-env-name

## Run Thermal Tracks
Thermal Tracks is primarily run through Jupyter notebooks. You can find example analysis notebooks in the [Analysis_notebook](./Analysis_notebooks) folder.
A full tutorial on how to choose parameters and run the analysis will be added shortly.
## License
GNU General Public License v3.0

## Project status
This project is currently **under active development**.
- Resource allocation is not yet optimized; currently, all (available) resources are used, which is suboptimal.
- A comprehensive tutorial is forthcoming to guide users through installation and usage.
- The final, polished version of the code will be publicly released alongside the manuscript shortly.

## Funding
This work was developed during my postdoctoral research in the [Bertozzi Lab](https://bertozzigroup.stanford.edu/) at Stanford University. Development was supported by an [EMBO Postdoctoral Fellowship](https://www.embo.org/funding/fellowships-grants-and-career-support/postdoctoral-fellowships/) and a [Life Sciences Research Foundation (LSRF)](https://lsrf.org/) postdoctoral fellowship funded by the [Howard Hughes Medical Institute (HHMI)](https://www.hhmi.org/).

## Authors and Acknowledgements
This code has been written by Johannes F. Hevler with invaluable input from Shivam Verma.

## Contact
If you have any questions, suggestions, or want to collaborate, feel free to reach out:
- Email: jfhevler@stanford.edu  
