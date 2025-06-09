**Global Sensitivity Analysis of nitrate leaching (L) by NLES5 Model**
This repository contains the scripts for a Global Sensitivity Analysis (GSA) of NLES5, the empirical model that underpins Denmark's national nitrogen leaching assessment tool.

The analysis identifies the relative contribution of each input variable to the model's prediction of nitrate leaching. We applied the Sobol' method to compute first-order and total-order sensitivity indices, which account for both direct and interaction effects. The results were further validated through the use of surrogate models.

Methodological Workflow and Scripts ⚙️
The analysis is structured in the following steps, with each script corresponding to a key part of the workflow.

1. Model Implementation in R

The core NLES5 model functions are coded in 1_fgkNLES5.R.
2. Input Data Management

The script 2_run_franNLES5.R runs the model for the original NLES5 training dataset: Scenarier20190909B4_found0325.xls.
3. Sobol' Sensitivity Analysis

The sampling strategy and calculation of Sobol' indices are performed by the 3_sobol_predus.R script.
4. Validation via Model Surrogation

The surrogate model creation and validation analysis is contained in 4_surrog_predus.R.
