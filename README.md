**Global Sensitivity Analysis of nitrate leaching (L) by NLES5 Model**
This repository contains the scripts for a Global Sensitivity Analysis (GSA) of NLES5, the empirical model that underpins Denmark's national nitrogen leaching assessment tool.

The analysis identifies the relative contribution of each input variable to the model's prediction of nitrate leaching. We applied the Sobol' method to compute first-order and total-order sensitivity indices, which account for both direct and interaction effects. The results were further validated through the use of surrogate models.

The report can be downloaded from Aarhus University Pure system: 

[Giannini-Kurina, F., Ugilt Larsen, S., & Eriksen, J., (2025). Sensitivity Analysis of the NLES5 Model for the NUAR Calculator, Nr. 2025-0832229, 33 s., jun. 18, 2025.](https://pure.au.dk/portal/da/publications/sensitivity-analysis-of-the-nles5-model-for-the-nuar-calculator)


Methodological Workflow and Scripts ⚙️
The analysis is structured in the following steps, with each script corresponding to a key part of the workflow.

* **1. Model Implementation in R**
    * The core NLES5 model functions are coded in [`1_fgkNLES5.R`](https://github.com/francagiannini/vca_predAnnualN/blob/main/1_fgkNLES5.R)

* **2. Input Data Management**
    * The script [`2_run_franNLES5.R`](https://github.com/francagiannini/vca_predAnnualN/blob/main/2_run_franNLES5.R) runs the model for the training data

* **3. Sobol' Sensitivity Analysis**
    * The sampling strategy is performed by the [`3_sobol_predus.R`](https://github.com/francagiannini/vca_predAnnualN/blob/main/3_sobol_predus.R) script.

* **4. Validation via Model Surrogation**
    * The surrogate models and results is contained in [`4_surrog_predus.R`](https://github.com/francagiannini/vca_predAnnualN/blob/main/4_surrog_predus.R).
