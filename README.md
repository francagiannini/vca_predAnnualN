This document presents an initial Global Sensitivity Analysis focused on NLES5, the empirical model that underpins the Danish national nitrogen leaching calculator, NUAR. The main aim was to illustrate the rela-tive contribution of the various input variables to the model on the variation in the predicted nitrate leach-ing. For Global Sensitivity Analysis, the Sobol' method was applied to compute first- and total-order indi-ces, that account for both direct and interaction effects. Furthermore, to ensure consistency we explored the use of surrogate models for sensitivity analysis. 

Methodological workflow and scripts 

Model Implementation in R: 
fgkNLES5. R here are the core NLES5 function coded

Input data management:
run_franNLES5.R runs the model for the NLES5 training dataset Scenarier20190909B4_found0325.xls

Sampling strategy for Sobol' sensitivity analysis and Sobol' Sensitivity analysis indices calculation:
sobol_predus.R srcipt

Methodological validation trough model surrogation: 
surrog_predus.R
