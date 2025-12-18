# Exploring the Drivers of Crimean-Congo Haemorrhagic Fever in Southern Iraq  

Here, we explore the environmental and socio-economic drivers of CCHF, primarily transmitted via _Hyalomma_ ticks, in southeastern Iraq, an area which has seen frequent CCHF resurgences, especially in recent years. 

The repository takes you through the full modelling approach, including data cleaning and harmonisation, preliminary covariate exploration, model fitting and performance evaluation and visualisation of the results. 

All analysis was completed in R. 

The repository is structured in the following way: 
1. 01_data_processing - seperate files for processing each group of covariates
2. 02_data_harmonisation - a script which brought all the individual covariate data files together and harmonised them across a common spatial (admin unit 1) and temporal (weeks) scale, along with the shapefiles used to make the spatial harmonisation (provided by GADM).
3. 03_prelim_analysis - An analysis script which shows the Pearson correlation analysis (to assess multicoliniarity) and LASSO/Elastic Net analysis (to assess covariate performance). A subset of covaraites were selected for model inclusion.
4. 04_main_analyis - fitting of the generalised additive mixed effect models and how the best models were selected
5. 05_analysis_viz - several methods of visualising the best fit models and their results including effect size and residual diagnostics

Developer: Gina E. C. Charnley
