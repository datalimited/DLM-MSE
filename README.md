# DLM-MSE
##Data-limited MSE R project
MSE code for testing catch-only fisheries models and related HCRs

#Developers: Jessica Walsh, Coilin Minto, Ernesto Jardim, Sean Anderson and co-authors on paper

Code for Paper: Walsh et al. Fish and Fisheries 2018 
Title: Trade-offs for data-limited fisheries when using harvest strategies based on catch-only models
Developed as part of the Data Limited Fisheries Working Group 
funded by Conservation International and Gordon and Betty Moore Foundation.


#1 Simulations script - creates the simulated species 
  # 01-sim-setup.R

#2 Runs MSE using perfect informtion
  # 02a-perfect-info-select-MSE-run.R   #selects which scenario to run, and draws code from 02b_base code.
  # 02b-pefect-info-MSE-base.R  #contains base code for MSE
  # 02c-perfect-info-performance
  # calls functions from scripts 3, 4 and 4a.

#3 Conduct assessment using catch only models - select scenarios to run assessment, and 05b is base code
  # 03a-preassess-run.R
  # 03b-dl-preassess-base.R

#4 Conduct MSE with COM assessments
  # 04a-MSE-run.R 
  # 04b-MSE-base.R (run using 02a-MSE-run.R) 

#5 Diagnostic plots to check simulations and MSE running
  # 05a-plots-simulations.Rmd
  # 05b_mse_diagnostic_all_sp.Rmd

#6 Performance metrics
  # 06a-performance-metrics.Rmd 
  # 06b-performance-metrics-UR.R #with Underreporting
  

