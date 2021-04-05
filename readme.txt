This repository contains the code needed to perform our finite element and Sobol sensitivity analysis for cone photoreceptors in bright light. We highlight the following for making use of the code as in our paper.

1.) main.m
This is the main driver script for most of the modeling tasks. These include: the activation signaling cascade on discs (Act), the conversion of E* as a surface quantity to a volumic one (Evol), the 2nd messenger diffusion of cGMP and Ca2+ in the cytosol. To run basic flash simulations, these three cases should be run in sequence. In addition to these basic cases, main.m can also perform several sensitivity and fitting tasks: gradient based local sensitivity about a choice of parameters (SA_pd), Monte Carlo simulation for Sobol, and the MCMC fitting of experimental flash responses (RanW). For Sobol analysis, supercomputing resources were required for adequate numbers of simulation for Monte Carlo and processing.

2.) 'common/data_set.m'
This is where model parameters for basic photoreceptor simulations are specified. For fitting and sensitivity analysis, parameters that are not varied are taken from this file. 

3.) 'sensitivity analysis/param fit'
This directory contains the files for where the MCMC random walk parameters are specified: 'mcmc_alpha.m', and 'mcmc_ranw.m'. In addition 'fitness.m' specifies the experimental data one is trying to fit as an approximating polynomial.  These polynomials were extracted from image data. 

4.) 'sensitivity analysis/Sobol/mc_sample_Sbl.m'
This script is where uniform ranges for parameters, to be used in Sobol analysis, may be specified.

5.) prepr_Sbl.m
This script should be run before running main.m with the Sobol option to prepare the random parameter samples which will be used in Monte Carlo simulation

6.) 'post process/postpr_SASbleval.m,postpr_SASbl.m'
The first script should be run immediately after main.m has concluded to generate output required to run next the second script.  The second script generates Sobol csv files and figures. Whether to plot QQ-plots as well may be done by setting flag_plot in 'post process/mc_plot.m'. These scripts require that the files Sbl_Evals.mat and trial_master.mat both be in the post process folder before running.  These files may be downloaded from this paper's associated Dryad repository.
