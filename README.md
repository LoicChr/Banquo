# Banquo

This code is associated with the article: 'Integrating traits, competition and abiotic filtering improves biodiversity predictions' by Loïc Chalmandrier*, Daniel B. Stouffer*, Adam S. T. Purcell, William G. Lee, Andrew J. Tanentzap and Daniel C. Laughlin*
 
 *code developpers
  
The code models species abundances using functional traits, community data and environmental data. It assumes that species abundances are constrained by abiotic filtering (Traitspace model) and then competition (Banquo model). First species carrying capacities are modeled through the Traitspace model (using a traits-environment multivariate linear model) and then pairwise competitive interaction are calibrated using the Banquo model through a nloptr algorithm.


System requirements
R 4.1.0

The repository contains the following files and folders.
./data/community.cover
./lib/abgFunctions.R: functions to calculate community alpha and beta diversity
     /interactionMatrix.R: functions to calculate the pairwise interaction matrix
     /Likelihood.R: functions to calculate the likelihood of the Banquo models
     /traitspace_and_banquo.R: functions to compute Traitspace and to compute Banquo
     
./main/analysis_results.R: script using the final output to produce the table and figures of the article
      /assembly_models.R: Main script. It uses the output of traitspace_gen.R and compute the different Banquo models. 
      /data_prep.R: script to prepare the data.
      /priors.R: script to create the prior density function, the prior sampling function and parameter bounds
      /traitspace_gen.R: script to compute Traitspace from the data.
./results/*.eps: Figures of the main article.
         /*.jpg: Figures of the supplementary material.
         /tab.csv: Table 1 of the main article
         /abio_TRUE/*/posterior_objs.Rdata: contains the output of the script assembly_models.R. It contains most notably the posterior distribution of the abiotic and the abiotic + biotic models. The objects are saved in subfolder named after the traits used to calibrate biotic interactions (including no interactions in the folder called 'none')
         /abio_TRUE/*/posterior_objs.Rdata: contains the output of the script assembly_models.R. It contains most notably the posterior distribution of the null model and the biotic models. The objects are saved in subfolder named after the traits used to calibrate biotic interactions (including no interactions in the folder called 'none')
         /post_models/post_Model_*.jpg: graphics of the posterior distribution of each assembly model.