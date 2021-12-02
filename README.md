# Banquo

This code is associated with the article: 'Integrating traits, competition and abiotic filtering improves biodiversity predictions' by Loïc Chalmandrier\*, Daniel B. Stouffer, Adam S. T. Purcell, William G. Lee, Andrew J. Tanentzap and Daniel C. Laughlin
 
\* main code developper and contact. (loic.chalmandrier@biologie.uni-regensburg.de)

## Overview  
The code models species abundances using functional traits, community data and environmental data. It assumes that species abundances are constrained by abiotic filtering (Traitspace model) and then competition (Banquo model). First species carrying capacities are modeled through the Traitspace model (using a trait-environment multivariate linear model) and then pairwise competitive interaction are calibrated using the Banquo model through a nloptr algorithm.


### System requirements
R 4.1.0

### Libraries
BayesianTools, corpcor, DEoptim, extraDistr, ggpubr, mclust, mgcv, plyr, pROC, purrr, tidyverse


## Repository content

<p> './data/community.cover.23sp.csv': community matrix - site by species matrix <br>
'data/cover_class.txt': contains the names and the minimum and maximum cover associated to each cover class <br>
'data/data_ready.Rdata': contains the formatted data (output of the script data_prep.R) <br>
'kettleholedata.23sp.csv': contains the functional trait data in each site and for each species <br>
'traitspace_objs.Rdata': contains the output of the script traitspace_gen.R </p>
<p>  './lib/abgFunctions.R': functions to calculate community alpha and beta diversity <br>
'./lib/interactionMatrix.R': functions to calculate the pairwise interaction matrix <br>
'./lib/Likelihood.R': functions to calculate the likelihood of the Banquo models<br> 
'./lib/traitspace_and_banquo.R': functions to compute Traitspace and to compute Banquo</p>
<p> './main/analysis_results.R': script using the final output to produce the table and figures of the article <br>
'./main/assembly_models.R': Main script. It uses the output of traitspace_gen.R and compute the different Banquo models.<br>
'./main/data_prep.R': script to prepare the data.<br>
'./main/priors.R': script to create the prior density function, the prior sampling function and parameter bounds <br>
'./main/traitspace_gen.R': script to compute Traitspace from the data.</p>
<p>  './results/\*.eps': Figures of the main article.<br>
'./results/\*.jpg': Figures of the supplementary material.<br>
'./results/tab.csv': Table 1 of the main article<br>
'./results/abio_TRUE/\*/posterior_objs.Rdata': it contains most notably the posterior distribution of the abiotic and the abiotic + biotic models. The objects are saved in subfolder named after the traits used to calibrate biotic interactions (including no interactions in the folder called 'none')<br>
'./results/abio_TRUE/\*/posterior_objs.Rdata': it contains most notably the posterior distribution of the null model and the biotic models. The objects are saved in subfolder named after the traits used to calibrate biotic interactions (including no interactions in the folder called 'none')<br>
'./results/post_models/post_Model_\*.jpg': graphics of the posterior distribution of each assembly model.</p>
         
## How to use the repository
1. set the working directory to the root of the repository.
2. run data_prep.R: from the raw data, this script prepares the 'data/data_ready.Rdata' object
3. run traitspace_gen.R: this will use the content of 'data/data_ready.Rdata' to compute Traitspace and create 'data/traitspace_objs.Rdata' 
4. run assembly_models.R: this script will run the Banquo models (more details of its functioning in the header of the script). after the models are run, the script saves the output in the 'posterior_objs.Rdata' objects
5. run analysis_results.R: this script uses the posterior_objs.Rdata, the function in ./lib/ and the objects in ./data/ to create the figures and table from the article and the supplementary materials.
