# Banquo2

This code is associated with the article: 'Integrating traits, competition and abiotic filtering improves biodiversity predictions' by Loïc Chalmandrier*, Daniel B. Stouffer*, Adam S. T. Purcell, William G. Lee, Andrew J. Tanentzap and Daniel C. Laughlin*
 
 *code developpers
  
The code models species abundances using functional traits, community data and environmental data. It assumes that species abundances are constrained by abiotic filtering (Traitspace model) and then competition (Banquo model). First species carrying capacities are modeled through the Traitspace model (using a traits-environment multivariate linear model) and then pairwise competitive interaction are calibrated using the Banquo model through a nloptr algorithm.

main.R is the main script that show the procedure to treat the dataset and reproduces the model used in the original article. Details explanation of each function is provided in their associated script.
