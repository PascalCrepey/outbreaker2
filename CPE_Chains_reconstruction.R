#############################
#### Loading of packages ####
#############################
library(outbreaker2)
library(data.table)
library(fitdistrplus)
library(parallel)

###################################
#### Loading of simulated data ####
###################################
real_data <- readRDS("./Real_data_CPE.rds")

cpe <- real_data$cpe
transfer_matrix <- real_data$transfer

source("./Functions_chains_reconstruction.R")

###########################
#### Global parameters ####
###########################
cores = 1
runs = 1

n_iter_mcmc <- 10000
n_sample <- n_iter_mcmc*0.0001
burning <- n_iter_mcmc*0.01

# Compute or not priors for alpha (ancestors) #
prior_alpha <- TRUE

# Initialization of poisson scale #
init_poisson_scale <- 1
move_poisson_scale <- TRUE

# Other parameters #
move_sigma <- TRUE
init_sigma <- 0.99
move_pi <- TRUE
init_pi <- 1

#############################
#### Preparation of data ####
#############################
# Dates of detection #
dates <- cpe[, date]
# IDs of hospitals infected or colonized #
ids <- cpe[, finess]
# Number of cases who can move in the network #
n_cases <- cpe[, as.numeric(no_cases)]

# Preparation of transfer matrix #
transfer_matrix <- transfer_matrix / rowSums(transfer_matrix)
# To take into account hospitals doing 0 transfer but admitting transferred
# patients
transfer_matrix[is.nan(transfer_matrix)] <- 0 

#################################################
#### Generating generation time distribution ####
#################################################
w <- dgamma(1:70, shape = 7, rate = 1)

################################
#### Chains' reconstruction ####
################################
output <- RealChainsReconstruction(dates = dates, 
                                   w = w, 
                                   n_cases = n_cases, 
                                   transfers = transfer_matrix, 
                                   ids = ids,
                                   imported = cpe[, imported], 
                                   n_iter_mcmc = n_iter_mcmc, 
                                   n_sample = n_sample, 
                                   burning = burning,
                                   prior_alpha = prior_alpha,
                                   move_sigma = move_sigma,
                                   init_sigma = init_sigma,
                                   move_pi = move_pi,
                                   init_pi = init_pi,
                                   init_poisson_scale = init_poisson_scale, 
                                   move_poisson_scale = move_poisson_scale)




