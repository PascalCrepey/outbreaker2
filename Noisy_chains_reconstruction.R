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
load("./fakeMat.Rdata")
load("./data1.Rdata")
load("./chains1.Rdata")

source("./Functions_chains_reconstruction.R")

###########################
#### Global parameters ####
###########################
cores = 2
runs = 100

n_iter_mcmc <- 1000
n_sample <- 10 # n_iter_mcmc*0.0001
burning <- 100 # n_iter_mcmc*0.01

# Compute or not priors for alpha (ancestors) #
prior_alpha <- TRUE

# Minimal support #
min.support <- 10^(-seq(0, 4, by = 0.05))

# Initialization of poisson scale #
init_poisson_scale <- 1
move_poisson_scale <- TRUE

# Adding noise on dates of infection #
adding_noise <- TRUE
lambda_noise <- 3.5

# Other parameters #
move_sigma <- FALSE
init_sigma <- 0.99
move_pi <- TRUE
init_pi <- 1

#############################
#### Preparation of data ####
#############################
detect100 <- data$detect100
chains_detect100 <- chains$detect100

## Use of data with 100% of detection ##
# Dates of detection #
dates <- detect100[, t]
# IDs of hospitals infected or colonized #
ids <- detect100[, hospID]
# Number of cases who can move in the network #
n_cases <- detect100[, colonized + infected]

# Preparation of transfer matrix #
fakeMat <- fakeMat / rowSums(fakeMat)
# To take into account hospitals doing 0 transfer but admitting transferred
# patients
fakeMat[is.nan(fakeMat)] <- 0 

#######################################################
#### Preparation of real transmission chains' data ####
#######################################################
chains_detect100_bind <- rbindlist(chains_detect100)
setnames(chains_detect100_bind, 
         c("hospID","origin"),
         c("to","from"))

# Keeping in mind the chains #
counter <- 0
for(i in 1:chains_detect100_bind[,.N]){
  if(chains_detect100_bind[i, is.na(from)])
    counter <- counter + 1
  chains_detect100_bind[i, chain := counter]
}

# Length of chains #
lengths_chains <- chains_detect100_bind[, .N, by = "chain"]
chains_detect100_bind <- merge(chains_detect100_bind,
                               lengths_chains[, .(chain,
                                                  length_chain = N-1)],
                               by = "chain")

## Preparation of chains_detect100_bind to facilitate the computation of parameters ##
detect100[,num := seq_len(.N)]

chains_detect100_bind.2 <- merge(chains_detect100_bind,
                                 detect100[,.(hospID, t_descendant = t,
                                              t_detect_descendant = t_detect, 
                                              num)],
                                 by.x = c("to", "t_detect"),
                                 by.y = c("hospID","t_detect_descendant"))
chains_detect100_bind.2 <- merge(chains_detect100_bind.2,
                                 detect100[,.(hospID, 
                                              t_ancestor = t, 
                                              num)],
                                 by.x = c("from"),
                                 by.y = c("hospID"),
                                 all.x=TRUE)
chains_detect100_bind.2 <- chains_detect100_bind.2[t_ancestor<t | is.na(from)]
chains_detect100_bind.2 <- chains_detect100_bind.2[order(to, t, t_ancestor)]
chains_detect100_bind.2 <- chains_detect100_bind.2[,.SD[.N],
                                                   by = c("from", "to", "t")]
# chains_detect100_bind.2[, t.y := NULL]
setnames(chains_detect100_bind.2, c("num.x","num.y"), c("ids_to","ids_from"))

detect100[, num := NULL]

chains_detect100_bind <- chains_detect100_bind.2

#################################################
#### Estimating generation time distribution ####
#################################################
chains_detect100_bind[, generation_time := t - t_ancestor]

hist(chains_detect100_bind$generation_time)

generation_time.dist <- as.data.table(prop.table(table(chains_detect100_bind$generation_time)))
setnames(generation_time.dist, 
         c("V1", "N"),
         c("generation_time", "p"))

w <- generation_time.dist[, p]

################################
#### Chains' reconstruction ####
################################
out = mclapply(1:runs, mc.cores = cores, FUN = function(line) {
  output = ChainsReconstruction(dates = dates, 
                                w = w, 
                                n_cases = n_cases, 
                                fakeMat = fakeMat, 
                                ids = ids,
                                detect100 = detect100, 
                                chains_detect100_bind = chains_detect100_bind, 
                                n_iter_mcmc = n_iter_mcmc, 
                                n_sample = n_sample, 
                                burning = burning,
                                min.support = min.support,
                                prior_alpha = prior_alpha,
                                adding_noise = adding_noise, 
                                lambda_noise = lambda_noise,
                                move_sigma = move_sigma,
                                init_sigma = init_sigma,
                                move_pi = move_pi,
                                init_pi = init_pi,
                                init_poisson_scale = init_poisson_scale, 
                                move_poisson_scale = move_poisson_scale)
  return(output)
})

saveRDS(out, file = paste0("1-results_",n_iter_mcmc,"_",n_sample,"_prioralpha.rds"))






