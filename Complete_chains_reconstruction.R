#############################
#### Loading of packages ####
#############################
library(outbreaker2)
library(data.table)
library(fitdistrplus)

###################################
#### Loading of simulated data ####
###################################
load("../StageCRENet/hackathon/outbreaker/fakeMat.Rdata")
load("../StageCRENet/hackathon/outbreaker/data/data1.Rdata")
load("../StageCRENet/hackathon/outbreaker/chains/chains1.Rdata")

source("./Functions_chains_reconstruction.R")

###########################
#### Global parameters ####
###########################
# # Allowed difference between time of infection found in bayesian output and time of detection #
# allowed.diff <- c(0, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.0)
# Minimum support #
# min.support <- c(0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90)
# Number of MCMC iterations #
# n_iter_mcmc <- c(1000,2500,5000,10000,20000,30000,50000)
n_iter_mcmc <- 5000
n_sample <- 1

# Compute or not priors for alpha (ancestors) #
prior_alpha <- TRUE

# Minimal support #
min.support <- 10^(-seq(0, 8, by = 0.05))

# Initialization of poisson scale #
init_poisson_scale <- 1

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
n_cases <- detect100[, colonized]

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
# Data #
data_outbreaker <- outbreaker_data(dates = dates,
                                   w_dens = w,
                                   n_cases = n_cases,
                                   hosp_matrix = fakeMat,
                                   ids = ids)

## Identification of imported cases ##
check.importation <- merge(detect100[, n_line := seq_len(.N)],
                           chains_detect100_bind,
                           by.x = c("hospID", "t_detect"),
                           by.y = c("to", "t_detect"))
detect100[check.importation[is.na(from) & imported == 0, n_line], 
          imported := 1]
detect100[, n_line := NULL]
imported <- ifelse(detect100[, imported] == 1 | 
                     detect100[, t_detect] == 1, NA_integer_, 1)

if(prior_alpha == T){
  ## Computing priors for alpha ##
  imported <- sapply(seq_len(length(imported)), 
                     FUN = prior_ancestor,
                     imported = imported,
                     data_outbreaker = data_outbreaker,
                     fakeMat = fakeMat)
}

# Config parameters #
config <- create_config(prior_poisson_scale = c(5, 2),
                        move_poisson_scale = TRUE,
                        # init_potential_colonised = rep(2, data_outbreaker$N),
                        # sd_potential_colonised = 5,
                        pb = TRUE,
                        find_import = FALSE,
                        outlier_threshold = 5,
                        data = data_outbreaker,
                        init_tree = imported,
                        n_iter = n_iter_mcmc, 
                        sample_every = n_sample,
                        init_poisson_scale = init_poisson_scale)

# Reconstruction of chains #
results_mcmc <- ComputeBayesian(outbreaker_data = data_outbreaker, 
                                n_iter_mcmc = n_iter_mcmc,
                                ids = ids, 
                                config = config)


##################################
#### Estimation of parameters ####
##################################
parameters <- ComputeParameters(results_bayesian = results_mcmc,
                                data_outbreaker = data_outbreaker,
                                real_data = chains_detect100_bind,
                                min.support = min.support,
                                burning = 5000,
                                init_alpha = imported)

## Shannon entropy ##
results_bayesian_burning <- CreateOutputBayesian(results_mcmc$res, ids, 
                                                 burning = 5000, 
                                                 init_alpha = imported)

summary(results_bayesian_burning$aa[to %in% chains_detect100_bind[,ids_to],
                            .(result=-sum(support*log(support))), by="to"])

save(results_mcmc, parameters,
     file = paste0("./tmp_results/20190809/results_mcmc_",n_iter_mcmc,"_",n_sample,"_prioralpha_poisson_scale1.RData"))

############################
#### Plot of ROC curves ####
############################
PlotROC(se = parameters$parameters_links_aa$se_links.aa,
        sp = parameters$parameters_links_aa$sp_links.aa,
        type.plot = "ggplot")

Plot_se_ppv(se = parameters$parameters_links_aa$se_links.aa,
            ppv = parameters$parameters_links_aa$ppv_links.aa,
            type.plot = "ggplot")

# PlotChains(parameters_chains = parameters$parameters_chains_aa,
#            minimal_support = parameters$minimal_support,
#            type.plot = "ggplot")





