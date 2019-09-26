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
n_iter_mcmc <- 250000
n_sample <- n_iter_mcmc*0.0001
burning <- n_iter_mcmc*0.01

# Compute or not priors for alpha (ancestors) #
prior_alpha <- TRUE

# Minimal support #
min.support <- 10^(-seq(0, 4, by = 0.05))

# Initialization of poisson scale #
init_poisson_scale <- 1

# Adding noise on dates of infection #
adding_noise <- TRUE
lambda_noise <- 3.5

detect100 <- data$detect100
chains_detect100 <- chains$detect100

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


################################################
#### Identification of episodes to undetect ####
################################################
## Identifying the depth of each episode in its propagation chain ##
setkey(chains_detect100_bind,chain,t)

chains_detect100_bind[is.na(from),
                      depth := 0]
index <- 1

while(chains_detect100_bind[is.na(depth), .N] > 0){
  chains_detect100_bind<-merge(chains_detect100_bind,
                               chains_detect100_bind[, .(chain, ids_to, 
                                                         depth_ancestor = depth)],
                               by.x = c("chain", "ids_from"),
                               by.y = c("chain", "ids_to"),
                               all.x = TRUE)
  chains_detect100_bind[depth_ancestor == index - 1,
                        depth := depth_ancestor + 1]
  chains_detect100_bind[, depth_ancestor := NULL]
  index <- index + 1
}

## Identifying the list of ancestors and descendants ##
list_ancestors <- chains_detect100_bind[!is.na(from), unique(ids_from)]
list_descendants <- chains_detect100_bind[!is.na(from), unique(ids_to)]

list_both <- list_ancestors[which(list_ancestors %in% list_descendants)]

## Deletion of lines to undetect ##
# Identification in chains_detect100_bind #
deleted_ids <- chains_detect100_bind[depth %% 2 == 1 &
                                       ids_to %in% list_both, 
                                     ids_to]

# Deletion in detect100 #
detect100 <- detect100[!num %in% deleted_ids]

# Isolation in chains_detect100_bind #
delete_descendant <- chains_detect100_bind[ids_to %in% deleted_ids]
delete_ancestor <- chains_detect100_bind[ids_from %in% deleted_ids]
# Deletion in chains_detect100_bind #
chains_detect100_bind <- chains_detect100_bind[!(ids_to %in% deleted_ids |
                                                 ids_from %in% deleted_ids)]

## Preparation of new links to add ##
new_chains_to_add <- merge(delete_ancestor[, .(chain, ids_from, ids_to, to, t, t_detect,
                                               infected, detected, colonized, length_chain,
                                               t_descendant, generation_time, depth)],
                           delete_descendant[, .(chain, ids_to, ids_from, from)],
                           by.x = c("chain", "ids_from"),
                           by.y = c("chain", "ids_to"))
new_chains_to_add[,ids_from := NULL]
setnames(new_chains_to_add, "ids_from.y", "ids_from")

## Merging the two databases ##
chains_detect100_bind <- rbind(chains_detect100_bind,
                               new_chains_to_add, 
                               fill = TRUE)

## Compute new ids ##
detect100[, new_ids := seq_len(.N)]

chains_detect100_bind <- merge(chains_detect100_bind,
                               detect100[, .(num, new_ids_from = new_ids)],
                               by.x = "ids_from",
                               by.y = "num",
                               all.x = TRUE)
chains_detect100_bind <- merge(chains_detect100_bind,
                               detect100[, .(num, new_ids_to = new_ids)],
                               by.x = "ids_to",
                               by.y = "num",
                               all.x = TRUE)
chains_detect100_bind[, c("ids_to", "ids_from") := NULL]
setnames(chains_detect100_bind, 
         c("new_ids_from", "new_ids_to"),
         c("ids_from", "ids_to"))

#############################
#### Preparation of data ####
#############################
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

################################
#### Chains' reconstruction ####
################################
## Adding noise on dates if needed ##
if(adding_noise){
  dates <- round(dates + rpois(length(dates), 
                               lambda = lambda_noise), 0)
  dates[which(dates < 0)] <- 0
}

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
config <- create_config(prior_poisson_scale = c(1, 1),
                        move_poisson_scale = TRUE,
                        init_potential_colonised = n_cases*init_poisson_scale,
                        # sd_potential_colonised = 5,
                        pb = TRUE,
                        find_import = FALSE,
                        outlier_threshold = 5,
                        data = data_outbreaker,
                        init_tree = imported,
                        n_iter = n_iter_mcmc, 
                        sample_every = n_sample,
                        init_poisson_scale = init_poisson_scale,
                        move_sigma = TRUE,
                        init_sigma = 0.9,
                        move_pi = TRUE,
                        init_pi = 1)

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
                                burning = burning,
                                init_alpha = imported)

save(results_mcmc, parameters,
     file = paste0("./tmp_results/20190902/1-results_mcmc_",n_iter_mcmc,"_",n_sample,"_prioralpha_1.RData"))

############################
#### Plot of ROC curves ####
############################
PlotROC(se = parameters$parameters_links_aa$se_links.aa,
        sp = parameters$parameters_links_aa$sp_links.aa,
        min.support = min.support,
        type.plot = "ggplot")

Plot_se_ppv(se = parameters$parameters_links_aa$se_links.aa,
            ppv = parameters$parameters_links_aa$ppv_links.aa,
            type.plot = "ggplot")
















