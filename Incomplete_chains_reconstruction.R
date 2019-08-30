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
n_iter_mcmc <- 500000
n_sample <- 50

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

# Preparation of transfer matrix #
fakeMat <- fakeMat / rowSums(fakeMat)
# To take into account hospitals doing 0 transfer but admitting transferred
# patients
fakeMat[is.nan(fakeMat)] <- 0 

########################################################
#### Estimating the distribution of generation time ####
########################################################
## Using the chains results ##
chains_detect100 <- lapply(chains_detect100, generation_time)
# chains_detect100 <- lapply(chains_detect100, get_origin)

chains_detect100_bind <- rbindlist(chains_detect100)
setnames(chains_detect100_bind, 
         c("hospID","origin"),
         c("to","from"))

hist(chains_detect100_bind$w)

dist.w <- fitdist(chains_detect100_bind[!is.na(w),w],"nbinom")

w <- dnbinom(1:50, 
             size = dist.w$estimate[1],
             mu = dist.w$estimate[2])

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


################################
#### "Undetecting" episodes ####
################################
chains_detect100_bind[, no_line := .I]

data_undetected <- merge(chains_detect100_bind[, .(chain, ancestor = from, 
                                                   t, t_detect, hospID = to,
                                                   no_line_hospID = no_line)],
                         chains_detect100_bind[, .(chain, hospID = from, t_descendant = t,
                                                   t_detect_descendant = t_detect,
                                                   descendant = to,
                                                   no_line_descendant = no_line)],
                         by.x = c("chain", "hospID"),
                         by.y = c("chain", "hospID"),
                         all.x = TRUE)
data_undetected <- data_undetected[t_descendant > t &
                                     !is.na(ancestor) & 
                                     descendant != ancestor,
                                   .SD[.N],
                                   by=c("chain","hospID","descendant","ancestor")]
to_undetect <- data_undetected[,unique(no_line_hospID)]

detect100[, no_line := seq_len(.N)]
tmp <- merge(detect100,
             chains_detect100_bind[no_line %in% to_undetect,
                                   .(chain, to, t, t_detect, colonized, infected, detected, from, t_lag, w, length_chain)],
             by.x = c("hospID", "t", "t_detect"),
             by.y = c("to", "t", "t_detect"))

detect100.modified <- detect100[!no_line %in% tmp[,no_line]]
detect100.modified[, no_line := NULL]

detect100.modified[, no_line := seq_len(.N)]

# Dates of detection #
dates <- detect100.modified[, t_detect]
# IDs of hospitals infected or colonized #
ids <- detect100.modified[, hospID]
# Number of cases who can move in the network #
n_cases <- detect100.modified[, colonized]

##############################################################
# ## Use of data with 100% of detection ##
# setkey(chains_detect100_bind, chain, t)
# 
# chains_detect100_bind[, no_line := .I]
# # First and last line of a chain #
# list_first_last <- chains_detect100_bind[,
#                                          .SD[c(1,.N)],
#                                          by = "chain"][, unique(no_line)] 
# chains_detect100_bind[, I_first_last := (no_line %in% list_first_last)]
# 
# # Lines to "undetect" #
# chains_detect100_bind[, to_lag := c(NA, 
#                                       chains_detect100_bind[1:(.N-1),to])]
# chains_detect100_bind[, I_to_undetect := FALSE]
# chains_detect100_bind[to_lag == from &
#                         length_chain > 1 &
#                         I_first_last == FALSE, 
#                       I_to_undetect := TRUE]
# chains_detect100_bind[, no_line := NULL]
# 
# detect100[, no_line := seq_len(.N)]
# tmp <- merge(detect100,
#              chains_detect100_bind[I_to_undetect == TRUE,
#                                    .(chain, to, t, t_detect, colonized, infected, detected, from, t_lag, w, length_chain)],
#              by.x = c("hospID", "t", "t_detect"),
#              by.y = c("to", "t", "t_detect"))
# 
# detect100.modified <- detect100[!no_line %in% tmp[,no_line]]
# detect100.modified[, no_line := NULL]
# 
# detect100.modified[, no_line := seq_len(.N)]
# 
# # Dates of detection #
# dates <- detect100.modified[, t_detect]
# # IDs of hospitals infected or colonized #
# ids <- detect100.modified[, hospID]
# # Number of cases who can move in the network #
# n_cases <- detect100.modified[, colonized]
#########################################################################



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
check.importation <- merge(detect100.modified,
                           chains_detect100_bind[,-12],
                           by.x = c("hospID", "t_detect"),
                           by.y = c("to", "t_detect"))
detect100.modified[check.importation[is.na(from) & imported == 0, no_line], 
          imported := 1]
imported <- ifelse(detect100.modified[, imported] == 1 | 
                     detect100.modified[, t_detect] == 1, NA_integer_, 1)

if(prior_alpha == T){
  ## Computing priors for alpha ##
  imported <- sapply(seq_len(length(imported)), 
                     FUN = prior_ancestor,
                     imported = imported,
                     data_outbreaker = data_outbreaker,
                     fakeMat = fakeMat)
}

# Config parameters #
config <- create_config(init_potential_colonised = rep(2, data_outbreaker$N),
                        sd_potential_colonised = 5,
                        prior_poisson_scale = c(5, 2),
                        move_poisson_scale = FALSE,
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

## Preparation of chains_detect100_bind to facilitate the computation of parameters ##
# chains_tomerge <- chains_detect100_bind[!no_line %in% 
#                                           data_undetected[,
#                                                           c(no_line_hospID,no_line_descendant)],
#                                         -c(9,10,11,12)]
# chains_tomerge <- rbind(chains_tomerge,
#                         data_undetected[, .(from = ancestor,
#                                             to = descendant,
#                                             chain, 
#                                             t = t_descendant,
#                                             t_detect = t_detect_descendant,
#                                             detected = 2)],
#                         fill = TRUE)
# setkey(chains_tomerge, chain, t)
# 
chains.undetected <- merge(data_undetected,
                           detect100.modified[,.(hospID, t_detect, no_line)],
                           by.x = c("descendant", "t_detect_descendant"),
                           by.y = c("hospID","t_detect"))
chains.undetected <- merge(chains.undetected,
                           detect100.modified[,.(hospID, t_ancestor = t, no_line)],
                           by.x = c("ancestor"),
                           by.y = c("hospID"))
chains.undetected <- chains.undetected[t_ancestor<t_descendant]
chains.undetected <- chains.undetected[order(descendant, t_descendant, t_ancestor)]
chains.undetected <- chains.undetected[,.SD[.N],
                                       by = c("ancestor", "descendant", "t_descendant")]
chains.undetected[, t_ancestor := NULL]
setnames(chains.undetected, c("no_line.x","no_line.y"), c("ids_descendant","ids_ancestor"))

comparison <- merge(results_mcmc$res_aa,
                    chains.undetected[,.(ids_ancestor, ids_descendant, 
                                         ancestor, descendant, t_descendant, chain)],
                    by.x = "to",
                    by.y = "ids_descendant")

comparison[ids_ancestor == from, .N]

##################################
#### Estimation of parameters ####
##################################
parameters <- ComputeParameters(results_bayesian = results_mcmc,
                                data_outbreaker = data_outbreaker,
                                real_data = chains_detect100_bind,
                                min.support = min.support,
                                burning = 100,
                                init_alpha = imported)

## Shannon entropy ##
results_bayesian_burning <- CreateOutputBayesian(results_mcmc$res, ids, 
                                                 burning = 100, 
                                                 init_alpha = imported)

summary(results_bayesian_burning$aa[to %in% chains_detect100_bind[,ids_to],
                                    .(result=-sum(support*log(support))), by="to"])

# save(results_mcmc, parameters, 
#      file = paste0("./tmp_results/20190809/results_mcmc_",n_iter_mcmc,"_",n_sample,"_prioralpha_poisson_scale1.RData"))

############################
#### Plot of ROC curves ####
############################
PlotROC(se = parameters$parameters_links_aa$se_links.aa,
        sp = parameters$parameters_links_aa$sp_links.aa,
        type.plot = "ggplot")

PlotChains(parameters_chains = parameters$parameters_chains_aa,
           minimal_support = parameters$minimal_support,
           type.plot = "ggplot")





