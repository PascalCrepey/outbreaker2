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
cores = 10
runs = 100

n_iter_mcmc <- 50000
n_sample <- n_iter_mcmc*0.0001
burning <- n_iter_mcmc*0.01

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
move_sigma <- TRUE
init_sigma <- 0.99
move_pi <- TRUE
init_pi <- 1

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

saveRDS(out, file = paste0("1-results_incomplete_",n_iter_mcmc,"_",n_sample,"_prioralpha.rds"))




