#########################
#### Loading of data ####
#########################
results <- readRDS("./tmp_results/Genouest/1-results_50000_5_parLapply.rds")

##########################
#### Global variables ####
##########################
burning <- max(results[[1]]$results_mcmc$res$step) * 0.01

###################################
#### Computation of parameters ####
###################################
### Global parameters ####
# Pi #
pi <- unlist(lapply(results, 
                    function(r) {
                      r$results_mcmc$res[r$results_mcmc$res$step>burning,]$pi
                      }    
                    )
             )
quantile(pi, c(0.025,0.5,0.975))

# Sigma #
sigma <- unlist(lapply(results, 
                    function(r) {
                      r$results_mcmc$res[r$results_mcmc$res$step>burning,]$sigma
                      }    
                    )
                )
quantile(sigma, c(0.025,0.5,0.975))

# Poisson scale #
poisson_scale <- unlist(lapply(results, 
                               function(r) {
                                 r$results_mcmc$res[r$results_mcmc$res$step>burning,]$poisson_scale
                                 }   
                               )
                        )
quantile(poisson_scale, c(0.025,0.5,0.975))

#### A revoir en s'inspirant ddes calculs sur pi, sigma et poisson scale ######
# Median Shannon entropy #
# Already computed using burning period #
median_shannon_entropy <- unlist(lapply(results, 
                                    function(r) {
                                      as.numeric(r$parameters$shannon_entropy[3])
                                      }    
                                    )
                             )

### Consensus ancestor ###
# Already computed using burning period #
# Number of true positive links #
tp_links.consensus <- unlist(lapply(results, 
                                    function(r) {
                                      r$parameters$parameters_links_consensus$tp_links.consensus
                                      }                   
                                    )
                             )

# Sensitivity #
se_links.consensus <- unlist(lapply(results, 
                                    function(r) {
                                      r$parameters$parameters_links_consensus$se_links.consensus
                                      }  
                                    )
                             )

# Percentage of global chains reconstructed #
global_chains_consensus <- unlist(lapply(results, 
                                    function(r) {
                                      r$parameters$parameters_chains_consensus$global_chains_consensus
                                      }  
                                    )
                                  )






