#############################
#### Loading of packages ####
#############################
library(outbreaker2)
library(data.table)
library(ggplot2)
library(cowplot)

#########################
#### Loading of data ####
#########################
results <- readRDS("./tmp_results/Genouest/RealChainsReconstruction.rds")

source("./Functions_chains_reconstruction.R")

min_support <- 0.063
gamma_shape <- 4
max_episodes <- 3054

###########################################
#### Modification of the burnin period ####
###########################################
results <- lapply(results, function(r, burning){
  outputBayesian <- CreateOutputBayesian(results_bayesian = r$output$res,
                                         burning = burning,
                                         init_alpha = r$output$res_consensus$init_alpha,
                                         ids = r$output$res_consensus$to_label)
  names(outputBayesian) <- c("res_consensus", "res_aa")
  outputBayesian$res <- r$output$res
  
  return(list(output = outputBayesian,
              rate = r$rate))
  },
  burning = 5000
  )

###################################
#### Computation of parameters ####
###################################
gamma_rate <- unlist(lapply(results, 
                           function(r) {
                             r$rate
                             }
                           )
                    )

generation_time_mean <- gamma_shape / gamma_rate
  
no_links_retrieved <- unlist(lapply(results, 
                                    function(r) {
                                      r$output$res_consensus[is.na(from) | (support >= min_support &
                                                                              kappa %in% 1:3), .N]
                                      }   
                                    ))

community_episodes <- unlist(lapply(results, 
                                    function(r) {
                                      r$output$res_consensus[is.na(from) & 
                                                               !is.na(init_alpha), .N]
                                      }   
                                    ))

imported_episodes <- unlist(lapply(results, 
                                   function(r) {
                                     r$output$res_consensus[is.na(from) & 
                                                              is.na(init_alpha), .N]
                                     }   
                                   ))

data_plot <- data.table(generation_time_mean = generation_time_mean,
                        no_links_retrieved = no_links_retrieved,
                        community_episodes = community_episodes,
                        imported_episodes = imported_episodes)

FigPlot_nolinks <- ggplot(data_plot, aes(x = generation_time_mean)) +
  geom_line(aes(y = no_links_retrieved)) +
  geom_point(aes(y = no_links_retrieved)) + 
  geom_vline(xintercept = generation_time_mean[which.max(no_links_retrieved)],
             linetype="dashed", 
             color = "black",
             size = 1.1) +
  theme_minimal() +
  xlab("Mean of CPE generation interval") + 
  ylab("#Episodes with links retrieved") +
  theme(axis.title.x = element_text(size = 50),
        axis.title.y = element_text(size = 50),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40))
FigPlot_nolinks

# FigPlot_community <- ggplot(data_plot, aes(x = generation_time_mean)) +
#   geom_line(aes(y = community_episodes)) +
#   geom_point(aes(y = community_episodes)) + 
#   geom_vline(xintercept = generation_time_mean[which.max(no_links_retrieved)],
#              linetype="dotted", 
#              color = "blue",
#              size = 1.05) +
#   theme_minimal()
# FigPlot_community

FigPlot <- ggplot(data_plot, aes(x = generation_time_mean)) +
  geom_ribbon(aes(ymin = 0, 
                  ymax = imported_episodes / max_episodes * 100, 
                  fill = "imported"))+
  geom_ribbon(aes(ymin = imported_episodes / max_episodes * 100, 
                  ymax = (imported_episodes + community_episodes) / max_episodes * 100, 
                  fill= "community")) +
  geom_ribbon(aes(ymin = (imported_episodes + community_episodes) / max_episodes * 100, 
                  ymax = no_links_retrieved / max_episodes * 100, 
                  fill = "links_retrieved")) +
  geom_ribbon(aes(ymin = no_links_retrieved / max_episodes * 100, 
                  ymax = 100, 
                  fill = "not_retrieved")) +
  scale_fill_manual(name="",
                    values=c('imported'='#0056A3', 'community'='#1F83AE',
                             'links_retrieved'='#A5C3DE', 'not_retrieved'='#F0F8FF'),
                    labels=c('imported'='Imported episodes',
                             'community'='Community-imported episodes',
                             'links_retrieved'='Episodes with enough support',
                             'not_retrieved'='Episodes without enough support')) +
  geom_vline(xintercept = generation_time_mean[which.max(no_links_retrieved)],
             linetype="dashed", 
             color = "black",
             size = 1.1) + 
  xlab("Mean of CPE generation interval") +
  ylab("% of links") +
  theme_minimal() +
  theme(legend.position="bottom") +
  theme(axis.title.x = element_text(size = 50),
        axis.title.y = element_text(size = 50),
        axis.text.x = element_text(size = 40),
        axis.text.y = element_text(size = 40),
        legend.text = element_text(size = 40))

plot_grid(FigPlot_nolinks,FigPlot,
          labels=c("A","B"), nrow=2,
          align = "v", label_x = 0.95, 
          label_size = 50)
