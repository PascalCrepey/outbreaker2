source("./Functions_chains_reconstruction.R")

#### Function to change the burning period of an already computed Bayesian analysis ####
BurninChange <- function(res, burning, init_alpha, ids, real_data, min.support){
  outputBayesian <- CreateOutputBayesian(results_bayesian = res,
                                         burning = burning,
                                         init_alpha = init_alpha,
                                         ids = ids)
  
  names(outputBayesian) <- c("res_consensus", "res_aa")
  outputBayesian$res <- res
  
  parameters <- ComputeParameters(results_bayesian = outputBayesian,
                                  real_data = real_data,
                                  min.support = min.support,
                                  burning = 0,
                                  init_alpha = init_alpha,
                                  ids = ids)
  
  return(list(results_mcmc = outputBayesian,
              parameters = parameters))
}




#### Function to edit the parameters according to the pattern: median [Q1-Q3] ####
ParametersEditing <- function(ResultVector, no_digits, percent){
  quantiles <- round(quantile(ResultVector, 
                              c(0.025,0.5,0.975), 
                              na.rm = TRUE),
                     no_digits) * ifelse(percent, 100, 1)
  return(paste0(quantiles[2], 
                " [", quantiles[1], "-",
                quantiles[3], "]"))
}

#### Function to compute the 95%CI of sensibility and specificity according to minimal support ####
ParametersAccordingSupport <- function(i, results, min.support){
  se_links <- quantile(unlist(lapply(results, 
                                     function(r) {
                                       r$parameters$parameters_links_aa$se_links.aa[i]
                                     }  
  )),
  c(0.025,0.5,0.975))
  sp_links <- quantile(unlist(lapply(results,
                                     function(r) {
                                       r$parameters$parameters_links_aa$sp_links.aa[i]
                                     }  
  )),
  c(0.025,0.5,0.975))
  return(list(min.support = min.support[i],
              se_2.5 = se_links[1],
              se_50 = se_links [2],
              se_97.5 = se_links[3],
              sp_2.5 = sp_links[1],
              sp_50 = sp_links [2],
              sp_97.5 = sp_links[3]))
}

#### Main function to synthetise the different parameters computed on the results of the BDLF ####
ParametersSynthesis <- function(index, vec_burning, location, type,
                                NotRetrieved = FALSE, chains_detect100_bind = NULL,
                                StudyKappa = FALSE, new_chains_added = NULL,
                                NewBurning = FALSE, min.support = NULL){
  cat(paste0("Scenario ", index, "\n"))
  results <- readRDS(paste0("./tmp_results/",location,"/1-results",
                            ifelse(type == "complete", "", paste0("_",type)),
                            "_50000_5_parLapply_scenario", index, ".rds"))
  burning <- vec_burning[index]
  
  ## Computation of parameters using burning period to be sure of burning period used ##
  if(NewBurning){
    results <- lapply(results, 
                              function(r) {
                                BurninChange(res = r$results_mcmc$res,
                                             real_data = chains_detect100_bind,
                                             min.support = min.support,
                                             burning = burning,
                                             init_alpha = r$results_mcmc$res_consensus$init_alpha,
                                             ids = r$results_mcmc$res_consensus$to_label)
                              }   
    )
  }

  ###################################
  #### Computation of parameters ####
  ###################################
  ######
  ### Global parameters ###
  ######
  # Pi #
  pi <- unlist(lapply(results, 
                      function(r) {
                        r$results_mcmc$res[r$results_mcmc$res$step>burning,]$pi
                      }    
  )
  )
  
  # Sigma #
  sigma <- unlist(lapply(results, 
                         function(r) {
                           r$results_mcmc$res[r$results_mcmc$res$step>burning,]$sigma
                         }  
  )
  )
  
  # Poisson scale #
  poisson_scale <- unlist(lapply(results, 
                                 function(r) {
                                   r$results_mcmc$res[r$results_mcmc$res$step>burning,]$poisson_scale
                                 }   
  )
  )
  
  # Median Shannon entropy #
  # Already computed using burning period #
  median_shannon_entropy <- unlist(lapply(results, 
                                          function(r) {
                                            r$results_mcmc$res_aa[to %in%
                                                                    chains_detect100_bind[,ids_to],
                                                                  .(result=-sum(support*log(support))),
                                                                  by="to"][,result]
                                          }    
  )
  )
  
  ######
  ### Consensus ancestor ###
  ######
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
  
  ######
  ### All ancestors ###
  ######
  # Sensitvity #
  se_links.aa <- unlist(lapply(results, 
                               function(r) {
                                 r$parameters$parameters_links_aa$se_links.aa
                               }  
  )
  )
  
  # Specificity #
  sp_links.aa <- unlist(lapply(results, 
                               function(r) {
                                 r$parameters$parameters_links_aa$sp_links.aa
                               }  
  )
  )
  
  # Number of links correctly retrieved #
  tp_links.aa <- unlist(lapply(results, 
                               function(r) {
                                 r$parameters$parameters_links_aa$tp_links.aa
                               }  
  )
  )
  
  # Positive predictive value #
  ppv_links.aa <- unlist(lapply(results, 
                                function(r) {
                                  r$parameters$parameters_links_aa$ppv_links.aa
                                } 
  )
  )
  
  # Negative predictive value #
  npv_links.aa <- unlist(lapply(results, 
                                function(r) {
                                  r$parameters$parameters_links_aa$npv_links.aa
                                }  
  )
  )
  
  # Percentage of entire chains reconstructed #
  global_chains_aa <- unlist(lapply(results, 
                                    function(r) {
                                      r$parameters$parameters_chains_aa$global_chains_aa
                                    }  
  )
  )
  
  ######
  ### Maximum sensitivity ###
  ######
  # Maximum sensitivity #
  se_links.aa_max <- unlist(lapply(results, 
                                   function(r) {
                                     max(r$parameters$parameters_links_aa$se_links.aa)
                                   }  
  )
  )
  
  # ID of the maximum sensitivity #
  se_max_id <- unlist(lapply(results, 
                             function(r) {
                               which(r$parameters$parameters_links_aa$se_links.aa == 
                                       max(r$parameters$parameters_links_aa$se_links.aa))[1]
                             }  
  )
  )
  
  # Specificity related to maximum sensitivity #
  sp_links.aa_related_se <- mapply(function(r, id) {
    r$parameters$parameters_links_aa$sp_links.aa[id]
  },
  results,
  se_max_id)
  
  # Number of links correctly retrieved related to maximum sensitivity #
  tp_links.aa_related_se <- mapply(function(r, id) {
    r$parameters$parameters_links_aa$tp_links.aa[id]
  },
  results,
  se_max_id)
  
  # Minimal support related to maximum sensitivity #
  minimal_support.aa_related_se <- mapply(function(r, id) {
    r$parameters$minimal_support[id]
  },
  results,
  se_max_id)
  
  # Positive predictive value retrieved related to maximum sensitivity #
  ppv_links.aa_related_se <- mapply(function(r, id) {
    r$parameters$parameters_links_aa$ppv_links.aa[id]
  },
  results,
  se_max_id)
  
  # Negative predictive value related to maximum sensitivity #
  npv_links.aa_related_se <- mapply(function(r, id) {
    r$parameters$parameters_links_aa$npv_links.aa[id]
  },
  results,
  se_max_id)
  
  # Percentage of entire chains reconstructed related to maximum sensitivity #
  global_chains_aa_related_se <- mapply(function(r, id) {
    r$parameters$parameters_chains_aa$global_chains_aa[id]
  },
  results,
  se_max_id)
  
  ######
  ### Youden index ###
  ######
  # Youden index #
  youden_index <- unlist(lapply(results, 
                                function(r) {
                                  r$parameters$parameters_links_aa$se_links.aa +
                                    r$parameters$parameters_links_aa$sp_links.aa -
                                    1
                                }   
  )
  )
  
  # Index of Youden index max #
  youden_max <- unlist(lapply(results, 
                              function(r) {
                                which.max(r$parameters$parameters_links_aa$se_links.aa +
                                            r$parameters$parameters_links_aa$sp_links.aa -
                                            1)
                              }   
  )
  )
  
  # Sensitivity linked to Youden max #
  se_youden <- mapply(function(r, id) {
    r$parameters$parameters_links_aa$se_links.aa[id]
  },
  results,
  youden_max)
  
  # Specificity linked to Youden max #
  sp_youden <- mapply(function(r, id) {
    r$parameters$parameters_links_aa$sp_links.aa[id]
  },
  results,
  youden_max)
  
  # Number of links correctly retrieved related to Youden max #
  tp_youden <- mapply(function(r, id) {
    r$parameters$parameters_links_aa$tp_links.aa[id]
  },
  results,
  youden_max)
  
  # Positive predictive value retrieved related to Youden max #
  ppv_youden <- mapply(function(r, id) {
    r$parameters$parameters_links_aa$ppv_links.aa[id]
  },
  results,
  youden_max)
  
  # Negative predictive value related to Youden max #
  npv_youden <- mapply(function(r, id) {
    r$parameters$parameters_links_aa$npv_links.aa[id]
  },
  results,
  youden_max)
  
  # Minimal support related to Youden max #
  minimal_support_youden <- mapply(function(r, id) {
    r$parameters$minimal_support[id]
  },
  results,
  youden_max)
  
  # Percentage of entire chains reconstructed related to Youden max #
  global_chains_youden <- mapply(function(r, id) {
    r$parameters$parameters_chains_aa$global_chains_aa[id]
  },
  results,
  youden_max)
  
  ######
  ### Creation of result data.table ###
  ######
  output <- data.table(pi = ParametersEditing(pi,2,FALSE),
                       sigma = ParametersEditing(sigma,2,FALSE),
                       poisson_scale = ParametersEditing(poisson_scale,2,FALSE),
                       shannon_entropy = ParametersEditing(median_shannon_entropy,2,FALSE),
                       tp_links.consensus = ParametersEditing(tp_links.consensus,1,FALSE),
                       se_links.consensus = ParametersEditing(se_links.consensus,4,TRUE),
                       global_chains_consensus = ParametersEditing(global_chains_consensus,4,TRUE),
                       se_links.aa = ParametersEditing(se_links.aa,4,TRUE),
                       sp_links.aa = ParametersEditing(sp_links.aa,4,TRUE),
                       tp_links.aa = ParametersEditing(tp_links.aa,4,FALSE),
                       ppv_links.aa = ParametersEditing(ppv_links.aa,4,TRUE),
                       npv_links.aa = ParametersEditing(npv_links.aa,4,TRUE),
                       global_chains_aa = ParametersEditing(global_chains_aa,4,TRUE),
                       se_links.aa_max = ParametersEditing(se_links.aa_max,4,TRUE),
                       minimal_support.aa_related_se = ParametersEditing(minimal_support.aa_related_se,4,FALSE),
                       tp_links.aa_related_se = ParametersEditing(tp_links.aa_related_se,4,FALSE),
                       sp_links.aa_related_se = ParametersEditing(sp_links.aa_related_se,4,TRUE),
                       ppv_links.aa_related_se = ParametersEditing(ppv_links.aa_related_se,4,TRUE),
                       npv_links.aa_related_se = ParametersEditing(npv_links.aa_related_se,4,TRUE),
                       global_chains_aa_related_se = ParametersEditing(global_chains_aa_related_se,4,TRUE),
                       youden_index_max = ParametersEditing(youden_max,0,FALSE),
                       minimal_support_youden = ParametersEditing(minimal_support_youden,4,FALSE),
                       se_youden = ParametersEditing(se_youden,4,TRUE),
                       sp_youden = ParametersEditing(sp_youden,4,TRUE),
                       tp_youden = ParametersEditing(tp_youden,2,FALSE),
                       ppv_youden = ParametersEditing(ppv_youden,4,TRUE),
                       npv_youden = ParametersEditing(npv_youden,4,TRUE),
                       global_chains_youden = ParametersEditing(global_chains_youden,4,TRUE)
  )
  
  ###########################
  #### Plot of ROC curve ####
  ###########################
  # Preparation of data #
  data_ROC <- rbindlist(lapply(seq_len(length(min.support)), 
                               FUN = ParametersAccordingSupport, 
                               results = results, 
                               min.support = min.support))
  
  
  
  data_ROC[, min.support_label := as.character(round(min.support, 4))]
  data_ROC[, no_line := .I]
  data_ROC[!no_line %in% seq(1,.N,10),
           min.support_label := ""]
  
  FigPlot <- ggplot(data = data_ROC, aes(x = 1-sp_50, y = se_50)) + 
    geom_rect(aes(xmin = 1-sp_97.5, xmax = 1-sp_2.5,
                  ymin = se_2.5, ymax = se_97.5), fill = '#FF4040') +
    geom_line() +
    geom_point() +
    theme_minimal() +
    xlab("1-Specificity") +
    ylab("Sensitivity") +
    ylim(c(0,1)) + 
    xlim(c(0,1)) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                 colour = "red", linetype = 2, size = 0.75) +
    geom_segment(aes(x = 1 - data_ROC[no_line == median(youden_max), sp_50] + 0.015, 
                     xend = 1 - data_ROC[no_line == median(youden_max), sp_50] - 0.015,
                     y = data_ROC[no_line == median(youden_max), se_50] - 0.02, 
                     yend = data_ROC[no_line == median(youden_max), se_50] + 0.02), 
                 colour = "#228B22", linetype = 1, size = 2) +
    geom_text_repel(label = data_ROC$min.support_label,
                    point.padding = 1)
    #                 force = 3, 
    #                 size = 7) +
    # theme(axis.title.x = element_text(size = 20),
    #       axis.title.y = element_text(size = 20),
    #       axis.text.x = element_text(size = 15),
    #       axis.text.y = element_text(size = 15),
    #       legend.text = element_text(size = 15))
  
  ## Complete the output ##
  list_output <- list(parameters = output,
                      ROC = FigPlot)
  
  ######################################
  #### Study of not retrieved links ####
  ######################################
  if(NotRetrieved){
    links_not_retrieved <- rbindlist(lapply(results, 
                                            function(r) {
                                              merge(r$results_mcmc$res_aa,
                                                    chains_detect100_bind[,.(ids_from,ids_to,new)],
                                                    by.x = c("to", "from"),
                                                    by.y = c("ids_to","ids_from"),
                                                    all.y = TRUE)[is.na(support)]
                                            }    
    )
    )
    
    links_not_retrieved[,from_to := paste(from, to, sep = "_")]
    from_to.dt <- as.data.table(table(links_not_retrieved$from_to, links_not_retrieved$new))
    setnames(from_to.dt, c("V1","V2"), c("from_to","new_link"))         
    from_to.dt <- from_to.dt[order(N, decreasing = TRUE)]
    from_to.dt <- from_to.dt[N != 0]
    
    ## Complete the output ##
    list_output$not_retrieved_links <- from_to.dt[, scenario := index]
  }
  
  ################################################
  #### Study of kappa (number of generations) ####
  ################################################
  if(StudyKappa){
    kappa_new_links <- rbindlist(lapply(results, 
                                        function(r) {
                                          merge(r$results_mcmc$res_aa,
                                                new_chains_added[,.(ids_from,ids_to)],
                                                by.x = c("to", "from"),
                                                by.y = c("ids_to","ids_from"))
                                        }    
    )
    )
    
    kappa_new_links[,from_to := paste(from, to, sep = "_")]
    kappa_new_links.dt <- as.data.table(table(kappa_new_links$from_to, kappa_new_links$kappa))
    setnames(kappa_new_links.dt, c("V1","V2"), c("from_to","kappa"))         
    kappa_new_links.dt <- kappa_new_links.dt[order(N, decreasing = TRUE)]
    kappa_new_links.dt <- kappa_new_links.dt[N != 0]
    
    ## Complete the output ##
    list_output$kappa_new_links <- kappa_new_links.dt[, scenario := index]
  }
  
  ##################################################
  #### Deletion of results to free memory space ####
  ##################################################
  rm(results)
  gc()
  
  return(list_output)
}