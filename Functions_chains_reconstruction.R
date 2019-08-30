library(plotly)
library(crosstalk)

##################################################
#### Functions for estimating generation time ####
##################################################
# # Estimation of generation time #
# generation_time <- function(x){
#   if(x[,.N] > 1){
#     x[, t_lag := c(NA, x[1:(.N-1), t])]
#     x[, w := t - t_lag]
#   }
#   else
#     x[, ":=" (t_lag = NA,
#                w = NA)]
# 
# }
# 
# # Finding the origin hospital #
# get_origin <- function(x)
# {
#   x[, from := c(NA, x[1:(.N-1), hospID])]
#   setnames(x, c("t","hospID"), c("time","to"))
# }

####################################
#### Function to compute priors ####
####################################
# Function to compute priors for alpha #
prior_ancestor <- function(x, imported, data_outbreaker, fakeMat){
  if(is.na(imported[x]))
    ancestor <- NA
  else{
    # Identification of potential ancestors #
    line <- which(rownames(data_outbreaker$hosp_matrix) %in% ids[data_outbreaker$can_be_ances[,x]])
    # Identification of considered establishment #
    col <- which(rownames(data_outbreaker$hosp_matrix)==ids[x])
    # Looking for the best ancestor according to transfer probability #
    max_value <- which(fakeMat[line, col] == max(fakeMat[line, col]))
    max_value <- sample(max_value,1)
    # Return best ancestor possible #
    ancestor <- which(data_outbreaker$can_be_ances[,x]==TRUE)[max_value]
  }
  return(ancestor)
}


###############################################################
#### Functions related to the creation of bayesian results ####
###############################################################
## Main function to compute bayesian results ##
ComputeBayesian <- function(outbreaker_data, n_iter_mcmc, ids, config){
  cat(paste0("--- Number of iterations of MCMC = ", n_iter_mcmc," ---"))
  
  # Reconstruction of chains using bayesian statistics #
  res <- outbreaker(outbreaker_data, config)
  
  cat("--- Preparation of results ---\n")
  # Comparison with real data #
  tmp <- CreateOutputBayesian(res, ids, burning = 0, 
                              init_alpha = config$init_alpha)
  
  return(list(res = res,
              res_consensus = tmp$consensus,
              res_aa = tmp$aa))
}

## Auxiliary functions ##
# Main function computing the results, according to the type of ancester selected, of bayesian results #
CreateOutputBayesian <- function(results_bayesian, ids, init_alpha, burning = 0){
  # Creating output with all ancesters (aa) #
  res.aa <- rbindlist(lapply(seq_len(length(ids)), 
                                    FUN = get_links, 
                                    x = as.data.table(results_bayesian)[step > burning]))
  
  # Creating output with most posterior ancestor (consensus) #
  res.consensus <- res.aa[order(to, support),]
  res.consensus <- res.consensus[, .SD[.N], by="to"]
  res.consensus <- res.consensus[, init_alpha := init_alpha]
  
  res.aa <- merge(res.aa,
                  res.consensus[, .(to, init_alpha)],
                  by = "to")
  
  # # Getting labels #
  res.aa[, from_label := sapply(as.numeric(res.aa$from), get_labels, labels=ids)]
  res.aa[, to_label := sapply(as.numeric(res.aa$to), get_labels, labels=ids)]
  res.consensus[, from_label := sapply(as.numeric(res.consensus$from), get_labels, labels=ids)]
  res.consensus[, to_label := sapply(as.numeric(res.consensus$to), get_labels, labels=ids)]
  
  return(list(consensus = res.consensus,
              aa = res.aa))
}

# Function to get links #
get_links <- function(x, name)
{
  selection <- x[, .SD,
                 .SDcols = c(paste0("alpha_",name), 
                             paste0("t_inf_",name), 
                             paste0("kappa_",name))]
  names(selection) <- c("alpha", "t_inf", "kappa")
  
  result <- selection[, .(support = .N,
                          time = median(t_inf, na.rm = TRUE),
                          kappa = median(kappa, na.rm = TRUE)),
                      by = "alpha"]
  result[, support := support / selection[,.N]]
  setnames(result, "alpha", "from")
  result[, to := name]
  
  return(result)
}

# Function to get labels #
get_labels <- function(x, labels = NULL) {
  if(!is.na(x)) x <- labels[x]
  return(x)
}

# # Auxiliary function making the comparison (links and chains) of bayesian results and real data #
# MergeData <- function(set, real_data){
#   compare.set <- merge(set,
#                        real_data[, .(from, to, 
#                                      time_inf = t, 
#                                      t_detect, 
#                                      chain)],
#                        by = c("from","to"),
#                        all = TRUE)
# 
#   return(compare.set)
# }

######################################################################################
#### Functions related to the estimation of parameters (sensitivity, specificity) ####
######################################################################################
## Main function to compute sensitivity and specificity ##
ComputeParameters <- function(results_bayesian, data_outbreaker, real_data, min.support,
                              burning = 0, init_alpha){
  ## Preparing data (merging with real data) ##
  cat("--- Preparing data ---")
  if(burning < 0) 
    stop("Burning period must be positive")
  else if(max(results_bayesian$res$step) < burning)
    stop("Burning period is greater than the total number of MCMC iterations")
  else if(burning == 0){ # No burning period
    merged_consensus <- merge(results_bayesian$res_consensus,
                              real_data[, .(from, to, 
                                            ids_from, ids_to,
                                            time_inf = t,
                                            t_detect,
                                            chain,
                                            length_chain)],
                              by.x = c("from","to"),
                              by.y = c("ids_from","ids_to"),
                              all = TRUE)
    
    merged_aa <- merge(results_bayesian$res_aa,
                       real_data[, .(from, to, 
                                     ids_from, ids_to,
                                     time_inf = t,
                                     t_detect,
                                     chain,
                                     length_chain)],
                       by.x = c("from","to"),
                       by.y = c("ids_from","ids_to"),
                       all = TRUE)
  }
  else{ # Considering a burning period
    # Calculating consensus and all ancestors with burning period #
    tmp <- CreateOutputBayesian(results_bayesian$res, ids, burning = burning,
                                init_alpha = init_alpha)
    merged_consensus <- merge(tmp$consensus,
                              real_data[, .(from, to, 
                                            ids_from, ids_to,
                                            time_inf = t,
                                            t_detect,
                                            chain,
                                            length_chain)],
                              by.x = c("from","to"),
                              by.y = c("ids_from","ids_to"),
                              all = TRUE)
    
    merged_aa <- merge(tmp$aa,
                       real_data[, .(from, to, 
                                     ids_from, ids_to,
                                     time_inf = t,
                                     t_detect,
                                     chain,
                                     length_chain)],
                       by.x = c("from","to"),
                       by.y = c("ids_from","ids_to"),
                       all = TRUE)
  } 
  
  ## We only keep the links in real_data ##
  merged_consensus <- merged_consensus[to %in% real_data[,ids_to]]
  merged_aa <- merged_aa[to %in% real_data[,ids_to]]
  
  ## Computing parameters ##
  cat("\n--- Compute parameters using consensus ancestor ---")
  # Consensus ancestor #
  parameters_links_consensus <- ComputeParametersLinks(data = merged_consensus, 
                                           min.support = NULL, 
                                           type = "consensus", 
                                           outbreaker_data = data_outbreaker)
  
  global_chains_consensus <- ComputeParametersChains_global(data = merged_consensus,
                                                                real_data = real_data)
  
  bylength_chains_consensus <- ComputeParametersChains_chainlength(data = merged_consensus,
                                                   real_data = real_data)
  
  cat("\n--- Compute parameters using all ancestors ---")
  # All ancestors #
  parameters_links_aa <- ComputeParametersLinks(data = merged_aa, 
                                          min.support = min.support, 
                                          type = "aa", 
                                          outbreaker_data = data_outbreaker)
  
  global_chains_aa <- sapply(min.support, 
                                  FUN = ComputeParametersChains_global,
                                  data = merged_aa,
                                  real_data = real_data)
  
  bylength_chains_aa <- sapply(min.support, 
                             FUN = ComputeParametersChains_chainlength,
                             data = merged_aa,
                             real_data = real_data)

  bylength_chains_aa.dt <- cbind(t(t(bylength_chains_aa[1,1]$length_chain)),
                              rbindlist(list(bylength_chains_aa[2,])))
  names(bylength_chains_aa.dt) <- c("length_chains",
                                    min.support)
  
  return(list(minimal_support = min.support,
              parameters_links_consensus = parameters_links_consensus,
              parameters_chains_consensus = list(global_chains_consensus = global_chains_consensus,
                                                 bylength_chains_consensus = bylength_chains_consensus),
              parameters_links_aa = parameters_links_aa,
              parameters_chains_aa = list(global_chains_aa = global_chains_aa,
                                          bylength_chains_aa = bylength_chains_aa.dt)))
}

## Auxiliary functions ##
# Auxiliary function to compute parameters regarding links #
ComputeParametersLinks <- function(data, min.support = NULL, type, outbreaker_data){
  if(is.null(min.support)){
    # True positive #
    assign(paste0("tp_links.",type),true_positive_links(data = data))
    
    # False positive #
    assign(paste0("fp_links.",type),false_positive_links(data = data))
    
    # False negative #
    assign(paste0("fn_links.",type),false_negative_links(data = data))
  }
  else{
    # True positive #
    assign(paste0("tp_links.",type),sapply(min.support, 
                                           FUN = true_positive_links,
                                           data = data))
    # False positive #
    assign(paste0("fp_links.",type),sapply(min.support, 
                                           FUN = false_positive_links,
                                           data = data))
    
    # False negative #
    assign(paste0("fn_links.",type),sapply(min.support, 
                                           FUN = false_negative_links,
                                           data = data))
  }
    
    # True negative #
    assign(paste0("tn_links.",type),true_negative_links(data = data,
                                                        tp = get(paste0("tp_links.",type)),
                                                        fp = get(paste0("fp_links.",type)),
                                                        fn = get(paste0("fn_links.",type))))
    
    # Sensitivity #
    assign(paste0("se_links.",type),sensitivity_links(tp = get(paste0("tp_links.",type)),
                                                      fn = get(paste0("fn_links.",type))))      
    
    # Specificity #
    assign(paste0("sp_links.",type),specificity_links(tn = get(paste0("tn_links.",type)),
                                                      fp = get(paste0("fp_links.",type))))
    
    # Positive predictive value #
    assign(paste0("ppv_links.",type),ppv_links(tp = get(paste0("tp_links.",type)),
                                                      fp = get(paste0("fp_links.",type))))      
    
    # Negative predictive value #
    assign(paste0("pnv_links.",type),pnv_links(tn = get(paste0("tn_links.",type)),
                                                      fn = get(paste0("fn_links.",type))))
    
    returnList <- list(get(paste0("tp_links.",type)),
                       get(paste0("fp_links.",type)),
                       get(paste0("fn_links.",type)),
                       get(paste0("tn_links.",type)),
                       get(paste0("se_links.",type)),
                       get(paste0("sp_links.",type)),
                       get(paste0("ppv_links.",type)),
                       get(paste0("pnv_links.",type)))
    
    names(returnList) <- c(paste0("tp_links.",type),
                           paste0("fp_links.",type),
                           paste0("fn_links.",type),
                           paste0("tn_links.",type),
                           paste0("se_links.",type),
                           paste0("sp_links.",type),
                           paste0("ppv_links.",type),
                           paste0("pnv_links.",type))
  
  return(returnList)
}

# Auxiliary function to compute the number of true positive links #
true_positive_links <- function(min_support = 0, data){
  return(data[support >= min_support & 
                !is.na(time_inf), .N])
} 

# Auxiliary function to compute the number of false positive links #
false_positive_links <- function(min_support = 0, data){
  return(data[support >= min_support & 
                is.na(time_inf), .N])
} 

# Auxiliary function to compute the number of false negative links #
false_negative_links <- function(min_support = 0, data){
  return(data[(is.na(time) |
                (!is.na(time_inf) & support < min_support)), .N])
} 

# Auxiliary function to compute the number of true negative links #
true_negative_links <- function(data, tp, fp, fn){
  result <- rep(data[,.N],
                length(tp))
  result <- result - tp - fp - fn
  return(result)
} 

# Auxiliary function to compute the sensitivity regarding links reconstruction #
sensitivity_links <- function(tp, fn){
  return(tp/(tp + fn))
}

# Auxiliary function to compute the specificity regarding links reconstruction #
specificity_links <- function(tn, fp){
  return(tn/(tn + fp))
}

# Auxiliary function to compute the predictive positive value regarding links reconstruction #
ppv_links <- function(tp, fp){
  return(tp/(tp + fp))
}

# Auxiliary function to compute the predictive negative value regarding links reconstruction #
pnv_links <- function(tn, fn){
  return(tn/(tn + fn))
}

# Auxiliary function to compute parameters regarding chains #
ComputeParametersChains_global <- function(data, min_support = NULL, real_data){
  # Deleting cycles in real data #
  # real_data_decycled <- unique(real_data[, .(from, to, chain)])
  
  if(is.null(min_support)){
    results.chains <- merge(data[!is.na(from) & !is.na(time), .N, by = "chain"],
                            real_data[!is.na(from), .N, by = "chain"],
                            by="chain",
                            all.y = TRUE)
    results.chains[, I_match := (N.x == N.y)]
    results.chains <- merge(results.chains,
                            unique(real_data[, 
                                             .(chain, length_chain)]),
                            by="chain")
    
    global_percentage <- results.chains[I_match == 1, .N]/
      results.chains[,.N]
  }
  else{
    # Percentage of corresponding chains #
    # tmp.data <- data[!is.na(from) & !is.na(time) & 
    #                    support >= min_support, ]
    # tmp.data <- unique(tmp.data[!is.na(chain), .(from, to, chain)])
    results.chains <- merge(data[!is.na(from) & !is.na(time) &
                                   support >= min_support, .N, by = "chain"],
                            real_data[!is.na(from), .N, by = "chain"],
                            by="chain",
                            all.y = TRUE)
    results.chains[, I_match := (N.x == N.y)]
    results.chains <- merge(results.chains,
                            unique(real_data[, 
                                             .(chain, length_chain)]),
                            by="chain")
    
    global_percentage <- results.chains[I_match == 1, .N]/
      results.chains[,.N]
  }
  return(global_percentage)
}

ComputeParametersChains_chainlength <- function(data, min_support = NULL, real_data){
  # Deleting cycles in real data #
  # real_data_decycled <- unique(real_data[, .(from, to, chain)])
  
  if(is.null(min_support)){
    results.chains <- merge(data[!is.na(from) & !is.na(time), .N, by = "chain"],
                            real_data[!is.na(from), .N, by = "chain"],
                            by="chain",
                            all.y = TRUE)
    results.chains[, I_match := (N.x == N.y)]
    results.chains <- merge(results.chains,
                            unique(real_data[, 
                                             .(chain, length_chain)]),
                            by="chain")
    
    chains_by_length <- merge(results.chains[, sum(I_match, na.rm = TRUE),
                                             by = "length_chain"],
                              results.chains[, .N, 
                                             by = "length_chain"],
                              by = "length_chain")
    chains_by_length[, percentage_reconstructed := V1 / N]
  }
  else{
    # Percentage of corresponding chains #
    # tmp.data <- data[!is.na(from) & !is.na(time) & 
    #                    support >= min_support, ]
    # tmp.data <- unique(tmp.data[!is.na(chain), .(from, to, chain)])
    results.chains <- merge(data[!is.na(from) & !is.na(time) &
                                   support >= min_support, .N, by = "chain"],
                            real_data[!is.na(from), .N, by = "chain"],
                            by="chain",
                            all.y = TRUE)
    results.chains[, I_match := (N.x == N.y)]
    results.chains <- merge(results.chains,
                            unique(real_data[, 
                                             .(chain, length_chain)]),
                            by="chain")
    
    chains_by_length <- merge(results.chains[, sum(I_match, na.rm = TRUE),
                                             by = "length_chain"],
                              results.chains[, .N, 
                                             by = "length_chain"],
                              by = "length_chain")
    chains_by_length[, percentage_reconstructed := V1 / N]
  }
  return(chains_by_length[,.(length_chain, percentage_reconstructed)])
}


###########################################
#### Functions to represent parameters ####
###########################################
PlotROC <- function(se, sp, type.plot){
  data <- data.table(x = 1 - sp, y = se)
  
  FigPlot<-ggplot(data = data, aes(x = x, y = y)) + 
    geom_line() +
    geom_point() +
    theme_minimal() +
    xlab("1-specificity") +
    ylab("Sensitivity") +
    ylim(c(0,1)) + 
    xlim(c(0,1)) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                 colour = "red", linetype = 2, size = 0.75)
  
  if(type.plot == "ggplot"){
    FigPlot
  }
  else if(type.plot == "plotly"){
    ggplotly(FigPlot)
  }
}

Plot_se_ppv <- function(se, ppv, type.plot){
  data <- data.table(x = se, y = ppv)
  
  FigPlot<-ggplot(data = data, aes(x = x, y = y)) + 
    geom_line() +
    geom_point() +
    theme_minimal() +
    xlab("Sensitivity") +
    ylab("Positive predictive value") +
    ylim(c(0,1)) + 
    xlim(c(0,1)) +
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=17,face="bold"))
  
  if(type.plot == "ggplot"){
    FigPlot
  }
  else if(type.plot == "plotly"){
    ggplotly(FigPlot)
  }
}

PlotChains <- function(parameters_chains, minimal_support, type.plot){
  data <- data.table(x = minimal_support, y = parameters_chains)
  
  FigPlot<-ggplot(data = data, aes(x = x, y = y*100)) + 
    geom_line() +
    geom_point() +
    theme_minimal() +
    xlab("Minimal support") +
    ylab("Percentage of reconstructed chains")
  
  if(type.plot == "ggplot"){
    FigPlot
  }
  else if(type.plot == "plotly"){
    ggplotly(FigPlot)
  }
} 







