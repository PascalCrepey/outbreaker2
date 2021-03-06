## #' Stores MCMC samples for outbreaker
## #'
## #' This function creates stores MCMC samples for outbreaker: augmented data and parameter states, likelihood, priors and posterior.
## #'
## #' @author Thibaut Jombart (\email{thibautjombart@@gmail.com})
## #'
## #' @inheritParams outbreaker
## #'
## #' @param param a list of output items as returned by \code{create_param}
## #'
## #' @param step an integer indicating the MCMC iteration being stored
## #'
## #' @export
## #'
outbreaker_mcmc_store <- function(param_current, param_store, data, config,
                                  likelihoods, priors, step) {
  ## UPDATE COUNTER
  counter <- param_store$counter <- param_store$counter + 1

  ## STORE STEP
  param_store$step[counter] <- step

  ## STORE LIKELIHOOD, PRIOR, POSTERIOR
  param_store$like[counter] <- cpp_ll_all(data, param_current, NULL, likelihoods)
  param_store$prior[counter] <- cpp_prior_all(param_current, config, priors)
  param_store$post[counter] <- param_store$like[counter] + param_store$prior[counter]

  ## PARAMETERS AND AUGMENTED DATA
  param_store$mu[counter] <- param_current$mu
  param_store$pi[counter] <- param_current$pi
  param_store$eps[counter] <- param_current$eps
  param_store$lambda[counter] <- param_current$lambda
  param_store$sigma[counter] <- param_current$sigma
  param_store$poisson_scale[counter] <- param_current$poisson_scale
  
  param_store$alpha[[counter]] <- param_current$alpha
  param_store$t_inf[[counter]] <- param_current$t_inf
  param_store$kappa[[counter]] <- param_current$kappa
  param_store$potential_colonised[[counter]] <- param_current$potential_colonised

  return(param_store)
}

