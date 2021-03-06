#' Process input data for outbreaker
#'
#' This function performs various checks on input data given to outbreaker.  It
#' takes a list of named items as input, performs various checks, set defaults
#' where arguments are missing, and return a correct list of data input. If no
#' input is given, it returns the default settings.
#'
#' Acceptables arguments for ... are:
#' \describe{
#' \item{dates}{dates a vector indicating the collection dates, provided either as
#' integer numbers or in a usual date format such as \code{Date} or
#' \code{POSIXct} format. By convention, zero will indicate the oldest date. If
#' the vector is named, the vector names will be used for matching cases to
#' contact tracing data and labelled DNA sequences.}
#' 
#' \item{dna}{the DNA sequences in \code{DNAbin} format (see
#' \code{\link[ape]{read.dna}} in the ape package); this can be imported from a
#' fasta file (extension .fa, .fas, or .fasta) using \code{adegenet}'s function
#' \link[adegenet]{fasta2DNAbin}.}
#'
#' \item{ctd}{the contact tracing data provided as a matrix/dataframe of two
#' columns, indicating a reported contact between the two individuals whose ids
#' are provided in a given row of the data, or an epicontacts object. In the case
#' of the latter, linelist IDs will be used for matching dates and DNA
#' sequences}
#'
#' \item{w_dens}{a vector of numeric values indicating the generation time
#' distribution, reflecting the infectious potential of a case t = 1, 2, ...
#' time steps after infection. By convention, it is assumed that
#' newly infected patients cannot see new infections on the same time step. If not
#' standardized, this distribution is rescaled to sum to 1.}
#'
#' \item{f_dens}{similar to \code{w_dens}, except that this is the distribution
#' of the colonization time, i_e. time interval during which the pathogen can
#' be sampled from the patient.}}
#' 
#' @param ... a list of data items to be processed (see description)
#'
#' @param data optionally, an existing list of data item as returned by \code{outbreaker_data}.
#'
#' @author Thibaut Jombart (\email{thibautjombart@@gmail.com})
#'
#' @export
#'
#' @examples
#'
#' x <- fake_outbreak
#' outbreaker_data(dates = x$sample, dna = x$dna, w_dens = x$w)
#'
outbreaker_data <- function(..., data = list(...)) {

  ## SET DEFAULTS ##
  defaults <- list(dates = NULL,
                   w_dens = NULL,
                   f_dens = NULL,
                   dna = NULL,
                   ctd = NULL,
                   N = 0L,
                   L = 0L,
                   D = NULL,
                   max_range = NA,
                   can_be_ances = NULL,
                   log_w_dens = NULL,
                   log_f_dens = NULL,
                   contacts = NULL,
                   C_combn = NULL,
                   C_nrow = NULL,
                   ids = NULL,
                   has_dna = logical(0),
                   id_in_dna = integer(0),
                   id_in_hosp_matrix = integer(0),
                   hosp_matrix = NULL,
                   n_cases = NULL)

  ## MODIFY DATA WITH ARGUMENTS ##
  data <- modify_defaults(defaults, data, FALSE)

  ## Set up case ids
  if(is.null(data$ids)) {
    if(!is.null(names(data$dates))) {
      data$ids <- names(data$dates)
    } else if(!is.null(data$ctd) & inherits(data$ctd, "epicontacts")){
      data$ids <- as.character(data$ctd$linelist$id)
    } else {
      data$ids <- as.character(seq_along(data$dates))
    }
  }

  ## CHECK DATA ##
  ## CHECK DATES
  if (!is.null(data$dates)) {
    if (inherits(data$dates, "Date")) {
      data$dates <- data$dates-min(data$dates)
    }
    if (inherits(data$dates, "POSIXct")) {
      data$dates <- difftime(data$dates, min(data$dates), units="days")
    }
    data$dates <- as.integer(round(data$dates))
    data$N <- length(data$dates)
    data$max_range <- diff(range(data$dates))
    ## get temporal ordering constraint:
    ## canBeAnces[i,j] is 'i' can be ancestor of 'j'
    data$can_be_ances <- outer(data$dates,
                               data$dates,
                               FUN="<") # strict < is needed as we impose w(0)=0
    diag(data$can_be_ances) <- FALSE
  }

  ## CHECK W_DENS
  if (!is.null(data$w_dens)) {
    if (any(data$w_dens<0)) {
      stop("w_dens has negative entries (these should be probabilities!)")
    }

    if (any(!is.finite(data$w_dens))) {
      stop("non-finite values detected in w_dens")
    }


    ## Remove trailing zeroes to prevent starting with -Inf temporal loglike
    if(data$w_dens[length(data$w_dens)] < 1e-15) {
      final_index <- max(which(data$w_dens > 1e-15))
      data$w_dens <- data$w_dens[1:final_index]
    }

    ## add an exponential tail summing to 1e-4 to 'w'
    ## to cover the span of the outbreak
    ## (avoids starting with -Inf temporal loglike)
    if (length(data$w_dens) < data$max_range) {
      length_to_add <- (data$max_range-length(data$w_dens)) + 10 # +10 to be on the safe side
      val_to_add <- stats::dexp(seq_len(length_to_add), 1)
      val_to_add <- 1e-4*(val_to_add/sum(val_to_add))
      val_to_add[val_to_add == 0] <- min(val_to_add[val_to_add != 0]) # In case there is a huge data range
      data$w_dens <- c(data$w_dens, val_to_add)
    }

    ## standardize the mass function
    data$w_dens <- data$w_dens / sum(data$w_dens)
    data$log_w_dens <- matrix(log(data$w_dens), nrow = 1)
  }

  ## CHECK F_DENS
  if (!is.null(data$w_dens) && is.null(data$f_dens)) {
    data$f_dens <- data$w_dens
  }
  if (!is.null(data$f_dens)) {
    if (any(data$f_dens<0)) {
      stop("f_dens has negative entries (these should be probabilities!)")
    }

    if (any(!is.finite(data$f_dens))) {
      stop("non-finite values detected in f_dens")
    }

    data$f_dens <- data$f_dens / sum(data$f_dens)
    data$log_f_dens <- log(data$f_dens)
  }
  
  ## CHECK DNA

  if (!is.null(data$dna)) {
    if (!inherits(data$dna, "DNAbin")) stop("dna is not a DNAbin object.")
    if (!is.matrix(data$dna)) data$dna <- as.matrix(data$dna)

    ## get matrix of distances

    data$L <- ncol(data$dna) #  (genome length)
    data$D <- as.matrix(ape::dist.dna(data$dna, model="N")) # distance matrix
    storage.mode(data$D) <- "integer" # essential for C/C++ interface

    ## get matching between sequences and cases

    if (is.null(rownames(data$dna))) {
      if (nrow(data$dna) != data$N) {
        msg <- sprintf(paste("numbers of sequences and cases differ (%d vs %d):",
                             "please label sequences"),
                       nrow(data$dna), data$N)
        stop(msg)
      }

      ## These need to be indices
      rownames(data$D) <- colnames(data$D) <- seq_len(data$N)

      ## These need to match dates/ctd ids
      rownames(data$dna) <- data$ids
    }
    data$id_in_dna <- match(data$ids, rownames(data$dna))
    if(any(is.na(match(rownames(data$dna), data$ids)))) {
      stop("DNA sequence labels don't match case ids")
    }

  } else {
    data$L <- 0L
    data$D <- matrix(integer(0), ncol = 0, nrow = 0)
    data$id_in_dna <- rep(NA_integer_, data$N)
  }
  data$has_dna <- !is.na(data$id_in_dna)


  ## CHECK CTD

  if (!is.null(data$ctd)) {
    ctd <- data$ctd
    if (inherits(ctd, c("matrix", "data.frame"))) {
      ## prevent factor -> integer conversion
      ctd <- apply(ctd, 2, as.character)
      if (!is.matrix(ctd)) {
        ctd <- as.matrix(ctd)
      }
      if(ncol(ctd) != 2) {
        stop("ctd must contain two columns")
      }
    } else if(inherits(ctd, "epicontacts")) {
      ## prevent factor -> integer conversion
      ctd <- apply(ctd$contacts[c("from", "to")], 2, as.character)
    } else {
      stop("ctd is not a matrix, data.frame or epicontacts object")
    }

    unq <- unique(as.vector(ctd[,1:2]))
    not_found <- unq[!unq %in% data$ids]
    if (length(not_found) != 0) {
      not_found <- sort(unique(not_found))
      stop(paste("Individual(s)", paste(not_found, collapse = ", "),
                 "in ctd are unknown cases (idx < 1 or > N")
           )
    }
    contacts <- matrix(0, data$N, data$N)
    mtch_1 <- match(ctd[,1], data$ids)
    mtch_2 <- match(ctd[,2], data$ids)
    contacts[cbind(mtch_2, mtch_1)] <- 1
    
    data$contacts <- contacts
    data$C_combn <- data$N*(data$N - 1)
    data$C_nrow <- nrow(ctd)
  } else {
    data$contacts <- matrix(integer(0), ncol = 0, nrow = 0)
  }

  ## CHECK HOSP_MATRIX
  if (!is.null(data$hosp_matrix) && sum(dim(data$hosp_matrix))>0) {
    if (!inherits(data$hosp_matrix, c("data.frame", "matrix"))) {
      stop("hosp_matrix is not a matrix or a data.frame")
    }
    if (!is.matrix(data$hosp_matrix)) {
      data$hosp_matrix <- as.matrix(data$hosp_matrix)
    }
    dims <- dim(data$hosp_matrix)
    if (!dims[1] == dims[2]) {
      stop("hosp_matrix must be a square matrix")
    }
    ## Check if colnames and rownames are identical
    col_names = colnames(data$hosp_matrix)
    row_names = rownames(data$hosp_matrix)
    if (is.null(col_names) | is.null(row_names)) {
      stop("hosp_matrix must have column and row names")
    }
    extra_cols = setdiff(col_names, row_names)
    extra_rows = setdiff(row_names, col_names)
    if (length(extra_cols)) {
      msg_cols = paste0("Column names ",
                        paste0(extra_cols, collapse = ', '),
                        " are missing in row names. ")
    } else {
      msg_cols = character(0)
    }
    if (length(extra_rows)) {
      msg_rows = paste0("Row names ",
                        paste0(extra_rows, collapse = ', '),
                        " are missing in column names.")
    } else {
      msg_rows = character(0)
    }
    if (length(msg_cols) | length(msg_rows)) {
      stop(msg_cols, msg_rows)
    }
    ## Check for negative values
    if (any(data$hosp_matrix < 0)) {
      stop("hosp_matrix should only contain-non negative values")
    }
    if (any(is.na(data$hosp_matrix))) {
      stop("hosp_matrix has missing values")
    }
    data$id_in_hosp_matrix <- match(data$ids, rownames(data$hosp_matrix))
    if(any(is.na(data$id_in_hosp_matrix))) {
      stop("Hospital ids missing in hosp_matrix")
    }
  } else {
    data$hosp_matrix <- matrix(integer(0), ncol = 0, nrow = 0)    
  }

  ## CHECK N_CASES
  if (!is.null(data$n_cases)) {
    if (!is.numeric(data$n_cases)) {
      stop("n_cases must be an integer vector")
    }
  }
  
  ## output is a list of checked data
  return(data)

}
