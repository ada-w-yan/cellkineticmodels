#' perform leave-one-out cross validation and comparison for multiple models
#'
#' @param fit_list a list of outputs from rstan::stan
#' @return a list with the elements
#' loo_list: list of outputs from loo::loo for each model fit
#' compare: output of loo::loo_compare across the model fits
#' @importFrom loo extract_log_lik relative_eff loo loo_compare
#' @export
calc_and_compare_loos <- function(fit_list) {
  calc_loo <- function(log_lik) {
    if(inherits(log_lik, "list") && "fit" %in% names(log_lik)) {
      log_lik <- extract_log_lik(log_lik$fit, merge_chains = FALSE)
    } else if(inherits(log_lik, "stanfit")) {
      log_lik <- extract_log_lik(log_lik, merge_chains = FALSE)
    }
    r_eff <- relative_eff(exp(log_lik))
    loo_results <- loo(log_lik, r_eff = r_eff)
    loo_results
  }
  loo_list <- lapply(fit_list, calc_loo)

  compare <- loo_compare(loo_list)
  # compare <- do.call(loo_compare, loo_list)
  return(list(loo_list = loo_list,
         compare = compare))
}

#' contatenate log likelihoods for fits of the same model across different data sets
#'
#' @param fit_filenames a vector of length M: filenames containing model fits (output of rstan::stan)
#' @return an array of dimension c(N_samples, n_chains, M) concatenating the log likelihood
#' arrays returned from loo::extract_log_lik for each filename
#' @importFrom loo extract_log_lik
#' @importFrom dplyr %>%
#' @importFrom abind abind
#' @export
read_and_concat_log_lik <- function(fit_filenames) {
  read_extract_log_lik <- function(fit_filename) {
    readRDS(fit_filename)$fit %>%
      extract_log_lik(merge_chains = FALSE)
  }
  log_liks <- lapply(fit_filenames, read_extract_log_lik)
  nrow_target <- vnapply(log_liks, function(x) dim(x)[1]) %>% min
  ncols <- vnapply(log_liks, function(x) dim(x)[2])
  stopifnot(all(ncols == ncols[1]))
  log_liks <- lapply(log_liks, function(x) x[seq_len(nrow_target),,])
  concat_log_lik <- do.call(abind, list(log_liks, along = 3))
  concat_log_lik
}

#' short summary of preferred model order
#'
#' @param compare_output output of loo::loo_compare
#' @return a string of form "model X > model Y > model Z", listing models
#' from most to least preferred
#' @importFrom dplyr %>%
#' @export
preferred_model_order <- function(compare_output) {
  model_names <- rownames(compare_output)
  compare_output <- compare_output[-1,, drop = FALSE]
  n_ses <- -compare_output[,"elpd_diff"]/compare_output[,"se_diff"]
  for(i in seq_along(n_ses)) {
    if(n_ses[i] > 4) {
      model_names[i + 1] <- paste0(model_names[i + 1], "*")
    }
  }
  model_names %>%
    list %>%
    c(list(collapse = " > ")) %>%
    do.call(paste0, .)
}

#' given names of .rds files outputted by fitting models to data using RStan,
#' calculate the median for a given parameter across all of them
#'
#' @param filenames character vector
#' @param par_name string: name of parameter to take the median of
#' @return numeric scalar: population median
#' @export
calc_population_median_across_individual_posteriors <- function(filenames, par_name) {
  extract_pars <- function(filename) {
    fit <- readRDS(filename)$fit %>%
      as.matrix
    fit[,par_name]
  }
  pars <- lapply(filenames, extract_pars)
  n_samples <- vnapply(pars, length)
  # make sure each fit has the same number of iterations to estimate unbiased median
  if(any(n_samples != n_samples[1])) {
    min_samples <- min(n_samples)
    pars <- lapply(pars, sample, size = min_samples)
  }
  pars %>%
    unlist %>%
    median
}

#' calculate median, 2.5th and 97.5th percentiles
#'
#' @param x numeric vector
#' @return named numeric vector of length 3: median, 2.5th and 97.5th percentiles
#' @export
calc_med_95 <- function(x) {
  quantile(x, probs = c(.025, .5, .975), na.rm = TRUE)
}

#' get number of divergences in each chain
#'
#' @param fit stanfit object
#' @return data frame with two columns, and number of rows equal to the number of chains.
#' column 1: divergences: the number of divergences in each chain,
#' column 2: samples:  the total number of samples
#' @export
get_divergence <- function(fit) {
  sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
  divergences <- vnapply(sampler_params, function(x) sum(x[,"divergent__"]))
  samples <- vnapply(sampler_params, nrow)
  data.frame(divergences = divergences,
             samples = samples)
}

#' get maximum value of Rhat across model parameters and generated quantities
#'
#' @param fit stanfit object
#' @return scalar: maximum value of Rhat
#' @export
get_max_Rhat <- function(fit) {
  max(rstan::summary(fit)$summary[,"Rhat"], na.rm = TRUE)
}

#' compute 95% CIs for ratios between lymphocyte parameters between cell populations
#'
#' @param fit output of fit_all_stan
#' @return named list with three elements: pb_w, dstar and delay.
#' Each element is a named list with two elements: quantiles and p_value.
#' quantiles is a matrix with 3 rows and p(p+1)/2 columns, where p is the number of cell
#' populations for the individual in question.  Eech column shows the 2.5th percentile,
#' 50th percentile and 97.5th percentile for the ratio between parameter values for
#' pairs of cell populations.
#' p_value is a vector of length p(p+1)/2.  Each element is 2 * min(q, 1 - q), where q is the
#' proportion of draws from the posterior whose ratios are less than 1, for a particular
#' pair of cell populations.
#' @export
within_individual_comparison <- function(fit) {
  par_names <- c(get_par_names()$L, "p")
  pop <- fit$data$pop
  pairs1 <- combn(seq_len(pop), 2)
  paired_comparison <- function(par_name) {
    if(par_name == "p") {
      pb_w_indexed <- vcapply(seq_len(pop), index_par_name, par_name = "pb_w")
      pb_w <- extract_fit(fit$fit, par_names = pb_w_indexed, N_samples = 0)
      b_w <- extract_fit(fit$fit, par_names = "b_w", N_samples = 0, drop = TRUE)
      pars <- pb_w / b_w
    } else {
      par_names_indexed <- vcapply(seq_len(pop), index_par_name, par_name = par_name)
      pars <- extract_fit(fit$fit, par_names = par_names_indexed, N_samples = 0)
    }
    ratios <- matrix(0, ncol = ncol(pairs1), nrow = nrow(pars))
    for(i in seq_len(ncol(pairs1))) {
      ratios[,i] <- pars[,pairs1[2,i]] / pars[,pairs1[1,i]]
    }
    quantiles <- apply(ratios, 2, calc_med_95)
    p_value <- apply(ratios, 2, function(x) 2 * min(mean(x > 1), mean(x < 1)))
    list(quantiles = quantiles, p_value = p_value)
  }
  comparisons <- lapply(par_names, paired_comparison)
  names(comparisons) <- par_names
  comparisons
}

#' compute 95% CIs for ratios between lymphocyte parameters between cell populations for a given id,
#' when inference was done across ids using a hierarchical model
#'
#' @param fit output of fit_all_stan
#' @param n integer. index of id in fit to compute for
#' @return named list with three elements: pb_w, dstar and delay.
#' Each element is a named list with two elements: quantiles and p_value.
#' quantiles is a matrix with 3 rows and p(p+1)/2 columns, where p is the number of cell
#' populations for the individual in question.  Eech column shows the 2.5th percentile,
#' 50th percentile and 97.5th percentile for the ratio between parameter values for
#' pairs of cell populations.
#' p_value is a vector of length p(p+1)/2.  Each element is 2 * min(q, 1 - q), where q is the
#' proportion of draws from the posterior whose ratios are less than 1, for a particular
#' pair of cell populations.
#' @export
within_individual_comparison_from_multi_id <- function(fit, n) {
  par_names <- c(get_par_names()$L, "p")
  pop <- fit$data$pop[n]
  pairs1 <- combn(seq_len(pop), 2)
  paired_comparison <- function(par_name) {
    if(par_name == "p") {
      pb_w_indexed <- vcapply(seq_len(pop), function(x) index_par_name("pb_w", c(n, x)))
      pb_w <- extract_fit(fit$fit, par_names = pb_w_indexed, N_samples = 0)
      b_w <- extract_fit(fit$fit, par_names = index_par_name("b_w", n), N_samples = 0, drop = TRUE)
      pars <- pb_w / b_w
    } else {
      par_names_indexed <- vcapply(seq_len(pop), function(x) index_par_name(par_name, c(n, x)))
      pars <- extract_fit(fit$fit, par_names = par_names_indexed, N_samples = 0)
    }
    ratios <- matrix(0, ncol = ncol(pairs1), nrow = nrow(pars))
    for(i in seq_len(ncol(pairs1))) {
      ratios[,i] <- pars[,pairs1[2,i]] / pars[,pairs1[1,i]]
    }
    quantiles <- apply(ratios, 2, calc_med_95)
    p_value <- apply(ratios, 2, function(x) 2 * min(mean(x > 1), mean(x < 1)))
    list(quantiles = quantiles, p_value = p_value)
  }
  comparisons <- lapply(par_names, paired_comparison)
  names(comparisons) <- par_names
  comparisons
}

#' compute 95% CIs for ratios between lymphocyte parameters for the same cell populations
#' between individuals
#'
#' @param cells_lic.status list of 2 character vectors: cell populations and licensing statuses for each individual.
#' @param fits 2 lists -- outputs of fit_all_stan for 2 individuals
#' @return named list with three elements: pb_w, dstar and delay.
#' Each element is a named list with two elements: quantiles and p_value.
#' quantiles is a matrix with 3 rows and p(p+1)/2 columns, where p is the number of cell
#' populations for the individual in question.  Eech column shows the 2.5th percentile,
#' 50th percentile and 97.5th percentile for the ratio between parameter values for
#' pairs of cell populations.
#' p_value is a vector of length p(p+1)/2.  Each element is 2 * min(q, 1 - q), where q is the
#' proportion of draws from the posterior whose ratios are less than 1, for a particular
#' pair of cell populations.
#' @export
between_individual_comparison <- function(cells_lic.status, fits) {
  stopifnot(length(fits) == 2)
  stopifnot(length(cells_lic.status) == 2)
  # find cell populations common beteween the two fits
  cells_lic.status_common <- do.call(intersect, unname(cells_lic.status))
  if(length(cells_lic.status_common) == 0) {
    return(NULL)
  }

  # extract parameters for the common cell populations only
  common_idx <- lapply(cells_lic.status, function(x) which (x %in% cells_lic.status_common))

  par_names <- c(get_par_names()$L, "p")

  get_pars <- function(common_idx, fit) {

    get_par <- function(par_name) {
      if(par_name == "p") {
        pb_w_indexed <- vcapply(common_idx, index_par_name, par_name = "pb_w")
        pb_w <- extract_fit(fit$fit, par_names = pb_w_indexed, N_samples = 0)
        b_w <- extract_fit(fit$fit, par_names = "b_w", N_samples = 0, drop = TRUE)
        pars <- pb_w / b_w
      } else {
        par_names_indexed <- vcapply(common_idx, index_par_name, par_name = par_name)
        pars <- extract_fit(fit$fit, par_names = par_names_indexed, N_samples = 0)
      }
      if(!is.matrix(pars)) {
        pars <- matrix(pars, ncol = 1)
      }
      pars
    }
    pars <- lapply(par_names, get_par)
    pars
  }
  pars <- Map(get_pars, common_idx, fits)

  # because fits are done separately, parameter ratios can be calculated by randomly sampling from
  # the two posteriors (unpaired) and dividing the two
  n_samples_compare <- 1e6
  N_samples_fit <- vnapply(pars, function(x) nrow(x[[1]]))
  samples <- vapply(N_samples_fit,
                    function(x) sample.int(x, size = n_samples_compare, replace = T),
                    numeric(n_samples_compare))

  comparison <- function(idx) {
    ratios <- pars[[2]][[idx]][samples[,2],,drop = FALSE] / pars[[1]][[idx]][samples[,1],,drop = FALSE]

    quantiles <- apply(ratios, 2, calc_med_95)
    p_value <- apply(ratios, 2, function(x) 2 * min(mean(x > 1), mean(x < 1)))
    if(!is.matrix(quantiles)) {
      quantiles <- matrix(quantiles, ncol = 1, dimnames = list(c("2.5%", "50%", "97.5%")))
    }
    colnames(quantiles) <- names(p_value) <- cells_lic.status_common
    list(quantiles = quantiles, p_value = p_value)
  }
  comparisons <- lapply(seq_along(par_names), comparison)
  names(comparisons) <- par_names
  comparisons
}

#' compute 95% CIs for ratios between lymphocyte parameters for the same cell populations
#' between individuals
#'
#' @param id_idx indices of 2 individuals in data set
#' @param cells_lic.status list of 2 character vectors: cell populations and licensing statuses for each individual.
#' @param fit output of fit_all_stan
#' @return tibble with columns:
#' each row is a cell population/parameter name combination.
#' 2.5%, 50%, 97.5%: double: 2.5th percentile, 50th percentile and 97.5th percentile for the ratio between parameter values for
#' pairs of cell populations.
#' p_value: double: Each element is 2 * min(q, 1 - q), where q is the
#' proportion of draws from the posterior whose ratios are less than 1, for a particular
#' pair of cell populations.
#' cells_lic.status_common: character: cell population and licensing status
#' par_name: character: one of pb_w, dstar, delay or p
#' id1: character: id of first participant
#' id2: character: id of second participant.  ratios are the parameter value for id2 divided by the
#' parameter value for id1.
#' @export
between_individual_comparison_from_multi_id <- function(id_idx, cells_lic.status, fit) {
  stopifnot(length(id_idx) == 2)
  stopifnot(length(cells_lic.status) == 2)
  # find cell populations common beteween the two fits
  cells_lic.status_common <- do.call(intersect, unname(cells_lic.status))
  if(length(cells_lic.status_common) == 0) {
    return(NULL)
  }

  # extract parameters for the common cell populations only
  common_idx <- lapply(cells_lic.status, function(x) which (x %in% cells_lic.status_common))

  par_names <- c(get_par_names()$L, "p")

  get_pars <- function(n, common_idx) {

    get_par <- function(par_name) {
      if(par_name == "p") {
        pb_w_indexed <- vcapply(common_idx, function(x) index_par_name("pb_w", c(n, x)))
        pb_w <- extract_fit(fit$fit, par_names = pb_w_indexed, N_samples = 0)
        b_w <- extract_fit(fit$fit, par_names = index_par_name("b_w", n), N_samples = 0, drop = TRUE)
        pars <- pb_w / b_w
      } else {
        par_names_indexed <- vcapply(common_idx, function(x) index_par_name(par_name, c(n, x)))
        pars <- extract_fit(fit$fit, par_names = par_names_indexed, N_samples = 0)
      }
      if(!is.matrix(pars)) {
        pars <- matrix(pars, ncol = 1)
      }
      pars
    }
    pars <- lapply(par_names, get_par) %>%
      do.call(cbind, .)
    pars
  }
  pars <- Map(get_pars, id_idx, common_idx)
  ratios <- pars[[2]] / pars[[1]]
  par_name_grid <- expand.grid(cells_lic.status_common = cells_lic.status_common,
                               par_name = par_names, stringsAsFactors = FALSE)
  colnames(ratios) <- paste(par_name_grid$par_names, par_name_grid$cells_lic.status_common)
  quantiles <- apply(ratios, 2, calc_med_95)
  p_value <- apply(ratios, 2, function(x) 2 * min(mean(x > 1), mean(x < 1)))
  if(!is.matrix(quantiles)) {
    quantiles <- matrix(quantiles, ncol = 1, dimnames = list(c("2.5%", "50%", "97.5%")))
  }
  quantiles <- quantiles %>%
    t %>%
    as_tibble %>%
    mutate(p_value = p_value) %>%
    bind_cols(par_name_grid)
  quantiles
}

#' check convergence diagnostics
#'
#' @param fit output from rstan::stan
#' @return a named numeric vector with the elements
#' total_iterations: total NUTS iterations
#' divergent_iterations: number of divergent iterations
#' exceed_max_treedepth_iterations: number of iterations which exceed the maximum treedepth
#' min_n_eff: minimum effective sample size across parameters
#' max_Rhat: maximum Rhat across parameters
#' @importFrom rstan get_divergent_iterations get_max_treedepth_iterations
#' @importFrom dplyr %>%
#' @export
check_diagnostics <- function(fit) {
  total_iterations <- length(rstan::get_divergent_iterations(fit))
  divergent_iterations <- sum(rstan::get_divergent_iterations(fit))
  exceed_max_treedepth_iterations <- sum(rstan::get_max_treedepth_iterations(fit))
  min_n_eff <- min(summary(fit)$summary[,"n_eff"], na.rm = TRUE)
  max_Rhat <- max(summary(fit)$summary[,"Rhat"], na.rm = TRUE)
  c(total_iterations = total_iterations,
    divergent_iterations = divergent_iterations,
    exceed_max_treedepth_iterations = exceed_max_treedepth_iterations,
    min_n_eff = min_n_eff,
    max_Rhat = max_Rhat)
}
