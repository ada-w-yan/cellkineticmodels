library(dplyr)
library(ggplot2)
library(tidyr)
library(deSolve)
library(latex2exp)
# load R package with model fitting functions
devtools::load_all("~/git_repos/cellkineticmodels/") # change this to the directory of the git repository
# working directory should be scripts/ relative to the directory of the git repository
# i.e. the folder this vignette is in
options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)
set.seed(1)
dir_name <- "/data/projects/punim0053/ayan/git_repos/cellkineticmodels/scripts/two_compartment_vs_kh_files/"
# dir_name <- "two_compartment_vs_kh_files/"
if(!dir.exists(dir_name)) dir.create(dir_name)
quantiles_filename <- paste0(dir_name, "bar_p_quantiles.rds")
loo_filename <- paste0(dir_name, "loo.rds")
n_compartments <- c(1, 1.5, 2)
make_plots <- TRUE

id_pars <- list(frac = 0.032, delta = 0.064, b_w = 4.18)

model <- function(t, x, pars) {
  if(t < 0) {
    return(list(0))
  }
  with(as.list(c(x, pars)), {
    calc_U <- function(t) {
      if(t <= label_end) {
        U <- frac * (1 - exp(-delta * t))
      } else {
        t_shift <- t - label_end
        U <- calc_U(label_end) * exp(-delta * t_shift)
      }
      U
    }
    
    U <- calc_U(t)
    dLdt <- p * b_w * U - p * L
    
    return(list(dLdt))
  })
}

solve_label_one_pop <- function(pars) {
  yini <- c(L = 0)
  t1 <- c(seq(0, 150))
  out <- ode(yini, t1, model, pars) %>% as.data.frame %>% as_tibble
  out
}

solve_label_multi_pop <- function(pars) {
  solve_label_inner <- function(pop_idx) {
    pars$p <- pars$p[pop_idx]
    solve_label_one_pop(pars) %>%
      mutate(pop_idx = as.character(pop_idx))
  }
  out <- lapply(seq_along(pars$p), solve_label_inner) %>%
    do.call(rbind, .)
  mean_label <- out %>%
    group_by(time) %>%
    summarise(L = sum(pars[["frac_pop"]] * L)) %>%
    ungroup() %>%
    mutate(pop_idx = "mean")
  out <- bind_rows(out, mean_label)
  out
}

p1_vec <- c(0.720, 0.360, .24, .12, 0.036, 0.018, 0.0072)
p2_vec <- c(0.0036, 0.016, 0.0072, 0.0108) %>% sort
alpha1_vec <- seq(.1, .9, by = .2)

sim_pars_grid <- expand_grid(p1 = p1_vec, p2 = p2_vec, alpha1 = alpha1_vec) %>%
  filter(p1 > p2) %>%
  mutate(idx = seq_len(n()))
sim_data_wrapper <- function(p1, p2, alpha1, eval = TRUE) {
  sim_filename <- paste0(dir_name, "sim_data_", num2str(p1), "_",
                         num2str(p2), "_", num2str(alpha1), "__body_water.rds")
  if(!eval) {
    return(sim_filename)
  }
  
  sim_pars <- list(p = c(p1, p2), frac_pop = c(alpha1, 1 - alpha1), label_end = 49)
  sim_pars <- c(sim_pars, id_pars, list(sigma_L = 0.005))
  
  t_vec <- seq_len(14) * 7
  sol <- solve_label_multi_pop(sim_pars)
  
  sim_data <- sol %>%
    filter(pop_idx == "mean" & time %in% t_vec) %>%
    mutate(Lhat = vnapply(L, function(x) rnorm(1, x, sim_pars$sigma_L)))
  saveRDS(list(pars = sim_pars, sol = sol, sim_data = sim_data), sim_filename)
  sim_data
}

p1_vec <- p1_vec %>% sort

fit_wrapper <- function(n_compartments, p1, p2, alpha1, prior_width, eval = TRUE) {
  sim_filename <- sim_data_wrapper(p1, p2, alpha1, eval = FALSE)
  
  if(n_compartments >= 2) {
    model_name <- paste0(n_compartments, "_compartment")
  } else if(n_compartments == 1.5){
    model_name <- "kinetic_heterogeneity"
  } else {
    model_name <- "single"
  }
  
  fit_filename <- sim_filename %>% gsub("sim_data", paste0("fit_", model_name, "_prior_width_", prior_width), .)
  
  if(file.exists(fit_filename) | !eval) {
    return(fit_filename)
  }
  
  sim_data <- sim_filename %>%
    readRDS
  if(n_compartments >= 2) {
    fit_data_list <- c(list(T = nrow(sim_data$sim_data),
                            ts = sim_data$sim_data$time,
                            L = sim_data$sim_data$Lhat,
                            sigma = sim_data$pars$sigma_L,
                            n_pop = n_compartments,
                            perfect_labelling = 0,
                            prior_width = prior_width),
                       sim_data$pars[c("label_end", "delta", "frac", "b_w")])
    # model <- "lymphocyte_model_fixed_U_n_compartments_variable_prior"
  } else if(n_compartments == 1.5) {
    fit_data_list <- c(list(T = nrow(sim_data$sim_data),
                            ts = sim_data$sim_data$time,
                            L = sim_data$sim_data$Lhat,
                            sigma = sim_data$pars$sigma_L,
                            prior_width = prior_width),
                       sim_data$pars[c("label_end", "delta", "frac", "b_w")])
    
    model <- "lymphocyte_model_fixed_U_b_w_variable_prior"
  } else {
    fit_data_list <- c(list(T = nrow(sim_data$sim_data),
                            ts = sim_data$sim_data$time,
                            L = sim_data$sim_data$Lhat,
                            sigma = sim_data$pars$sigma_L,
                            n_pop = n_compartments,
                            perfect_labelling = 0,
                            prior_width = prior_width),
                       sim_data$pars[c("label_end", "delta", "frac", "b_w")])
    model <- "lymphocyte_model_fixed_U_homogeneous_variable_prior"
  }
  
  fit <- fit_stan_defaults(model, fit_data_list, adapt_delta = 0.99)
  saveRDS(fit, fit_filename)
  invisible(fit)
}

prior_width_vec <- c(2, 5, 10)
fit_grid <- expand_grid(n_compartments = n_compartments, 
                        sim_pars_grid, 
                        prior_width = prior_width_vec, eval = TRUE) %>%
  select(-idx)

model <- "lymphocyte_model_fixed_U_n_compartments_variable_prior"
model <- rstan::stan_model(get_model_filename(model))
apply_named_args(fit_grid[seq(127, nrow(fit_grid)),], 1, fit_wrapper)