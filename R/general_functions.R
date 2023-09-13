#' get model filename
#' @param model_name string
#' @return string: model filename
#' @export
get_model_filename <- function(model_name) {
  git_repo_dirname <- "~/git_repos/cellkineticmodels/"
  paste0(git_repo_dirname, "model_stan/", model_name, ".stan")
}

#' short for \code{vapply(X, FUN, logical(1), ...)}
#' 
#' short for \code{vapply(X, FUN, logical(1), ...)}
#' 
#' @param X first input to vapply
#' @param FUN second input to vapply
#' @return output of \code{vapply(X, FUN, logical(1), ...)}: logical vector of 
#' same length as X
vlapply <- function(X, FUN, ...) {
  vapply(X, FUN, logical(1), ...)
}

#' short for \code{vapply(X, FUN, integer(1), ...)}
#' 
#' short for \code{vapply(X, FUN, integer(1), ...)}
#' 
#' @param X first input to vapply
#' @param FUN second input to vapply
#' @return output of \code{vapply(X, FUN, integer(1), ...)}: integer vector of 
#' same length as X
viapply <- function(X, FUN, ...) {
  vapply(X, FUN, integer(1), ...)
}

#' short for \code{vapply(X, FUN, numeric(1), ...)}
#' 
#' short for \code{vapply(X, FUN, numeric(1), ...)}
#' 
#' @param X first input to vapply
#' @param FUN second input to vapply
#' @return output of \code{vapply(X, FUN, numeric(1), ...)}: numeric vector of 
#' same length as X
vnapply <- function(X, FUN, ...) {
  vapply(X, FUN, numeric(1), ...)
}

#' short for \code{vapply(X, FUN, character(1), ...)}
#' 
#' short for \code{vapply(X, FUN, character(1), ...)}
#' 
#' @param X first input to vapply
#' @param FUN second input to vapply
#' @return output of \code{vapply(X, FUN, character(1), ...)}: character vector
vcapply <- function(X, FUN, ...) {
  vapply(X, FUN, character(1), ...)
}

#' a version of apply. if MARGIN = 1, the values in each column of X are passed to FUN
#' as named arguments according to the column name. if MARGIN = 2, the values in
#' each row of X are passed to FUN
#' as named arguments according to the row name
#' @param X see arguments for apply.  If MARGIN = 1, colnames(X) must match the named
#' arguments of FUN.
#' @param MARGIN see arguments for apply
#' @param FUN see arguments for apply.  Must have named arguments
#' @return see apply
apply_named_args <- function(X, MARGIN, FUN, keep_class = TRUE) {
  if(MARGIN == 2 | !keep_class) {
    return(apply(X, MARGIN, function(X) do.call(FUN, as.list(X))))
  }
  args <- formalArgs(FUN)
  X <- c(list(FUN), as.list(X[,args]))
  do.call(Map, X)
}

#' turn number(s) into string for filename, replacing decimal points with "point"
#' 
#' @param x numeric vector
#' @return character vector of same length as x
#' @export
num2str <- function(x) {
  gsub("\\.", "point", as.character(x))
}

#' concatenate strings for filename
#' 
#' @param strs list of characters
#' @param ext1 string: extension
#' @return string
#' @export
make_filename <- function(strs, ext1 = "rds") {
  vcapply(strs, as.character) %>%
  paste0(collapse = "_") %>%
    gsub(" ", "_", .) %>%
    paste(ext1, sep = ".")
}

#' subsample elements of vector with constant spacing
#' 
#' @param vec vector
#' @param length_out number of samples
#' @return vector of length length_out
#' @export
subsample <- function(vec, length_out) {
  vec[round(seq(1, length(vec), length.out = length_out))]
}

#' get the number of iterations post warm-up for each chain in a stanfit object
#' 
#' @param fit stanfit object
#' @return vector: iterations post warm-up for each chain
#' @export
get_n_iter <- function(fit) {
  vnapply(fit@stan_args, function(x) x$iter - x$warmup)
}

#' repeat elements of vector x, y times each element-wise, where y is a vector
#' 
#' @param x vector
#' @param y vector of same length as x
#' @return vector of length sum(y)
#' @importFrom dplyr %>%
#' @export
rep_each <- function(x, y) {
  Map(rep, x, y) %>%
    unlist %>%
    unname
}

#' linear interpolation between a set of points; return gradient of each line
#' 
#' @param x numeric vector: x-coords of points.  Must be strictly increasing.
#' @param y numeric vector: y-coords of points.  Must be same length as x.
#' @return numeric vector one shorter than x.
#' @export
calc_grad <- function(x, y) {
    stopifnot(length(y) == length(x))
    stopifnot(all(diff(x) > 0))
    diff(y) / diff(x)
}

#' evaluate a function for all cell populations and licensing statuses for a given id
#' 
#' @param f function to run.  First two arguments are cells and lic.status
#' @param id1 id to run for
#' @param ... other arguments to pass to f
#' @return list of outputs for all cell populations and licensing statuses
#' @export
run_for_id <- function(f, id1, ...) {
    id_data <- read.csv(paste0(get_input_dir(), "Labelling_Study_Participant_Lymphocyte_Data.csv")) %>%
        subset(id == id1 & days == 0)
    Map(function(x, y) f(x, y, ...), id_data$cells, id_data$lic.status)
}

#' evaluate a function for all cell populations and licensing statuses for given ids,
#' and label outputs with id, cell population and licensing status
#' 
#' @param f function to run.  First two arguments are cells and lic.status
#' @param ids ids to run for
#' @param ... other arguments to pass to f
#' @return data frame of outputs for ids, cell populations and licensing statuses
#' @export
run_ids_label <- function(f, ids) {
    run_and_label <- function(id, cells, lic.status) {
        res <- f(id, cells, lic.status) %>%
            as.data.frame
        res$cells <- cells
        res$lic.status <- lic.status
        res
    }
    wrapper <- function(id1) {
        id_data <- read.csv(paste0(get_input_dir(), "Labelling_Study_Participant_Lymphocyte_Data.csv")) %>%
            subset(id == id1 & days == 0)
        res <- Map(function(x, y) run_and_label(id1, x, y), id_data$cells, id_data$lic.status) %>%
            do.call(rbind, .)
        res$id <- id1
        res
    }
    lapply(ids, wrapper) %>%
        do.call(rbind, .)

}

#' Change fig.width and fig.height Dynamically Within an R Markdown Chunk
#' http://michaeljw.com/blog/post/subchunkify/
#' 
#' @param g a ggplot object
#' @param fig_height scalar: figure height
#' @param fig_width scalar: figure width
#' @return gpglot object with correct height and width
subchunkify <- function(g, fig_height=7, fig_width=5) {
  g_deparsed <- paste0(deparse(
    function() {g}
  ), collapse = '')
  
  sub_chunk <- paste0("
  `","``{r sub_chunk_", floor(runif(1) * 10000), ", fig.height=",
                      fig_height, ", fig.width=", fig_width, ", echo=FALSE}",
                      "\n(", 
                      g_deparsed
                      , ")()",
                      "\n`","``
  ")
  
  cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
}

#' adjust N_samples for predictions
#' 
#' @param N_samples original value for N_samples
#' @param chain_length length of NUTS chain
#' @return integer: chain_length if N_samples == 0, min(N_samples, chain_length) otherwise
#' @export
adjust_N_samples <- function(N_samples, chain_length) {
  if(N_samples == 0 || (N_samples > chain_length)) {
    if(N_samples > chain_length) {
      warning("number of samples requested exceeds length of chain, using length of chain instead")
    }
    N_samples <- chain_length
  }
  N_samples
}

#' append indices to parameter name
#' 
#' @param par_name string: parameter name
#' @param idx numeric vector
#' @return string.  e.g. if par_name == "abc" and idx == c(1,2), returns "abc[1,2]"
#' @export
index_par_name <- function(par_name, idx) {
  idx <- paste(idx, collapse = ",")
  paste0(par_name, "[", idx, "]")
}

#' make model predictions
#' 
#' @param model stan model object
#' @param data_list list of data and model parameters to feed to model
#' @return stanfit object
#' @export
make_pred <- function(model, data_list) {
  if(is.character(model)) {
    model <- rstan::stan_model(get_model_filename(model))
  }
  pred <- rstan::sampling(model,
                   data = data_list,
                   chains = 1, iter = 1,
                   algorithm = "Fixed_param")
  pred
}

#' extract model predictions as matrix
#' 
#' @param pred stanfit object
#' @param N_samples number of samples used to make model predictions
#' @param pred_name (optional) name of model prediction.  If not specified, return everything
#' @return matrix of model predictions
#' @importFrom dplyr %>%
#' @export
extract_pred <- function(pred, N_samples, pred_name) {
  pred <- as.matrix(pred)
  if(!missing(pred_name)) {
    pred <- pred[,grepl(pred_name, colnames(pred), fixed = TRUE)]
  } else {
    pred <- pred[1, -ncol(pred)]
  }
  pred <- matrix(pred, ncol = N_samples)
  pred
}

#' fit sta model with default parameters
#' 
#' @param model stan model object
#' @param data_list list of data to pass to model
#' @param init initialisation
#' @param adapt_delta scalar.  value of adapt_delta to use
#' @param evaluate_likelihood logical.  if TRUE, sample from posterior, if FALSE, sample from prior
#' @return list of two objects:
#' fit: stanfit object
#' data: data list used to fit stanfit object
#' @export
fit_stan_defaults <- function(model, data_list, init = "random", seed = 1,
                              adapt_delta = 0.8, max_treedepth = 10, iter = 2000, evaluate_likelihood = TRUE) {
  data_list$evaluate_likelihood <- as.numeric(evaluate_likelihood)
  chains <- 3 # number of NUTS chains to run in parallel
  set.seed(seed) # set random number generator seed for reproducibility
  # fit model to data using No U-Turn sampler
  if(is.character(model)) {
    model <- rstan::stan_model(get_model_filename(model))
  }
  if(seed > 1) {
    fit <- rstan::sampling(model, data=data_list,  chains=chains,verbose = TRUE, init = init,
                           iter = iter,
                           seed = seed,
                           control = list(adapt_delta = adapt_delta,
                                          max_treedepth = max_treedepth))
  } else {
    fit <- rstan::sampling(model, data=data_list,  chains=chains,verbose = TRUE, init = init,
                           iter = iter,
                           control = list(adapt_delta = adapt_delta,
                                          max_treedepth = max_treedepth))
  }

  list(fit = fit, data = data_list)
}

#' extract samples from fit object and thin
#' 
#' @param fit stanfit object
#' @param N_samples number of samples to extract.  if 0 or greater than the number of available samples, extract all
#' @param par_names (optional) character vector: namse of parameters to extract.  if missing, extract all
#' @param drop if par_names is specified and is of length 1, and drop = TRUE, return vector;
#' otherwise return matrix
#' @param by_chain if TRUE, return list of matrices with one matrix per chain.  Then N_samples is interpreted as per
#' chain.  if FALSE, return one matrix for the whole chain.
#' @return matrix with named columns containing samples
#' @export
extract_fit <- function(fit, N_samples, par_names, drop = FALSE, by_chain = FALSE) {
  ext_fit <- as.matrix(fit)
  if(missing(par_names)) {
    par_names <- colnames(ext_fit)
  }
  if(by_chain) {
    N_chains <- length(fit@inits)
  } else {
    N_chains <- 1
  }
  chain_length <- nrow(ext_fit) / N_chains
  N_samples <- adjust_N_samples(N_samples, chain_length)
  thinned_samples <- round(seq(1, chain_length, length.out = N_samples))
  if(by_chain) {
    ext_fit <- lapply(seq_len(N_chains), 
                      function(x) ext_fit[seq(chain_length * (x - 1) + 1, x * chain_length),])
    ext_fit <- lapply(ext_fit, function(x) x[thinned_samples, par_names, drop = drop])
  } else {
    ext_fit <- ext_fit[thinned_samples, par_names, drop = drop]
  }

  ext_fit
}

#' make dir for Rmd outputs
#' @return string: directory name
#' @importFrom dplyr %>%
#' @export
make_dir_Rmd <- function() {
  dir_name <- knitr::current_input() %>%
    sub(".Rmd", "_rds/", ., fixed = TRUE, ignore.case = TRUE)
  if(!dir.exists(dir_name)) {
    dir.create(dir_name)
  }
  dir_name
    
}

#' get the row with maximum likelihood from a stanfit object or 
#' 
#' @param fit_obj a stanfit object or matrix containing draws from posterior
#' @return named vector with parameters and generated quantities of maximum likelihood
#' @export
get_max_LL_row <- function(fit_obj) {
  if(inherits(fit_obj, "stanfit")) {
    fit_obj <- as.matrix(fit_obj)
  }
  max_LL_row <- fit_obj[which.max(fit_obj[,"lp__"]),]
}

#' summing log likelihood from saliva and each lymphocyte curve
#' 
#' @param data_list list with elements T_U (number of data points for saliva curve)
#' and T_L (vector of number of data points for each lymphocyte curve)
#' @param LL_row vector with elements named log_lik[i], each the log likelihood for a data point
#' @return vector of length 1 + length(T_L): log likelihood for each curve
#' @importFrom dplyr %>%
#' @export
calc_LL_curve <- function(data_list, LL_row) {
  LL_idx <- data_list$T_U + c(0, cumsum(data_list$T_L))
  LL_each <- cumsum(LL_row[grepl("log_lik", names(LL_row))])[LL_idx] %>%
    c(0, .) %>%
    diff %>%
    unname
  LL_each
}

#' wrapper to vectorise lymphocyte parameter names
#' 
#' @param pop number of lymphocyte populations
#' @param scalar_pars character vector.  names of scalar parameters
#' @param vector_pars character vector.  names of vector parameters
#' @return names of scarla and vectorised parameters
#' @importFrom dplyr %>%
#' @export
vectorise_lymphocyte_pars <- function(pop, scalar_pars = c("frac", "delta"),
                                      vector_pars = c("pb_w", "dstar", "delay")) {
  vector_pars_indexed <- lapply(vector_pars, 
                                function(x) vcapply(seq_len(pop), index_par_name, par_name = x)) %>%
    unlist
  all_par_names <- c(scalar_pars, vector_pars_indexed)
  all_par_names
}

#' specify values of fixed parameters
#' 
#' @param cells (optional) if "Granulocytes" or "Monocytes" specified, return R and deltaT for these
#' cell populations
#' @return named list of fixed parameters
#' @export
get_fixed_pars <- function(cells) {
  pars <- list(phase_one_end = 7, label_end = 49)
  if(!missing(cells)) {
    if(cells == "Granulocytes") {
      pars <- c(pars, list(R = 0.26, deltaT = 5.8))
    } else {
      pars <- c(pars, list(R = 0.87, deltaT = 1.6))
    }
  }
  pars
}

#' get vector of parameter names added by each layer of the model
#' @return named list of character vectors
#' @export
get_par_names <- function() {
  list(U =  c("frac", "delta"),
  Lb = c("z", "b_w"),
  L = c("pb_w", "dstar", "delay"))
}

suppress_final_line_warning <- function(w) {
  if(any(grepl("incomplete final line found on", w, fixed = TRUE))) {
    invokeRestart("muffleWarning")
  }
}

read_code_unspace <- function(filename) {
  code <- withCallingHandlers(readLines(filename), warning = suppress_final_line_warning)
  # remove spaces
  code <- gsub(" ", "", code)
}

#' return names of all functions in a .R file
#' 
#' @param filename location of .R file
#' @return character vector where each entry is the name of a function in the file
get_func_names <- function(filename) {
  # read the .R file, suppressing warnings about incomplete final line
  code <- read_code_unspace(filename)
  # identify lines which define functions
  func_str <- "<-function("
  potential_func_lines <- grep(func_str, code, fixed = TRUE)
  
  # determine whether each function is an outermost function
  is_outside_func <- function(code, potential_func_line) {
    # assume functions on first line are outside functions
    if(potential_func_line == 1) {
      return(TRUE)
    }
    # paste all code up to name of potential function
    pasted_code <- paste0(code[seq_len(potential_func_line - 1)], collapse = "")
    # potential_func_subline <- strsplit(code, split = func_str, fixed = TRUE)
    # potential_func_subline <- potential_func_subline[1]
    # pasted_code <- paste0(pasted_code, potential_func_subline, collapse = "")
    # count number of open and close curly brackets up to potential function name
    count_characters <- function(pasted_code, char_in) {
      n <- gregexpr(char_in, pasted_code, fixed = TRUE)
      length(n[[1]])
    }
    n_brackets <- vapply(c("{", "}"), function(x) count_characters(pasted_code, x), numeric(1))
    # the functino is an outermost function if the number of open and close brackets is the same
    n_brackets[1] == n_brackets[2]
  }
  
  func_lines <- potential_func_lines[vapply(potential_func_lines, 
                                            function(x) is_outside_func(code, x),
                                            logical(1))]
  # split off the function names
  func_lines <- strsplit(code[func_lines], split = func_str, fixed = TRUE)
  func_lines <- vapply(func_lines, function(x) x[1], character(1))
  func_lines
}

#' find the file containing a function
#' 
#' @param func_name name of function
#' @param dir_name directory in which to search (only searches in that directory
#' + subdirectories one level down)
#' @return the name of the file which has the function
#' @export
find_func_file <- function(func_name, dir_name = "~/git_repos/kirlabelling/R") {
  # list current directory + directories one level down
  dirs <- list.dirs(path = dir_name, recursive = FALSE)
  dirs <- c(dirs, dir_name)
  # exclude Rproj.user -- not a relevant directory
  dirs <- dirs[!grepl("Rproj.user", dirs, fixed = TRUE)]
  # list files in those directories
  filenames <- lapply(dirs, function(x) list.files(x, pattern="\\.R$", full.names=TRUE))
  # list functions in those files
  func_names <- lapply(filenames, function(x) lapply(x, get_func_names))
  # find index of file containing target function
  in_file <- lapply(func_names, function(x) vlapply(x, function(y) any(func_name == y)))
  # find directory of that file
  in_dir <- vlapply(in_file, any)
  # return filename
  
  if(length(filenames[in_dir]) == 0) {
    stop("function not found")
  }
  func_file <- filenames[in_dir][[1]][in_file[in_dir][[1]]]
  # open file in Rstudio
  file.edit(func_file)
  func_file
}

#' solve an ODE with error handling
#' 
#' if an error occurs, solve with a smaller time interval
#' for solution elements except the time vector, set to tol if the solution is below tol
#' 
#' @param mod an ode_system object generated by odin
#' @param solving_time times at which to solve the ODE
#' @param tol numeric vector of length 1: tolerance
#' @return a deSolve object: solution of the ODE
solve_ODE_error_handling <- function(mod, solving_time, tol = 1e-6) {
  
  ## maximum number of times to try to solve with smaller time interval
  max_fail_iter <- 5
  
  sol <- tryCatch(mod$run(solving_time), error = function(e) "error") # solve ODEs
  
  ## if error occurred while solving ODEs, try smaller solving interval
  if(!is.matrix(sol)) {
    fail_iter <- 0
    n_div <- 10
    solving_time_temp <- solving_time
    while(fail_iter < max_fail_iter && !is.matrix(sol)) {
      # subdivide each time step into n_div steps
      solving_time_temp <- interpolate_vector(solving_time_temp, n_div)
      # resolve
      sol <- tryCatch(mod$run(solving_time_temp), error = function(e) "error")
      fail_iter <- fail_iter + 1
    }
    # retrive solution for original time steps
    sol <- sol[solving_time_temp %in% solving_time,]
  }
  
  # if ODE solver exits early, pad with last value
  sol <- sol[sol[,1] %in% solving_time,]
  if(nrow(sol) < length(solving_time)) {
    last_row <- sol[nrow(sol),]
    rep_last_rows <- rep(list(last_row), length(solving_time) - nrow(sol))
    sol <- do.call(rbind, c(list(sol), rep_last_rows))
    sol[,1] <- solving_time
  }
  
  # for solution elements except the time vector, set to tol if the solution is below tol
  non_time_elements <- sol[,-1]
  non_time_elements[non_time_elements < tol] <- tol
  sol[,-1] <- non_time_elements
  
  sol
}

#' interpolate ordered numeric vector into n_div intervals betwen successive values
#' 
#' @param vec numeric vector to interpolate
#' @param n_div numeric vector of length 1: number of intervals into which to subdivide
#' successive values of the vector
#' @return numeric vector of length (length(vec) - 1) * n_div + 1: interpolated vector
interpolate_vector <- function(vec, n_div) {
  x <- seq(1, length(vec), by = 1/n_div)
  temp <- approx(seq_along(vec), vec, xout = x)
  vec_out <- temp$y
  vec_out
}

get_mod_info <- function(mod_filename) {
  # read the .R file, suppressing warnings about incomplete final line
  code <- read_code_unspace(mod_filename)
  
  get_pars <- function() {
    # identify lines which define parameters
    par_str <- "<-user("
    par_lines <- grep(par_str, code, fixed = TRUE)
    code <- code[par_lines] %>%
      lapply(., trimws)
    # exclude commented lines
    code <- code[substr(code, 1, 1) != "#"]
    
    extract_par_name_from_line <- function(par_str) {
      par_str <- strsplit(par_str, "<-user()", fixed = TRUE)[[1]][1]
      strsplit(par_str, "[", fixed = TRUE)[[1]][1]
    }
    
    vcapply(code, extract_par_name_from_line) %>% unique
  }
  
  get_compartments <- function() {
    compartment_str <- "initial("
    compartment_lines <- grep(compartment_str, code, fixed = TRUE)
    code <- code[compartment_lines] %>%
      lapply(., trimws)
    # exclude commented lines
    code <- code[substr(code, 1, 1) != "#"]
    
    extract_compartment_name_from_line <- function(compartment_str) {
      compartment_str <- strsplit(compartment_str, "initial(", fixed = TRUE)[[1]][2]
      compartment_str <- strsplit(compartment_str, "[", fixed = TRUE)[[1]][1]
      strsplit(compartment_str, ")", fixed = TRUE)[[1]][1]
    }
    
    vcapply(code, extract_compartment_name_from_line) %>% unique
  }
  
  list(pars = get_pars(), compartments = get_compartments())
  
}

#' get name of directory which package is in
#'
#' @return a string
get_package_dir <- function() {
  "~/git_repos/kirlabelling/"
}