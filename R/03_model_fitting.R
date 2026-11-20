#' @title Model Fitting Functions for mfeargen
#' @description Core functions for fitting Bayesian mixture models
#' @name model_fitting
#' @importFrom stats runif rnorm rgamma
NULL

#' Fit mfeargen Bayesian Model
#'
#' @description
#' Fits a Bayesian mixture model to fear generalization data using JAGS.
#' Supports both 1CS (simple conditioning) and 2CS (differential conditioning) designs.
#'
#' @param data Data frame with required columns: participant, trials, US_expect
#' @param model_type Character. "1cs" or "2cs" 
#' @param n.iter Integer. Number of MCMC iterations (default: 10000)
#' @param n.burnin Integer. Number of burn-in iterations (default: 5000)
#' @param n.chains Integer. Number of MCMC chains (default: 3)
#' @param n.thin Integer. Thinning interval (default: 2)
#' @param lambda_criterion Numeric. Criterion for overgeneralization boundary (default: 0.7)
#' @param priors mfeargen_priors object. Custom priors from set_priors() (optional)
#' @param prepared_data List. Pre-prepared JAGS data from prepare_mfeargen_data() (optional)
#' @param parallel Logical. Run chains in parallel (default: TRUE)
#' @param verbose Logical. Print progress messages (default: TRUE)
#' @param rating_scale Numeric vector c(min, max). Known response scale bounds.
#'   If omitted, uses simulation metadata or the observed range.
#' @param seed Integer. Optional seed for initialization and sampling.
#' @param observation_model Either `"normal"` for an uncensored Gaussian response
#'   or `"censored"` for ratings clipped at their instrument bounds. If NULL,
#'   uses the data's observation-model metadata, otherwise `"normal"`.
#' @param save_predictive Retain a full array of replicated ratings in the fit.
#'   Defaults to FALSE to reduce memory use; posterior_predict generates these
#'   ratings when needed.
#' @param ... Additional arguments passed to jagsUI
#' 
#' @return mfeargen object containing:
#'   - jags_fit: MCMC output from JAGS
#'   - model_type: Model specification used
#'   - data: Original data
#'   - convergence: Convergence diagnostics
#'   - group_membership: Predicted group assignments
#'   
#' @export
#' @import jagsUI
#' @examples
#' \dontrun{
#' # Simulate data
#' data <- simulate_mfeargen(n_participants = 20)
#' 
#' # Fit model
#' fit <- fit_mfeargen(data, model_type = "1cs", n.iter = 5000)
#' 
#' # View results
#' summary(fit)
#' plot(fit)
#' }
fit_mfeargen <- function(data,
                        model_type = c("1cs", "2cs"),
                        n.iter = 10000,
                        n.burnin = 5000,
                        n.chains = 3,
                        n.thin = 2,
                        lambda_criterion = 0.7,
                        priors = NULL,
                        prepared_data = NULL,
                        parallel = TRUE,
                        verbose = TRUE,
                        rating_scale = NULL,
                        seed = NULL,
                        observation_model = NULL,
                        save_predictive = FALSE,
                        ...) {
  model_type <- match.arg(model_type)
  if (!is.data.frame(data)) stop("data must be a data frame")
  if (is.null(observation_model)) {
    observation_model <- attr(prepared_data, "observation_model") %||%
      attr(data, "observation_model") %||% "normal"
  }
  observation_model <- match.arg(observation_model, c("normal", "censored"))
  counts <- list(n.iter = n.iter, n.burnin = n.burnin,
                 n.chains = n.chains, n.thin = n.thin)
  for (nm in names(counts)) {
    value <- counts[[nm]]
    minimum <- if (nm == "n.burnin") 0 else 1
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value != floor(value) || value < minimum) {
      stop(nm, " must be a whole number >= ", minimum)
    }
  }
  if (n.iter <= n.burnin) stop("n.iter must be greater than n.burnin")
  if (floor((n.iter - n.burnin) / n.thin) < 2) {
    stop("Sampling settings must retain at least two draws per chain")
  }
  for (nm in c("parallel", "verbose", "save_predictive")) {
    value <- get(nm)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(nm, " must be TRUE or FALSE")
    }
  }
  if (!is.null(priors)) validate_priors(priors)
  dots <- list(...)
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))) ||
                       anyDuplicated(names(dots)))) {
    stop("Additional JAGS arguments must have distinct names")
  }
  unknown <- setdiff(names(dots), names(formals(jagsUI::jags)))
  if (length(unknown)) stop("Unknown JAGS arguments: ", paste(unknown, collapse = ", "))
  reserved <- intersect(names(dots), c("data", "inits", "parameters.to.save", "model.file"))
  if (length(reserved)) stop("Reserved JAGS arguments: ", paste(reserved, collapse = ", "))
  if (isTRUE(dots$codaOnly)) stop("codaOnly = TRUE is incompatible with mfeargen post-processing")
  dic <- dots$DIC %||% TRUE
  if (!is.logical(dic) || length(dic) != 1L || is.na(dic)) stop("DIC must be TRUE or FALSE")
  dots$DIC <- NULL
  if (!is.null(dots$n.adapt)) validate_scalar(dots$n.adapt, "n.adapt", lower = 0, integer = TRUE)
  if (!is.null(dots$n.cores)) validate_scalar(dots$n.cores, "n.cores", lower = 1, integer = TRUE)
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) {
      stop("seed must be a non-negative integer")
    }
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  if (is.null(prepared_data)) {
    if (verbose) message("Preparing data for JAGS...")
    jags_data <- prepare_mfeargen_data(
      data, model_type, lambda_criterion, rating_scale = rating_scale,
      observation_model = observation_model
    )
  } else {
    jags_data <- prepared_data
    validate_prepared_data(jags_data, model_type, data, observation_model)
    if (!is.null(rating_scale) &&
        !isTRUE(all.equal(as.numeric(rating_scale), as.numeric(attr(jags_data, "rating_scale") %||% jags_data$rating_scale)))) {
      stop("rating_scale does not match prepared_data$rating_scale")
    }
  }
  validate_prepared_data(jags_data, model_type)
  attr(jags_data, "observation_model") <- observation_model
  require_jags_backend()

  model_file <- system.file("jags", paste0("model_", model_type, ".txt"),
                            package = "mfeargen")
  if (!nzchar(model_file) || !file.exists(model_file)) {
    stop("JAGS model file not found; reinstall mfeargen")
  }
  model_text <- readLines(model_file, warn = FALSE)
  if (!is.null(priors)) model_text <- apply_model_priors(model_text, priors)
  observation <- compile_observation_model(model_text, jags_data, observation_model)
  model_text <- observation$model_code
  if (!is.null(priors) || observation_model == "censored") {
    model_file <- tempfile("mfeargen-", fileext = ".txt")
    writeLines(model_text, model_file)
    on.exit(unlink(model_file), add = TRUE)
  }
  inits <- initialize_mcmc_parameters(jags_data, n.chains)
  parameters <- c("alpha", "lambda", "gp", "pi", "w0", "w1", "sigma",
                  "theta", "loglik", "alpha_mu", "alpha_kappa",
                  "lambda_mu", "lambda_sigma", "w0_mu", "w0_sigma", "w1_a", "w1_b")
  if (save_predictive) parameters <- c(parameters, "y_pre")
  # Prior metadata belongs to the R object, not JAGS' data environment.
  sampling_data <- observation$data[!grepl("^prior_", names(observation$data))]
  if (verbose) message("Running MCMC sampling...")
  fit <- tryCatch(do.call(jagsUI::jags, c(list(
    data = sampling_data, inits = inits, parameters.to.save = parameters,
    model.file = model_file, n.chains = n.chains, n.iter = n.iter,
    n.burnin = n.burnin, n.thin = n.thin, parallel = parallel,
    DIC = dic, verbose = verbose
  ), dots)), error = function(error) {
    stop("JAGS sampling failed: ", conditionMessage(error), call. = FALSE)
  })
  fit <- normalize_pointwise_output(fit, dim(jags_data$y))
  participants <- attr(jags_data, "participants")
  if (is.null(participants)) participants <- unique(data$participant)
  group_membership <- extract_group_membership(fit, participants)
  convergence <- check_convergence(fit, verbose = verbose)
  result <- list(
    jags_fit = fit, model_type = model_type, data = data, jags_data = jags_data,
    group_membership = group_membership, lambda_threshold = jags_data$lambda_threshold,
    convergence = convergence, rating_scale = attr(jags_data, "rating_scale") %||% jags_data$rating_scale,
    n.iter = n.iter, n.burnin = n.burnin,
    n.chains = n.chains, n.thin = n.thin, seed = seed, priors = priors,
    observation_model = observation_model, save_predictive = save_predictive,
    model_code = model_text, call = match.call()
  )
  class(result) <- "mfeargen"
  if (verbose) message("Model fitting complete (", jags_data$Nparticipants, " participants).")
  result
}

validate_prepared_data <- function(jags_data, model_type, data = NULL, observation_model = NULL) {
  if (!is.list(jags_data)) {
    stop("prepared_data must be a list created by prepare_mfeargen_data()")
  }
  required <- c("y", "Nparticipants", "Ntrials", "rating_scale", "lambda_threshold")
  required <- c(required, if (model_type == "1cs") {
    c("k", "r", "d_phy", "d_per_old")
  } else {
    c("k_plus", "k_minus", "r_plus", "r_minus", "d_p_phy", "d_m_phy", "d_p_per", "d_m_per")
  })
  missing <- setdiff(required, names(jags_data))
  if (length(missing)) stop("prepared_data is missing: ", paste(missing, collapse = ", "))
  for (nm in c("Nparticipants", "Ntrials")) {
    value <- jags_data[[nm]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != floor(value)) stop("Invalid prepared_data$", nm)
  }
  validate_rating_scale(jags_data$rating_scale)
  if (!identical(as.numeric(jags_data$rating_scale), c(0, 10))) {
    stop("prepared_data must use the normalized 0-10 scale; recreate it with prepare_mfeargen_data()")
  }
  declared_type <- attr(jags_data, "model_type")
  if (!is.null(declared_type) && declared_type != model_type) {
    stop("model_type does not match prepared_data")
  }
  expected <- c(jags_data$Nparticipants, jags_data$Ntrials)
  for (nm in setdiff(required, c("Nparticipants", "Ntrials", "rating_scale", "lambda_threshold"))) {
    value <- jags_data[[nm]]
    if (!is.numeric(value) || !identical(as.integer(dim(value)), as.integer(expected)) ||
        any(is.infinite(value)) || (nm != "y" && anyNA(value))) {
      stop("Invalid prepared_data$", nm, "; expected a numeric participant-by-trial matrix")
    }
  }
  scale <- jags_data$rating_scale
  if (!is.numeric(scale) || length(scale) != 2L || any(!is.finite(scale)) || scale[1] >= scale[2]) {
    stop("prepared_data$rating_scale must have two finite increasing bounds")
  }
  if (!is.numeric(jags_data$lambda_threshold) || length(jags_data$lambda_threshold) != 1L ||
      !is.finite(jags_data$lambda_threshold) || jags_data$lambda_threshold <= 0) {
    stop("prepared_data$lambda_threshold must be positive")
  }
  indicators <- if (model_type == "1cs") "k" else c("k_plus", "k_minus")
  for (name in indicators) {
    if (any(!jags_data[[name]] %in% c(0, 1))) stop("prepared_data$", name, " must be coded 0 or 1")
  }
  if (model_type == "2cs" && any(jags_data$k_plus + jags_data$k_minus > 1)) {
    stop("prepared_data learning indicators cannot select both cues on one trial")
  }
  outcomes <- if (model_type == "1cs") "r" else c("r_plus", "r_minus")
  for (name in outcomes) {
    allowed <- if (name == "r_minus") c(-1, 0, 1) else c(0, 1)
    if (any(!jags_data[[name]] %in% allowed)) stop("Invalid reinforcement targets in prepared_data$", name)
  }
  distances <- grep("^d_", required, value = TRUE)
  for (name in distances) {
    if (any(jags_data[[name]] < 0)) stop("prepared_data$", name, " must be non-negative")
  }
  if (!is.null(data)) {
    bounds <- attr(jags_data, "rating_scale")
    if (is.null(bounds)) stop("prepared_data is missing original rating_scale metadata; recreate it with prepare_mfeargen_data()")
    expected_data <- prepare_mfeargen_data(data, model_type, rating_scale = bounds,
                                          observation_model = observation_model)
    same_values <- vapply(setdiff(required, "lambda_threshold"), function(name) {
      isTRUE(all.equal(jags_data[[name]], expected_data[[name]], check.attributes = FALSE))
    }, logical(1))
    same_mapping <- vapply(c("participants", "row_index", "trial_counts"), function(name) {
      identical(attr(jags_data, name), attr(expected_data, name))
    }, logical(1))
    if (!all(same_values) || !all(same_mapping)) {
      stop("prepared_data does not match data; recreate it after changing ratings, trial order, IDs, or covariates")
    }
  }
  invisible(TRUE)
}

normalize_pointwise_output <- function(fit, panel_dimensions) {
  n_draws <- sum(vapply(fit$samples, nrow, integer(1)))
  summaries <- c("mean", "sd", "q2.5", "q25", "q50", "q75", "q97.5", "Rhat", "n.eff")
  for (node in intersect(c("theta", "loglik", "y_pre"), names(fit$sims.list))) {
    draws <- fit$sims.list[[node]]
    if (length(draws) != n_draws * prod(panel_dimensions)) {
      stop("Unexpected posterior dimensions for ", node)
    }
    fit$sims.list[[node]] <- array(draws, c(n_draws, panel_dimensions))
    for (summary in summaries) {
      values <- fit[[summary]][[node]]
      if (!is.null(values) && length(values) == prod(panel_dimensions)) {
        fit[[summary]][[node]] <- array(values, panel_dimensions)
      }
    }
  }
  fit
}

require_jags_backend <- function() {
  available <- tryCatch(requireNamespace("jagsUI", quietly = TRUE) &&
                          requireNamespace("rjags", quietly = TRUE),
                        error = function(error) FALSE)
  if (!available) {
    stop("Model fitting requires JAGS and the jagsUI/rjags R packages. Install JAGS 4.3 or later, restart R, then install.packages(c('rjags', 'jagsUI')).", call. = FALSE)
  }
  invisible(TRUE)
}

compile_observation_model <- function(model_text, jags_data, observation_model = "normal") {
  observation_model <- match.arg(observation_model, c("normal", "censored"))
  if (observation_model == "normal") return(list(model_code = model_text, data = jags_data))
  y <- jags_data$y
  dimensions <- length(dim(y))
  if (!dimensions %in% c(2L, 3L)) stop("Censored responses require a two- or three-dimensional response array")
  bounds <- jags_data$rating_scale %||% c(jags_data$A, jags_data$K)
  validate_rating_scale(bounds)
  if (any(y < bounds[1] | y > bounds[2], na.rm = TRUE)) stop("Censored ratings must be within their rating bounds")
  coordinates <- if (dimensions == 2L) c("i", "j") else c("i", "t", "j")
  index <- paste(coordinates, collapse = ",")
  node <- function(name) paste0(name, "[", index, "]")
  replace_relation <- function(name, relation, replacement) {
    rows <- which(startsWith(trimws(model_text), paste0(node(name), " ", relation)))
    if (length(rows) != 1L) stop("Cannot locate the ", name, " relation in the JAGS model")
    model_text[rows] <<- paste0("      ", replacement)
  }
  replace_relation("y", "~", paste0(node("y_cont"), " ~ dnorm(", node("theta"), ", tau_obs[i])"))
  replace_relation("loglik", "<-", paste0(node("loglik"), " <- ", node("observed"),
    " * ifelse(", node("censor_code"), " == 0, logdensity.norm(", node("response"), ", ",
    node("theta"), ", tau_obs[i]), log_tail[", node("censor_index"), "])") )
  replace_relation("y_pre", "~", paste0(node("y_pre_raw"), " ~ dnorm(", node("theta"), ", tau_obs[i])\n      ",
    node("y_pre"), " <- max(", bounds[1], ", min(", bounds[2], ", ", node("y_pre_raw"), "))"))

  observed <- !is.na(y)
  lower <- observed & y == bounds[1]
  upper <- observed & y == bounds[2]
  boundary <- lower | upper
  indices <- which(boundary, arr.ind = TRUE)
  count <- nrow(indices)
  if (!count) indices <- matrix(1L, nrow = 1L, ncol = dimensions)
  sampling_data <- jags_data
  sampling_data$y <- NULL
  sampling_data$y_cont <- y
  sampling_data$y_cont[boundary] <- NA_real_
  sampling_data$response <- y
  sampling_data$response[!observed] <- 0
  sampling_data$observed <- 1L * observed
  sampling_data$censor_code <- -1L * lower + 1L * upper
  sampling_data$censor_index <- array(1L, dim(y))
  if (count) sampling_data$censor_index[boundary] <- seq_len(count)
  sampling_data$N_censored <- max(1L, count)
  sampling_data$censor_active <- rep(as.integer(count > 0L), max(1L, count))
  sampling_data$censor_side <- if (count) sampling_data$censor_code[boundary] else -1L
  sampling_data$zeros <- rep(0L, max(1L, count))
  for (column in seq_along(coordinates)) {
    sampling_data[[paste0("censor_", coordinates[column])]] <- indices[, column]
  }
  censor_index <- paste0("censor_", coordinates, "[h]", collapse = ",")
  tail_mean <- paste0("theta[", censor_index, "]")
  # The normal-tail expansion avoids underflow without imposing a probability floor.
  expansion <- "1 - tail_inv[h] * (1 - tail_inv[h] * (3 - tail_inv[h] * (15 - tail_inv[h] * (105 - tail_inv[h] * (945 - tail_inv[h] * (10395 - tail_inv[h] * (135135 - tail_inv[h] * (2027025 - tail_inv[h] * 34459425))))))))"
  tail_code <- c(
    "  for (h in 1:N_censored) {",
    paste0("    tail_z[h] <- ifelse(censor_side[h] == -1, ", bounds[1], " - ", tail_mean,
           ", ", tail_mean, " - ", bounds[2], ") * sqrt(tau_obs[censor_i[h]])"),
    "    tail_x[h] <- max(-tail_z[h], 10)",
    "    tail_inv[h] <- pow(tail_x[h], -2)",
    paste0("    tail_series[h] <- ", expansion),
    "    log_tail[h] <- ifelse(tail_z[h] < -10, logdensity.norm(tail_x[h], 0, 1) - log(tail_x[h]) + log(tail_series[h]), log(phi(max(tail_z[h], -10))))",
    "    zeros[h] ~ dpois(-censor_active[h] * log_tail[h])",
    "  }"
  )
  model_text <- append(model_text, tail_code, after = length(model_text) - 1L)
  list(model_code = model_text, data = sampling_data)
}

apply_model_priors <- function(model_text, priors) {
  validate_priors(priors)
  nodes <- c(lambda_mu = "lambda_mu", lambda_sigma = "lambda_sigma",
             alpha_mu = "alpha_mu", alpha_kappa = "alpha_kappa",
             sigma1 = "sigma[1]", sigma2 = "sigma[2]", pi = "pi[1:4]",
             w0_mu = "w0_mu", w0_sigma = "w0_sigma", w1_a = "w1_a", w1_b = "w1_b")
  for (nm in names(nodes)) {
    node <- nodes[[nm]]
    line <- which(startsWith(trimws(model_text), paste0(node, " ~")))
    if (length(line) != 1L) stop("Cannot locate prior for ", nm)
    prior <- priors[[nm]]
    if (nm == "lambda_mu" && prior$dist == "norm" && is.null(prior$truncate)) {
      prior$truncate <- c(0, Inf)
    }
    if (nm == "alpha_mu") prior$truncate <- c(1e-9, 1 - 1e-9)
    model_text[line] <- paste0("  ", prior_to_jags(prior, node))
  }
  model_text
}

#' Initialize MCMC Parameters
#'
#' @param jags_data List. Prepared JAGS data.
#' @param n.chains Integer. Number of chains.
#' @return List of initial values for each chain.
#' @keywords internal
initialize_mcmc_parameters <- function(jags_data, n.chains) {
  n <- jags_data$Nparticipants
  lapply(seq_len(n.chains), function(chain) {
    list(
      alpha_1 = runif(n, 0.2, 0.8),
      lambda_1 = runif(n, 1.1, 3) * jags_data$lambda_threshold,
      lambda_2 = runif(n, 0.1, 0.9) * jags_data$lambda_threshold,
      w0 = rnorm(n, 0, 0.5), w1_1 = runif(n, 1.1, 3),
      gp = sample.int(4L, n, replace = TRUE), pi = rep(0.25, 4),
      .RNG.name = "base::Mersenne-Twister",
      .RNG.seed = sample.int(.Machine$integer.max, 1)
    )
  })
}


#' Extract Group Membership
#'
#' @description
#' Extracts the most probable group membership for each participant.
#'
#' @param fit JAGS output object
#' @param participants Original participant identifiers.
#' @return Data frame with group assignments and probabilities
#' @keywords internal
extract_group_membership <- function(fit, participants = NULL) {
  samples <- fit$sims.list$gp
  if (is.null(samples)) {
    warning("Group membership samples (gp) not found in model output")
    return(NULL)
  }
  if (is.null(dim(samples))) samples <- matrix(samples, ncol = 1L)
  n <- ncol(samples)
  if (is.null(participants)) participants <- seq_len(n)
  if (length(participants) != n) stop("Participant IDs do not match group samples")
  probs <- vapply(seq_len(4L), function(group) colMeans(samples == group), numeric(n))
  probs <- matrix(probs, nrow = n, ncol = 4L)
  groups <- max.col(probs, ties.method = "first")
  labels <- c("Non-Learner", "Overgeneralizer", "Physical", "Perceptual")
  membership <- data.frame(participant = participants, group = groups,
                           group_label = labels[groups], assigned_group = labels[groups],
                           max_probability = apply(probs, 1L, max))
  colnames(probs) <- paste0("prob_group_", seq_len(4L))
  cbind(membership, probs)
}


#' Check MCMC convergence and effective sample sizes
#'
#' @param fit A fitted mfeargen model or JAGS output containing separate chains.
#' @param threshold Upper limit for rank-normalized split Rhat.
#' @param verbose Print diagnostic messages.
#' @param min_ess_per_chain Minimum bulk and tail effective sample sizes per chain.
#' @return A list with an overall status, a parameter-level diagnostic table,
#'   Rhat and effective sample size summaries. `converged = NA` means that the
#'   available draws cannot establish convergence, including single-chain runs.
#' @details Uses rank-normalized split and folded Rhat, bulk ESS, and tail ESS.
#'   Diagnostics exclude trial expectations, pointwise likelihoods and predictive
#'   draws. Structurally zero group-specific coefficients are listed separately;
#'   other constant parameters remain unassessed. Passing these diagnostics is
#'   necessary but does not establish model adequacy or rule out unexplored modes.
#' @references Vehtari et al. (2021). Rank-normalization, folding, and
#'   localization: An improved Rhat for assessing convergence of MCMC.
#'   Bayesian Analysis, 16(2), 667-718. \doi{10.1214/20-BA1221}.
#' @export
check_convergence <- function(fit, threshold = 1.01, verbose = TRUE,
                              min_ess_per_chain = 100) {
  validate_scalar(threshold, "threshold", lower = 1, lower_open = TRUE)
  validate_scalar(min_ess_per_chain, "min_ess_per_chain", lower = 0,
                  lower_open = TRUE)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE")
  }
  if (inherits(fit, "mfeargen") || inherits(fit, "mfeargen_hierarchical")) fit <- fit$jags_fit
  if (!is.list(fit) && !is.numeric(fit)) stop("fit must contain JAGS output or numeric Rhat values")
  chains <- if (is.list(fit)) fit$samples else NULL
  n_chains <- length(chains)
  method <- "rank-normalized split Rhat, bulk ESS and tail ESS"
  reason <- NULL
  diagnostic_warnings <- character(0)
  if (n_chains) {
    columns <- colnames(chains[[1]])
    excluded <- grepl("^(theta|loglik|y_pre|y_pre_raw|y_cont|deviance|tau_obs|zeros|tail_[A-Za-z_]+|log_tail)(\\[|$)", columns)
    columns <- columns[!excluded]
    lengths <- vapply(chains, nrow, integer(1))
    if (!length(columns) || any(lengths != lengths[1])) {
      stop("Chain output must contain equally sized draws for model parameters")
    }
    draws <- array(NA_real_, c(lengths[1], n_chains, length(columns)),
                    dimnames = list(NULL, NULL, columns))
    for (chain in seq_len(n_chains)) {
      if (!all(columns %in% colnames(chains[[chain]]))) {
        stop("Monitored parameters differ across chains")
      }
      draws[, chain, ] <- as.matrix(chains[[chain]])[, columns, drop = FALSE]
    }
    structural <- vapply(seq_along(columns), function(index) {
      values <- draws[, , index]
      grepl("^(alpha|lambda|w1)(_stable)?(\\[|$)", columns[index]) &&
        all(is.finite(values)) && all(values == 0)
    }, logical(1))
    diagnostics <- data.frame(parameter = columns, rhat = NA_real_,
                               ess_bulk = NA_real_, ess_tail = NA_real_,
                               status = "unassessed", stringsAsFactors = FALSE)
    diagnostics$status[structural] <- "structural_constant"
    if (n_chains >= 2L && lengths[1] >= 4L && any(!structural)) {
      summaries <- withCallingHandlers(posterior::summarise_draws(
        posterior::as_draws_array(draws[, , !structural, drop = FALSE]),
        "rhat", "ess_bulk", "ess_tail"), warning = function(warning) {
          diagnostic_warnings <<- unique(c(diagnostic_warnings, conditionMessage(warning)))
          invokeRestart("muffleWarning")
        })
      diagnostics[!structural, c("rhat", "ess_bulk", "ess_tail")] <-
        as.data.frame(summaries[c("rhat", "ess_bulk", "ess_tail")])
    } else {
      reason <- if (n_chains < 2L) "At least two independent chains are required" else "Too few retained draws to assess convergence"
    }
  } else {
    rhat <- if (is.numeric(fit)) fit else if (!is.null(fit$Rhat)) {
      unlist(fit$Rhat)
    } else if (!is.null(fit$summary) && "Rhat" %in% colnames(fit$summary)) {
      fit$summary[, "Rhat"]
    } else numeric(0)
    parameters <- names(rhat)
    if (is.null(parameters)) parameters <- as.character(seq_along(rhat))
    diagnostics <- data.frame(parameter = parameters, rhat = as.numeric(rhat),
                               ess_bulk = NA_real_, ess_tail = NA_real_,
                               status = rep("unassessed", length(rhat)), stringsAsFactors = FALSE)
    structural <- rep(FALSE, length(rhat))
    method <- "supplied Rhat only"
    reason <- "Separate chains are required to calculate rank Rhat and effective sample sizes"
  }
  ess_threshold <- min_ess_per_chain * max(n_chains, 1L)
  assessed <- !structural & is.finite(diagnostics$rhat) &
    is.finite(diagnostics$ess_bulk) & is.finite(diagnostics$ess_tail)
  failed <- !structural & (
    (!is.na(diagnostics$rhat) & diagnostics$rhat >= threshold) |
    (!is.na(diagnostics$ess_bulk) & diagnostics$ess_bulk < ess_threshold) |
    (!is.na(diagnostics$ess_tail) & diagnostics$ess_tail < ess_threshold))
  diagnostics$status[assessed] <- "passed"
  diagnostics$status[failed] <- "failed"
  relevant <- !structural
  unknown <- relevant & diagnostics$status == "unassessed"
  converged <- if (any(failed)) FALSE else if (!any(relevant) || any(unknown)) NA else TRUE
  extreme <- function(values, fun) if (all(is.na(values))) NA_real_ else fun(values, na.rm = TRUE)
  result <- list(
    converged = converged,
    status = if (is.na(converged)) "unassessed" else if (converged) "passed" else "failed",
    max_rhat = extreme(diagnostics$rhat[relevant], max),
    min_neff = extreme(c(diagnostics$ess_bulk[relevant], diagnostics$ess_tail[relevant]), min),
    n_parameters = sum(relevant), n_poor_convergence = sum(failed),
    n_unassessed = sum(unknown), n_chains = n_chains,
    threshold = threshold, ess_threshold = ess_threshold, method = method,
    problematic_parameters = diagnostics$parameter[failed | unknown],
    structural_constants = diagnostics$parameter[structural],
    diagnostics = diagnostics, reason = reason, warnings = diagnostic_warnings
  )
  if (verbose && !isTRUE(converged)) {
    if (identical(converged, FALSE)) {
      message(sum(failed), " model parameters exceed the Rhat or effective sample size limits. Review fit$convergence$diagnostics.")
    } else {
      message(reason %||% "Some model parameters are constant or have unavailable diagnostics; convergence remains unassessed.")
    }
  }
  result
}


#' Extract model parameters
#'
#' @param fit A fitted mfeargen model or JAGS output.
#' @param parameter One parameter name, or "all".
#' @param summary_only Return posterior summaries when TRUE, or retained draws
#'   when FALSE.
#' @return For one parameter, a data frame of posterior summaries with named
#'   parameter rows, or an array of retained draws. For "all", a named list of
#'   parameter blocks, with group membership when available. Existing mean,
#'   standard-deviation, and quantile components are retained for compatibility.
#' @export
#' @examples
#' \dontrun{
#' fit <- fit_mfeargen(simulate_mfeargen(n_participants = 8))
#' extract_parameters(fit, "alpha")
#' extract_parameters(fit, "alpha", summary_only = FALSE)
#' }
extract_parameters <- function(fit, parameter = "all", summary_only = TRUE) {
  if (!is.character(parameter) || length(parameter) != 1L || is.na(parameter) || !nzchar(parameter)) {
    stop("parameter must be one parameter name or 'all'")
  }
  if (!is.logical(summary_only) || length(summary_only) != 1L || is.na(summary_only)) {
    stop("summary_only must be TRUE or FALSE")
  }
  fitted_model <- inherits(fit, "mfeargen") || inherits(fit, "mfeargen_hierarchical")
  output <- if (fitted_model) fit$jags_fit else fit
  if (!is.list(output)) stop("fit must be a fitted model or JAGS output list")
  summary <- output$summary
  if (!is.null(summary) && (!is.matrix(summary) && !is.data.frame(summary))) {
    stop("The JAGS summary must be a matrix or data frame with parameter row names")
  }
  if (summary_only && !is.null(summary)) {
    if (is.null(rownames(summary)) || !"mean" %in% colnames(summary)) {
      stop("The JAGS summary requires parameter row names and a mean column")
    }
    available <- unique(sub("\\[.*$", "", rownames(summary)))
  } else if (summary_only && !is.null(output$mean)) {
    available <- names(output$mean)
  } else if (!summary_only && !is.null(output$sims.list)) {
    available <- names(output$sims.list)
  } else {
    stop(if (summary_only) "Posterior summaries are unavailable; use summary_only = FALSE when retained draws are available" else "Retained posterior draws are unavailable")
  }
  parameters <- if (parameter == "all") available else parameter
  if (!all(parameters %in% available)) stop("Parameter '", parameter, "' not found in model output")
  results <- stats::setNames(vector("list", length(parameters)), parameters)
  for (name in parameters) {
    if (!summary_only) {
      results[[name]] <- output$sims.list[[name]]
      next
    }
    if (!is.null(summary)) {
      selected <- rownames(summary) == name | startsWith(rownames(summary), paste0(name, "["))
      block <- summary[selected, , drop = FALSE]
      column <- function(key) {
        if (key %in% colnames(block)) as.numeric(block[, key]) else rep(NA_real_, nrow(block))
      }
      results[[name]] <- data.frame(
        parameter = rownames(block), mean = column("mean"), sd = column("sd"),
        q2.5 = column("2.5%"), q50 = column("50%"), q97.5 = column("97.5%"),
        row.names = rownames(block))
    } else {
      values <- output$mean[[name]]
      if (!is.numeric(values) || !length(values)) stop("Posterior means are unavailable for ", name)
      labels <- if (!is.null(dim(values))) {
        indices <- arrayInd(seq_along(values), dim(values))
        paste0(name, "[", apply(indices, 1, paste, collapse = ","), "]")
      } else if (length(values) == 1L) name else paste0(name, "[", seq_along(values), "]")
      component <- function(key) {
        values <- output[[key]][[name]]
        if (is.null(values)) return(rep(NA_real_, length(labels)))
        if (length(values) != length(labels)) stop("Posterior summary dimensions do not match for ", name)
        as.vector(values)
      }
      results[[name]] <- data.frame(
        parameter = labels, mean = as.vector(values), sd = component("sd"),
        q2.5 = component("q2.5"), q50 = component("q50"), q97.5 = component("q97.5"),
        row.names = labels)
    }
  }
  if (parameter != "all") return(results[[parameter]])
  if (fitted_model && !is.null(fit$group_membership)) results$group_membership <- fit$group_membership
  results$mean <- output$mean
  results$sd <- output$sd
  results$quantiles <- list(q2.5 = output$q2.5, q50 = output$q50, q97.5 = output$q97.5)
  results
}
