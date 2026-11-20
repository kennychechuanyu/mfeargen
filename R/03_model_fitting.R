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
                        ...) {
  model_type <- match.arg(model_type)
  if (!is.data.frame(data)) stop("data must be a data frame")
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
  for (nm in c("parallel", "verbose")) {
    value <- get(nm)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(nm, " must be TRUE or FALSE")
    }
  }
  if (!is.null(priors)) validate_priors(priors)
  dots <- list(...)
  reserved <- intersect(names(dots), c("data", "inits", "parameters.to.save", "model.file"))
  if (length(reserved)) stop("Reserved JAGS arguments: ", paste(reserved, collapse = ", "))
  if (isTRUE(dots$codaOnly)) stop("codaOnly = TRUE is incompatible with mfeargen post-processing")
  dic <- dots$DIC %||% TRUE
  dots$DIC <- NULL
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
      data, model_type, lambda_criterion, rating_scale = rating_scale
    )
  } else {
    jags_data <- prepared_data
    validate_prepared_data(jags_data, model_type)
    if (!is.null(rating_scale) &&
        !isTRUE(all.equal(as.numeric(rating_scale), as.numeric(attr(jags_data, "rating_scale") %||% jags_data$rating_scale)))) {
      stop("rating_scale does not match prepared_data$rating_scale")
    }
  }
  validate_prepared_data(jags_data, model_type)

  model_file <- system.file("jags", paste0("model_", model_type, ".txt"),
                            package = "mfeargen")
  if (!nzchar(model_file) || !file.exists(model_file)) {
    stop("JAGS model file not found; reinstall mfeargen")
  }
  model_text <- readLines(model_file, warn = FALSE)
  if (!is.null(priors)) {
    model_text <- apply_model_priors(model_text, priors)
    model_file <- tempfile("mfeargen-", fileext = ".txt")
    writeLines(model_text, model_file)
    on.exit(unlink(model_file), add = TRUE)
  }
  inits <- initialize_mcmc_parameters(jags_data, n.chains)
  parameters <- c("alpha", "lambda", "gp", "pi", "w0", "w1", "sigma",
                  "theta", "loglik", "y_pre", "alpha_mu", "alpha_kappa",
                  "lambda_mu", "lambda_sigma", "w0_mu", "w0_sigma", "w1_a", "w1_b")
  # Prior metadata belongs to the R object, not JAGS' data environment.
  sampling_data <- jags_data[!grepl("^prior_", names(jags_data))]
  if (verbose) message("Running MCMC sampling...")
  fit <- do.call(jagsUI::jags, c(list(
    data = sampling_data, inits = inits, parameters.to.save = parameters,
    model.file = model_file, n.chains = n.chains, n.iter = n.iter,
    n.burnin = n.burnin, n.thin = n.thin, parallel = parallel,
    DIC = dic, verbose = verbose
  ), dots))
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
    model_code = model_text, call = match.call()
  )
  class(result) <- "mfeargen"
  if (verbose) message("Model fitting complete (", jags_data$Nparticipants, " participants).")
  result
}

validate_prepared_data <- function(jags_data, model_type) {
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
  invisible(TRUE)
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


#' Check MCMC Convergence
#'
#' @description
#' Checks convergence of MCMC chains using Rhat statistics.
#'
#' @param fit JAGS output object
#' @param threshold Numeric. Rhat threshold for convergence (default: 1.1)
#' @param verbose Logical. Print diagnostic messages
#' @return List with convergence information
#' @export
check_convergence <- function(fit, threshold = 1.1, verbose = TRUE) {
  validate_scalar(threshold, "threshold", 1)
  if (inherits(fit, "mfeargen")) fit <- fit$jags_fit
  if (is.numeric(fit)) {
    rhat <- fit
  } else if (!is.null(fit$Rhat)) {
    rhat <- unlist(fit$Rhat)
  } else if (!is.null(fit$summary) && "Rhat" %in% colnames(fit$summary)) {
    rhat <- fit$summary[, "Rhat"]
  } else if (!is.null(fit$BUGSoutput$summary) &&
             "Rhat" %in% colnames(fit$BUGSoutput$summary)) {
    rhat <- fit$BUGSoutput$summary[, "Rhat"]
  } else {
    rhat <- numeric(0)
  }
  rhat <- rhat[!is.na(rhat)]
  result <- list(
    converged = if (length(rhat)) all(rhat < threshold) else NA,
    max_rhat = if (length(rhat)) max(rhat) else NA_real_,
    n_parameters = length(rhat), n_poor_convergence = sum(rhat >= threshold)
  )
  if (verbose && identical(result$converged, FALSE)) {
    message(result$n_poor_convergence, " parameters have Rhat >= ", threshold,
            "; maximum Rhat = ", round(result$max_rhat, 3), ".")
  }
  if (verbose && is.na(result$converged)) {
    message("Rhat is unavailable; convergence cannot be assessed.")
  }
  result
}


#' Extract Model Parameters
#'
#' @description
#' Extracts posterior distributions of model parameters.
#'
#' @param fit mfeargen object or JAGS output
#' @param parameter Character. Which parameter to extract
#' @param summary_only Logical. Return summary statistics only
#' @return Data frame or matrix of parameter values
#' @export
#' @examples
#' \dontrun{
#' fit <- fit_mfeargen(data)
#' alpha_params <- extract_parameters(fit, "alpha")
#' }
extract_parameters <- function(fit, parameter = "all", summary_only = TRUE) {
  if (!is.character(parameter) || length(parameter) != 1L || is.na(parameter)) {
    stop("parameter must be one parameter name or 'all'")
  }
  
  # Handle mfeargen object
  if (inherits(fit, "mfeargen")) {
    jags_output <- fit$jags_fit
  } else {
    jags_output <- fit
  }
  
  # Available parameters - check what exists in the output
  if (!is.null(jags_output$mean)) {
    available_params <- names(jags_output$mean)
  } else if (!is.null(jags_output$summary)) {
    # Extract parameter names from summary
    param_names <- rownames(jags_output$summary)
    available_params <- unique(gsub("\\[.*\\]", "", param_names))
  } else {
    stop("Cannot extract parameters: model output format not recognized")
  }
  
  if (parameter == "all") {
    parameters <- available_params
  } else {
    if (!(parameter %in% available_params)) {
      stop("Parameter '", parameter, "' not found in model output")
    }
    parameters <- parameter
  }
  
  # Extract parameters
  results <- list()
  
  for (param in parameters) {
    if (summary_only) {
      # Return summary statistics
      results[[param]] <- data.frame(
        mean = as.vector(jags_output$mean[[param]]),
        sd = as.vector(jags_output$sd[[param]]),
        q2.5 = as.vector(jags_output$q2.5[[param]]),
        q50 = as.vector(jags_output$q50[[param]]),
        q97.5 = as.vector(jags_output$q97.5[[param]])
      )
    } else {
      # Return full posterior samples
      if (param %in% names(jags_output$sims.list)) {
        results[[param]] <- jags_output$sims.list[[param]]
      }
    }
  }
  
  # Add group membership if available and requested
  if (parameter == "all" && inherits(fit, "mfeargen") && !is.null(fit$group_membership)) {
    results$group_membership <- fit$group_membership
  }
  
  if (length(results) == 1) {
    return(results[[1]])
  } else {
    # For multiple parameters, add summary components for backwards compatibility
    results$mean <- jags_output$mean
    results$sd <- jags_output$sd
    results$quantiles <- list(
      q2.5 = jags_output$q2.5,
      q50 = jags_output$q50,
      q97.5 = jags_output$q97.5
    )
    return(results)
  }
}
