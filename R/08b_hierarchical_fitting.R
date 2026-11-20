#' Fit a Repeated-Session Fear Generalization Model
#'
#' Fits a single-cue mixture model with population parameters, stable participant
#' components, and session-specific parameters. Each session starts with zero
#' associative strength. Ratings are modeled on an internal 0-10 scale.
#'
#' @param data Trial-level data accepted by [prepare_mfeargen_data()], with an
#'   additional `session` column. Trial numbers must be unique within each
#'   participant and session.
#' @param model_type Only `"1cs"` is currently supported.
#' @param n.iter Total iterations per chain, including burn-in.
#' @param n.burnin Burn-in iterations per chain.
#' @param n.chains Number of MCMC chains.
#' @param n.thin Thinning interval.
#' @param lambda_criterion Similarity at the overgeneralization boundary.
#' @param priors Reserved for future use. Custom hierarchical priors are not
#'   supported; a non-NULL value produces an error.
#' @param parallel Run chains in parallel.
#' @param verbose Print sampler progress and diagnostics.
#' @param rating_scale Numeric vector containing the original rating bounds.
#' @param seed Optional integer seed for reproducible initialization and sampling.
#' @param ... Additional arguments passed to [jagsUI::jags()].
#' @return An `mfeargen_hierarchical` object containing posterior samples,
#'   participant and session mappings, convergence diagnostics, and variance
#'   summaries. Variance ratios describe the underlying normal distributions;
#'   truncation and group-specific parameter constraints affect realized variance.
#' @export
#' @examples
#' \dontrun{
#' first <- simulate_mfeargen(n_participants = 6, seed = 1)
#' second <- simulate_mfeargen(n_participants = 6, seed = 2)
#' first$session <- "first"
#' second$session <- "second"
#' repeated <- rbind(first, second)
#' fit <- fit_mfeargen_hierarchical(repeated, rating_scale = c(1, 10))
#' calculate_variance_decomposition(fit)
#' plot_hierarchical_structure(fit)
#' }
fit_mfeargen_hierarchical <- function(data,
                                     model_type = c("1cs", "2cs"),
                                     n.iter = 20000,
                                     n.burnin = 10000,
                                     n.chains = 3,
                                     n.thin = 5,
                                     lambda_criterion = 0.7,
                                     priors = NULL,
                                     parallel = TRUE,
                                     verbose = TRUE,
                                     rating_scale = NULL,
                                     seed = NULL,
                                     ...) {
  model_type <- match.arg(model_type)
  if (model_type != "1cs") stop("The repeated-session model currently supports only model_type = '1cs'")
  if (!is.null(priors)) stop("Custom priors are not supported by the repeated-session model")
  validate_scalar(n.iter, "n.iter", lower = 1, integer = TRUE)
  validate_scalar(n.burnin, "n.burnin", lower = 0, integer = TRUE)
  validate_scalar(n.chains, "n.chains", lower = 1, integer = TRUE)
  validate_scalar(n.thin, "n.thin", lower = 1, integer = TRUE)
  if (n.iter <= n.burnin || floor((n.iter - n.burnin) / n.thin) < 2) {
    stop("n.iter, n.burnin, and n.thin must retain at least two draws per chain")
  }
  for (argument in c("parallel", "verbose")) {
    value <- get(argument)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(argument, " must be TRUE or FALSE")
    }
  }
  if (!is.null(seed)) {
    validate_scalar(seed, "seed", lower = 0, upper = .Machine$integer.max, integer = TRUE)
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

  jags_data <- prepare_hierarchical_data(data, model_type, lambda_criterion, rating_scale)
  model_file <- system.file("jags", "model_1cs_hierarchical.txt", package = "mfeargen")
  if (!nzchar(model_file)) stop("The repeated-session JAGS model is missing from the package installation")
  parameters <- c(
    paste0("mu_", c("alpha", "lambda", "w0", "w1")),
    paste0("tau_", c("alpha", "lambda", "w0", "w1"), "_sq"),
    paste0("sigma_", c("alpha", "lambda", "w0", "w1"), "_sq"),
    paste0("stability_", c("alpha", "lambda", "w0", "w1")),
    paste0(c("alpha", "lambda", "w0", "w1"), "_stable"),
    "alpha", "lambda", "w0", "w1", "gp", "pi", "sigma_obs"
  )
  started <- Sys.time()
  jags_fit <- jagsUI::jags(
    data = jags_data,
    inits = initialize_hierarchical_parameters(jags_data, n.chains),
    parameters.to.save = parameters,
    model.file = model_file,
    n.iter = n.iter,
    n.burnin = n.burnin,
    n.chains = n.chains,
    n.thin = n.thin,
    parallel = parallel,
    verbose = verbose,
    ...
  )
  membership <- extract_group_membership(jags_fit)
  membership$participant <- attr(jags_data, "participants")
  results <- list(
    jags_fit = jags_fit,
    model_type = model_type,
    data = data,
    jags_data = jags_data,
    rating_scale = attr(jags_data, "rating_scale"),
    group_membership = membership,
    convergence = check_hierarchical_convergence(jags_fit),
    data_info = list(
      n_participants = jags_data$Nparticipants,
      total_sessions = sum(jags_data$Nsessions),
      total_trials = sum(jags_data$Ntrials),
      participants = attr(jags_data, "participants"),
      sessions = attr(jags_data, "sessions")
    ),
    population_parameters = extract_population_parameters(jags_fit),
    runtime = as.numeric(difftime(Sys.time(), started, units = "mins")),
    call = match.call()
  )
  class(results) <- c("mfeargen_hierarchical", "list")
  results$variance_components <- calculate_variance_decomposition(results)
  results
}

#' Prepare Repeated-Session Data
#' @keywords internal
prepare_hierarchical_data <- function(data, model_type = "1cs", lambda_criterion = 0.7,
                                      rating_scale = NULL) {
  if (!is.data.frame(data) || !nrow(data)) stop("data must be a nonempty data frame")
  missing <- setdiff(c("participant", "session", "trials", "US_expect"), names(data))
  if (length(missing)) stop("Missing hierarchical data columns: ", paste(missing, collapse = ", "))
  if (model_type != "1cs") stop("Only the single-cue repeated-session model is supported")
  if (anyNA(data$participant) || anyNA(data$session)) stop("participant and session IDs cannot be missing")

  participants <- unique(data$participant)
  sessions <- lapply(participants, function(id) unique(data$session[data$participant == id]))
  nsessions <- lengths(sessions)
  if (all(nsessions < 2)) stop("At least one participant must have repeated sessions")
  if (is.null(rating_scale)) rating_scale <- attr(data, "rating_scale")
  if (is.null(rating_scale)) rating_scale <- detect_rating_scale(data$US_expect)
  validate_rating_scale(rating_scale)

  session_data <- data
  session_data$participant <- integer(nrow(data))
  session_map <- matrix(0L, nrow = length(participants), ncol = max(nsessions))
  index <- 0L
  for (i in seq_along(participants)) {
    for (s in seq_along(sessions[[i]])) {
      index <- index + 1L
      rows <- data$participant == participants[i] & data$session == sessions[[i]][s]
      session_data$participant[rows] <- index
      session_map[i, s] <- index
    }
  }
  prepared <- prepare_mfeargen_data(session_data, "1cs", lambda_criterion,
                                    rating_scale = rating_scale)
  # Data preparation preserves first appearance, which can differ from the nested map.
  session_rows <- match(seq_len(index), attr(prepared, "participants"))
  counts <- attr(prepared, "trial_counts")
  dims <- c(length(participants), max(nsessions), prepared$Ntrials)
  result <- list(
    Nparticipants = length(participants),
    Nsessions = nsessions,
    Ntrials = matrix(0L, dims[1], dims[2]),
    y = array(NA_real_, dims),
    k = array(0, dims),
    r = array(0, dims),
    d_phy = array(0, dims),
    d_per = array(0, dims),
    lambda_threshold = prepared$lambda_threshold,
    A = 0,
    K = 10,
    pi_alpha = rep(1, 4)
  )
  for (i in seq_along(participants)) {
    for (s in seq_along(sessions[[i]])) {
      row <- session_rows[session_map[i, s]]
      result$Ntrials[i, s] <- counts[row]
      result$y[i, s, ] <- prepared$y[row, ]
      result$k[i, s, ] <- prepared$k[row, ]
      result$r[i, s, ] <- prepared$r[row, ]
      result$d_phy[i, s, ] <- prepared$d_phy[row, ]
      result$d_per[i, s, ] <- prepared$d_per_old[row, ]
    }
  }
  attr(result, "participants") <- participants
  attr(result, "sessions") <- sessions
  attr(result, "rating_scale") <- rating_scale
  result
}

#' Initialize Repeated-Session Chains
#' @keywords internal
initialize_hierarchical_parameters <- function(jags_data, n.chains) {
  lapply(seq_len(n.chains), function(chain) {
    list(
      gp = sample.int(4, jags_data$Nparticipants, replace = TRUE),
      pi = rep(0.25, 4),
      alpha_stable_raw = runif(jags_data$Nparticipants, 0.2, 0.8),
      lambda_stable_raw = rep(jags_data$lambda_threshold * 2, jags_data$Nparticipants),
      w0_stable = rnorm(jags_data$Nparticipants, -1, 0.5),
      w1_stable_raw = runif(jags_data$Nparticipants, 1, 3),
      .RNG.name = "base::Mersenne-Twister",
      .RNG.seed = sample.int(.Machine$integer.max, 1)
    )
  })
}

#' Check Repeated-Session Convergence
#' @keywords internal
check_hierarchical_convergence <- function(jags_fit) {
  rhat <- unlist(jags_fit$Rhat, use.names = TRUE)
  rhat <- rhat[is.finite(rhat)]
  n_parameters <- length(rhat)
  n_converged <- sum(rhat < 1.1)
  list(
    converged = n_parameters > 0 && n_converged == n_parameters,
    n_parameters = n_parameters,
    n_converged = n_converged,
    percent_converged = if (n_parameters) 100 * n_converged / n_parameters else 0,
    max_rhat = if (n_parameters) max(rhat) else NA_real_,
    problematic_parameters = names(rhat)[rhat >= 1.1]
  )
}

#' Extract Population Means
#' @keywords internal
extract_population_parameters <- function(jags_fit) {
  names <- paste0("mu_", c("alpha", "lambda", "w0", "w1"))
  stats::setNames(lapply(names, function(name) mean(jags_fit$sims.list[[name]])), names)
}

#' Print a Repeated-Session Model
#' @param x An `mfeargen_hierarchical` object.
#' @param ... Unused.
#' @export
print.mfeargen_hierarchical <- function(x, ...) {
  cat("Repeated-session mfeargen model (1cs)\n")
  cat(sprintf("Participants: %d; sessions: %d; trials: %d\n",
              x$data_info$n_participants, x$data_info$total_sessions,
              x$data_info$total_trials))
  cat(sprintf("Maximum Rhat: %.3f\n", x$convergence$max_rhat))
  print(x$variance_components)
  invisible(x)
}

#' Summarize a Repeated-Session Model
#' @param object An `mfeargen_hierarchical` object.
#' @param ... Unused.
#' @return A list of population summaries, variance components, group membership,
#'   and convergence diagnostics.
#' @export
summary.mfeargen_hierarchical <- function(object, ...) {
  list(
    population_parameters = object$population_parameters,
    variance_components = object$variance_components,
    group_membership = object$group_membership,
    convergence = object$convergence
  )
}

#' Plot a Repeated-Session Model
#' @param x An `mfeargen_hierarchical` object.
#' @param parameter Parameter to plot: `"alpha"`, `"lambda"`, `"w0"`, or `"w1"`.
#' @param ... Passed to [plot_hierarchical_structure()].
#' @return A ggplot object.
#' @export
plot.mfeargen_hierarchical <- function(x, parameter = "alpha", ...) {
  plot_hierarchical_structure(x, parameter = parameter, ...)
}
