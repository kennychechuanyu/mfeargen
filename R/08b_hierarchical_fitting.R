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
#' @param observation_model `"normal"` for Gaussian ratings or `"censored"`
#'   for ratings clipped to instrument bounds. NULL uses simulation metadata
#'   when available and otherwise selects `"normal"`.
#' @param save_predictive Store posterior predictive draws during fitting.
#'   FALSE saves memory; prediction draws can be generated from fitted parameters.
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
                                     observation_model = NULL,
                                     save_predictive = FALSE,
                                     ...) {
  model_type <- match.arg(model_type)
  if (is.null(observation_model)) observation_model <- attr(data, "observation_model") %||% "normal"
  observation_model <- match.arg(observation_model, c("normal", "censored"))
  if (model_type != "1cs") stop("The repeated-session model currently supports only model_type = '1cs'")
  if (!is.null(priors)) stop("Custom priors are not supported by the repeated-session model")
  validate_scalar(n.iter, "n.iter", lower = 1, integer = TRUE)
  validate_scalar(n.burnin, "n.burnin", lower = 0, integer = TRUE)
  validate_scalar(n.chains, "n.chains", lower = 1, integer = TRUE)
  validate_scalar(n.thin, "n.thin", lower = 1, integer = TRUE)
  if (n.iter <= n.burnin || floor((n.iter - n.burnin) / n.thin) < 2) {
    stop("n.iter, n.burnin, and n.thin must retain at least two draws per chain")
  }
  for (argument in c("parallel", "verbose", "save_predictive")) {
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

  require_jags_backend()
  jags_data <- prepare_hierarchical_data(data, model_type, lambda_criterion, rating_scale,
                                         observation_model = observation_model)
  model_file <- system.file("jags", "model_1cs_hierarchical.txt", package = "mfeargen")
  if (!nzchar(model_file)) stop("The repeated-session JAGS model is missing from the package installation")
  compiled <- compile_observation_model(readLines(model_file, warn = FALSE), jags_data,
                                        observation_model)
  if (observation_model == "censored") {
    model_file <- tempfile("mfeargen-repeated-", fileext = ".txt")
    writeLines(compiled$model_code, model_file)
    on.exit(unlink(model_file), add = TRUE)
  }
  parameters <- c(
    paste0("mu_", c("alpha", "lambda", "w0", "w1")),
    paste0("tau_", c("alpha", "lambda", "w0", "w1"), "_sq"),
    paste0("sigma_", c("alpha", "lambda", "w0", "w1"), "_sq"),
    paste0("stability_", c("alpha", "lambda", "w0", "w1")),
    paste0(c("alpha", "lambda", "w0", "w1"), "_stable"),
    "alpha", "lambda", "w0", "w1", "gp", "pi", "sigma_obs",
    "theta", "loglik"
  )
  if (save_predictive) parameters <- c(parameters, "y_pre")
  started <- Sys.time()
  jags_fit <- jagsUI::jags(
    data = compiled$data,
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
  jags_fit <- normalize_pointwise_output(jags_fit, dim(jags_data$y))
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
    n.iter = n.iter, n.burnin = n.burnin, n.chains = n.chains,
    n.thin = n.thin, seed = seed, observation_model = observation_model,
    save_predictive = save_predictive,
    call = match.call()
  )
  class(results) <- c("mfeargen_hierarchical", "list")
  results$variance_components <- calculate_variance_decomposition(results)
  results
}

#' Prepare Repeated-Session Data
#' @keywords internal
prepare_hierarchical_data <- function(data, model_type = "1cs", lambda_criterion = 0.7,
                                      rating_scale = NULL, observation_model = NULL) {
  if (!is.data.frame(data) || !nrow(data)) stop("data must be a nonempty data frame")
  missing <- setdiff(c("participant", "session", "trials", "US_expect"), names(data))
  if (length(missing)) stop("Missing hierarchical data columns: ", paste(missing, collapse = ", "))
  if (model_type != "1cs") stop("Only the single-cue repeated-session model is supported")
  for (column in c("participant", "session")) {
    values <- data[[column]]
    if (!(is.numeric(values) || is.character(values) || is.factor(values)) ||
        anyNA(values) || any(!nzchar(as.character(values))) ||
        (is.numeric(values) && any(!is.finite(values)))) {
      stop(column, " must contain non-missing numeric or text identifiers")
    }
  }

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
                                    rating_scale = rating_scale,
                                    observation_model = observation_model)
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
  row_index <- array(NA_integer_, dims)
  for (i in seq_along(participants)) {
    for (s in seq_along(sessions[[i]])) {
      row <- session_rows[session_map[i, s]]
      result$Ntrials[i, s] <- counts[row]
      result$y[i, s, ] <- prepared$y[row, ]
      result$k[i, s, ] <- prepared$k[row, ]
      result$r[i, s, ] <- prepared$r[row, ]
      result$d_phy[i, s, ] <- prepared$d_phy[row, ]
      result$d_per[i, s, ] <- prepared$d_per_old[row, ]
      row_index[i, s, ] <- attr(prepared, "row_index")[row, ]
    }
  }
  attr(result, "row_index") <- row_index
  attr(result, "observed") <- !is.na(result$y)
  attr(result, "observation_model") <- attr(prepared, "observation_model") %||% "normal"
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
  result <- check_convergence(jags_fit, verbose = FALSE)
  result$n_converged <- sum(result$diagnostics$status == "passed")
  result$percent_converged <- if (result$n_parameters) {
    100 * result$n_converged / result$n_parameters
  } else 0
  result
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
  cat("Convergence:", if (isTRUE(x$convergence$converged)) "passed" else "requires review", "\n")
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
    coefficients = stats::coef(object),
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

hierarchical_response_array <- function(object, values) {
  dimensions <- dim(object$jags_data$y)
  if (length(values) != prod(dimensions)) stop("Unexpected response dimensions in model output")
  result <- array(values, dimensions)
  row_index <- attr(object$jags_data, "row_index")
  if (!is.null(row_index)) result[is.na(row_index)] <- NA_real_
  dimnames(result) <- list(as.character(object$data_info$participants),
                          as.character(seq_len(dimensions[2])), NULL)
  attr(result, "sessions") <- object$data_info$sessions
  result
}

hierarchical_sigma_draws <- function(object, n_draws) {
  sigma <- object$jags_fit$sims.list$sigma_obs
  groups <- object$jags_fit$sims.list$gp
  if (is.null(sigma) || is.null(groups)) stop("Observation noise draws are unavailable; refit the model")
  sigma <- matrix(sigma, ncol = 2L)
  groups <- matrix(groups, ncol = object$jags_data$Nparticipants)
  if (nrow(sigma) != n_draws || nrow(groups) != n_draws) stop("Posterior draws have inconsistent dimensions")
  result <- matrix(NA_real_, n_draws, ncol(groups))
  for (i in seq_len(ncol(groups))) {
    result[, i] <- sigma[cbind(seq_len(n_draws), ifelse(groups[, i] == 1L, 2L, 1L))]
  }
  result
}

hierarchical_response_draws <- function(object, predictive = FALSE, scale = "response", seed = NULL) {
  node <- if (predictive) "y_pre" else "theta"
  draws <- object$jags_fit$sims.list[[node]]
  generate <- predictive && is.null(draws)
  if (generate) draws <- object$jags_fit$sims.list$theta
  expected <- dim(object$jags_data$y)
  if (is.null(draws)) stop("Posterior response draws are unavailable; refit the model")
  if (length(dim(draws)) != 4L || !identical(as.integer(dim(draws)[-1]), as.integer(expected))) {
    stop("Unexpected posterior dimensions for ", node)
  }
  if (generate) {
    noise <- hierarchical_sigma_draws(object, dim(draws)[1])
    draws <- with_prediction_seed(seed, {
      for (i in seq_len(expected[1])) {
        values <- as.vector(draws[, i, , ])
        deviations <- rep(noise[, i], length.out = length(values))
        present <- is.finite(values)
        values[present] <- stats::rnorm(sum(present), values[present], deviations[present])
        draws[, i, , ] <- values
      }
      draws
    })
  }
  if (response_model(object) == "censored") {
    if (predictive) draws <- pmin(pmax(draws, 0), 10) else {
      sigma <- hierarchical_sigma_draws(object, dim(draws)[1])
      for (i in seq_len(expected[1])) {
        draws[, i, , ] <- clipped_normal_mean(matrix(draws[, i, , ], nrow = dim(draws)[1]), sigma[, i])
      }
    }
  }
  row_index <- attr(object$jags_data, "row_index")
  flat <- matrix(draws, nrow = dim(draws)[1])
  flat[, is.na(as.vector(row_index))] <- NA_real_
  dim(flat) <- dim(draws)
  if (scale == "response") flat <- object$rating_scale[1] + flat * diff(object$rating_scale) / 10
  flat
}

#' Expected Responses for a Repeated-Session Model
#' @param object An `mfeargen_hierarchical` object.
#' @param scale `"response"` for original rating units, or `"normalized"` for 0-10.
#' @param ... Unused.
#' @return A participant by session by trial array. Absent sessions and padded
#'   trials are `NA`; the `sessions` attribute records original session IDs.
#' @export
fitted.mfeargen_hierarchical <- function(object, scale = c("response", "normalized"), ...) {
  scale <- match.arg(scale)
  if (response_model(object) == "censored") {
    draws <- hierarchical_response_draws(object, scale = scale)
    return(hierarchical_response_array(object, colMeans(matrix(draws, nrow = dim(draws)[1]))))
  }
  theta <- object$jags_fit$mean$theta
  if (is.null(theta)) stop("Expected responses were not monitored; refit the model")
  result <- hierarchical_response_array(object, theta)
  if (scale == "response") result <- result / 10 * diff(object$rating_scale) + object$rating_scale[1]
  result
}

#' Residuals for a Repeated-Session Model
#' @param object An `mfeargen_hierarchical` object.
#' @param type `"raw"` in original rating units or `"standardized"` by posterior
#'   mean conditional observation standard deviation.
#' @param ... Unused.
#' @return A participant by session by trial array; missing ratings are `NA`.
#' @export
residuals.mfeargen_hierarchical <- function(object, type = c("raw", "standardized"), ...) {
  type <- match.arg(type)
  result <- object$jags_data$y - stats::fitted(object, scale = "normalized")
  if (type == "raw") return(result * diff(object$rating_scale) / 10)
  theta <- object$jags_fit$sims.list$theta
  if (is.null(theta)) stop("Expected response draws are unavailable; refit the model")
  noise <- hierarchical_sigma_draws(object, dim(theta)[1])
  if (response_model(object) == "censored") {
    for (i in seq_len(object$jags_data$Nparticipants)) {
      theta[, i, , ] <- clipped_normal_sd(matrix(theta[, i, , ], nrow = dim(theta)[1]), noise[, i])
    }
    denominator <- hierarchical_response_array(object, colMeans(matrix(theta, nrow = dim(theta)[1])))
    return(result / denominator)
  }
  sweep(result, 1, colMeans(noise), "/")
}

#' Predictions for Repeated Sessions
#' @param object An `mfeargen_hierarchical` object.
#' @param newdata Reserved; must be `NULL`.
#' @param type `"groups"` or `"response"`.
#' @param interval `"none"`, `"credible"` for expected responses, or `"prediction"`
#'   for future observations under the same fitted trial histories.
#' @param level Interval probability between zero and one.
#' @param format `"array"` or `"data.frame"` in original input row order.
#' @param scale `"response"` for original rating units or `"normalized"` for 0-10.
#' @param seed Optional seed for generated predictive draws; preserves the caller's RNG state.
#' @param ... Arguments passed to [predict_groups()] when `type = "groups"`.
#' @return Group assignments or expected responses in original rating units.
#'   Array output has participant, session, and trial dimensions; intervals add
#'   a final dimension for estimate, lower, and upper. Prediction intervals from
#'   the normal observation model can extend beyond the instrument bounds;
#'   censored-model prediction intervals remain within them.
#' @export
predict.mfeargen_hierarchical <- function(object, newdata = NULL, type = c("groups", "response"),
                                           interval = c("none", "credible", "prediction"), level = 0.95,
                                           format = c("array", "data.frame"),
                                           scale = c("response", "normalized"), seed = NULL, ...) {
  if (!is.null(newdata)) stop("Prediction for newdata is unavailable; omit newdata to predict fitted sessions")
  type <- match.arg(type)
  if (type == "groups") return(predict_groups(object, ...))
  interval <- match.arg(interval)
  format <- match.arg(format)
  scale <- match.arg(scale)
  validate_scalar(level, "level", lower = 0, upper = 1, lower_open = TRUE, upper_open = TRUE)
  estimate <- stats::fitted(object, scale = scale)
  row_index <- attr(object$jags_data, "row_index")
  if (is.null(row_index)) stop("Trial mappings are unavailable; refit the model to obtain predictions")
  lower <- upper <- NULL
  if (interval != "none") {
    draws <- hierarchical_response_draws(object, predictive = interval == "prediction", scale = scale, seed = seed)
    cells <- length(estimate)
    draws <- matrix(draws, ncol = cells)
    present <- !is.na(row_index)
    limits <- apply(draws[, present, drop = FALSE], 2, stats::quantile,
                    probs = c((1 - level) / 2, (1 + level) / 2), names = FALSE)
    lower <- upper <- rep(NA_real_, cells)
    lower[present] <- limits[1, ]
    upper[present] <- limits[2, ]
    lower <- hierarchical_response_array(object, lower)
    upper <- hierarchical_response_array(object, upper)
  }
  if (format == "data.frame") {
    values <- rep(NA_real_, nrow(object$data))
    present <- !is.na(row_index)
    values[row_index[present]] <- estimate[present]
    result <- data.frame(participant = object$data$participant, session = object$data$session,
                         trials = object$data$trials, observed = object$data$US_expect,
                         estimate = values)
    if (scale == "normalized") {
      result$observed <- (result$observed - object$rating_scale[1]) / diff(object$rating_scale) * 10
    }
    if (!is.null(lower)) {
      result$lower <- result$upper <- NA_real_
      result$lower[row_index[present]] <- lower[present]
      result$upper[row_index[present]] <- upper[present]
    }
    return(result)
  }
  if (is.null(lower)) return(estimate)
  result <- array(c(estimate, lower, upper), c(dim(estimate), 3L),
                  dimnames = c(dimnames(estimate), list(c("estimate", "lower", "upper"))))
  attr(result, "sessions") <- object$data_info$sessions
  result
}

#' Coefficients for a Repeated-Session Model
#' @param object An `mfeargen_hierarchical` object.
#' @param ... Unused.
#' @param scale `"response"` converts observation noise to original rating units;
#'   `"normalized"` retains the internal 0-10 units. Other parameters retain
#'   their latent learning, response, and stimulus-distance units.
#' @return A data frame of posterior means, standard deviations, and 95 percent
#'   credible intervals. Participant and session IDs accompany indexed individual
#'   parameters; population and variance parameters have missing IDs.
#' @export
coef.mfeargen_hierarchical <- function(object, scale = c("response", "normalized"), ...) {
  scale <- match.arg(scale)
  summaries <- object$jags_fit$summary
  keep <- !grepl("^(theta|loglik|y_pre|gp|deviance)(\\[|$)", rownames(summaries)) &
    is.finite(summaries[, "mean"])
  summaries <- summaries[keep, , drop = FALSE]
  result <- data.frame(parameter = rownames(summaries), participant = NA_character_,
                        session = NA_character_, estimate = summaries[, "mean"],
                        posterior_sd = summaries[, "sd"], ci_lower = summaries[, "2.5%"],
                        ci_upper = summaries[, "97.5%"], row.names = NULL)
  families <- sub("\\[.*$", "", result$parameter)
  result$component <- ifelse(families == "sigma_obs", "observation_noise",
    ifelse(families == "pi", "group_probability",
    ifelse(grepl("^mu_", families), "population",
    ifelse(grepl("^(tau_|sigma_)", families), "parameter_variance",
    ifelse(grepl("^stability_", families), "variance_ratio",
    ifelse(grepl("_stable$", families), "participant", "session"))))))
  if (scale == "response") {
    noise <- families == "sigma_obs"
    result[noise, c("estimate", "posterior_sd", "ci_lower", "ci_upper")] <-
      result[noise, c("estimate", "posterior_sd", "ci_lower", "ci_upper")] * diff(object$rating_scale) / 10
  }
  for (row in seq_len(nrow(result))) {
    name <- result$parameter[row]
    if (name %in% paste0(c("alpha", "lambda", "w0", "w1"), "_stable") &&
        object$data_info$n_participants == 1L) {
      result$participant[row] <- as.character(object$data_info$participants[1])
      next
    }
    if (!grepl("^(alpha|lambda|w0|w1)(_stable)?\\[", name)) next
    indices <- as.integer(strsplit(sub(".*\\[([^]]+)\\]$", "\\1", name), ",", fixed = TRUE)[[1]])
    result$participant[row] <- as.character(object$data_info$participants[indices[1]])
    if (length(indices) == 2L && indices[2] <= length(object$data_info$sessions[[indices[1]]])) {
      result$session[row] <- as.character(object$data_info$sessions[[indices[1]]][indices[2]])
    }
  }
  result
}
