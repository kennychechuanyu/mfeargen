response_model <- function(object) {
  object$observation_model %||% attr(object$jags_data, "observation_model") %||% "normal"
}

response_bounds <- function(object) {
  bounds <- object$rating_scale %||% attr(object$jags_data, "rating_scale") %||% c(0, 10)
  validate_rating_scale(bounds)
  bounds
}

clipped_normal_mean <- function(location, sigma, lower = 0, upper = 10) {
  a <- (lower - location) / sigma
  b <- (upper - location) / sigma
  middle <- stats::pnorm(b) - stats::pnorm(a)
  lower * stats::pnorm(a) + location * middle +
    sigma * (stats::dnorm(a) - stats::dnorm(b)) +
    upper * stats::pnorm(b, lower.tail = FALSE)
}

clipped_normal_sd <- function(location, sigma, lower = 0, upper = 10) {
  a <- (lower - location) / sigma
  b <- (upper - location) / sigma
  pa <- stats::pnorm(a)
  pb <- stats::pnorm(b, lower.tail = FALSE)
  middle <- stats::pnorm(b) - pa
  da <- stats::dnorm(a)
  db <- stats::dnorm(b)
  a_da <- ifelse(is.finite(a), a * da, 0)
  b_db <- ifelse(is.finite(b), b * db, 0)
  centered_mean <- (lower - location) * pa + (upper - location) * pb + sigma * (da - db)
  centered_second <- (lower - location)^2 * pa + (upper - location)^2 * pb +
    sigma^2 * (middle + a_da - b_db)
  variance <- centered_second - centered_mean^2
  broad <- (location * 0 + sigma) > 10 * (upper - lower)
  if (any(broad)) {
    locations <- (location + sigma * 0)[broad]
    deviations <- (sigma + location * 0)[broad]
    # Integrating the bounded interval avoids cancellation when noise is very wide.
    nodes <- c(-0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
               -0.1834346424956498, 0.1834346424956498, 0.5255324099163290,
               0.7966664774136267, 0.9602898564975363)
    weights <- c(0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
                 0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
                 0.2223810344533745, 0.1012285362903763)
    half_width <- (upper - lower) / 2
    center <- (upper + lower) / 2
    first <- second <- numeric(length(locations))
    for (k in seq_along(nodes)) {
      value <- center + half_width * nodes[k]
      mass <- half_width * weights[k] * stats::dnorm(value, locations, deviations)
      first <- first + (value - locations) * mass
      second <- second + (value - locations)^2 * mass
    }
    left <- stats::pnorm(lower, locations, deviations)
    right <- stats::pnorm(upper, locations, deviations, lower.tail = FALSE)
    first <- first + (lower - locations) * left + (upper - locations) * right
    second <- second + (lower - locations)^2 * left + (upper - locations)^2 * right
    variance[broad] <- second - first^2
  }
  sqrt(pmin(pmax(variance, 0), (upper - lower)^2 / 4))
}

observation_sigma_draws <- function(object, n_draws) {
  sigma <- object$jags_fit$sims.list$sigma
  gp <- object$jags_fit$sims.list$gp
  n_participants <- object$jags_data$Nparticipants
  if (is.null(sigma) || is.null(gp)) stop("sigma and gp draws are required for response uncertainty.")
  sigma <- matrix(sigma, ncol = 2L)
  gp <- matrix(gp, ncol = n_participants)
  if (nrow(sigma) != n_draws || nrow(gp) != n_draws) {
    stop("Observation noise and parameter draws have different sample counts.")
  }
  result <- vapply(seq_len(n_participants), function(i) {
    sigma[cbind(seq_len(n_draws), ifelse(gp[, i] == 1L, 2L, 1L))]
  }, numeric(n_draws))
  matrix(result, nrow = n_draws, ncol = n_participants)
}

with_prediction_seed <- function(seed, code) {
  if (!is.null(seed)) {
    validate_scalar(seed, "seed", lower = 0, upper = .Machine$integer.max, integer = TRUE)
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
  force(code)
}

posterior_response_draws <- function(object, predictive = FALSE,
                                     scale = c("response", "normalized"),
                                     n_draws = NULL, seed = NULL) {
  scale <- match.arg(scale)
  cached <- predictive && !is.null(object$jags_fit$sims.list$y_pre)
  node <- if (cached) "y_pre" else "theta"
  draws <- object$jags_fit$sims.list[[node]]
  expected <- c(object$jags_data$Nparticipants, object$jags_data$Ntrials)
  if (is.null(draws)) stop("Posterior ", node, " draws are unavailable; refit the model.")
  if (is.null(dim(draws)) && prod(expected) == 1L) draws <- array(draws, c(length(draws), expected))
  if (length(dim(draws)) != 3L || !identical(as.integer(dim(draws)[-1]), as.integer(expected))) {
    stop("Unexpected posterior dimensions for ", node, ".")
  }
  total <- dim(draws)[1]
  index <- seq_len(total)
  if (!is.null(n_draws)) {
    validate_scalar(n_draws, "n_draws", lower = 1, upper = total, integer = TRUE)
    index <- unique(round(seq(1, total, length.out = n_draws)))
    draws <- draws[index, , , drop = FALSE]
  }
  sigma <- NULL
  if ((predictive && !cached) || (!predictive && response_model(object) == "censored")) {
    sigma <- observation_sigma_draws(object, total)[index, , drop = FALSE]
  }
  if (predictive && !cached) {
    draws <- with_prediction_seed(seed, {
      for (i in seq_len(expected[1])) {
        values <- matrix(draws[, i, ], nrow = length(index))
        noise <- matrix(stats::rnorm(length(values)), nrow = length(index))
        draws[, i, ] <- values + noise * sigma[, i]
      }
      draws
    })
  } else if (!is.null(seed)) {
    validate_scalar(seed, "seed", lower = 0, upper = .Machine$integer.max, integer = TRUE)
  }
  if (response_model(object) == "censored") {
    if (predictive) {
      draws <- pmin(pmax(draws, 0), 10)
    } else {
      for (i in seq_len(expected[1])) {
        values <- matrix(draws[, i, ], nrow = length(index))
        draws[, i, ] <- clipped_normal_mean(values, sigma[, i])
      }
    }
  }
  map <- attr(object$jags_data, "row_index")
  if (!is.null(map)) {
    flat <- matrix(draws, nrow = dim(draws)[1])
    flat[, is.na(as.vector(map))] <- NA_real_
    dim(flat) <- dim(draws)
    draws <- flat
  }
  if (scale == "response") {
    bounds <- response_bounds(object)
    draws <- bounds[1] + draws / 10 * diff(bounds)
  }
  dimnames(draws) <- list(NULL, as.character(attr(object$jags_data, "participants")), NULL)
  draws
}

#' Posterior Predictive Responses for Fitted Trials
#'
#' Return replicated ratings from the fitted model, including observation noise.
#' These draws describe repeated observations of fitted trial histories,
#' rather than responses from new participants or future trials.
#'
#' @param object A fitted mfeargen model.
#' @param n_draws Number of stored draws to return. NULL returns all draws;
#'   otherwise draws are selected evenly across stored samples.
#' @param scale "response" for the original rating scale, or "normalized" for 0-10.
#' @param seed Optional seed for generating ratings when predictive draws were
#'   not saved during fitting. Seeded calls preserve the caller's RNG state.
#' @return An array with dimensions draw, participant, trial. Participant IDs are
#'   retained and padded trials are NA. With censored observations, replicated
#'   responses are clipped to the rating bounds.
#' @export
posterior_predict <- function(object, n_draws = NULL,
                              scale = c("response", "normalized"), seed = NULL) {
  if (!inherits(object, "mfeargen")) stop("object must be a fitted mfeargen model.")
  posterior_response_draws(object, predictive = TRUE, scale = match.arg(scale),
                           n_draws = n_draws, seed = seed)
}

trial_prediction_table <- function(object, estimate, lower = NULL, upper = NULL) {
  map <- attr(object$jags_data, "row_index")
  if (is.null(map)) stop("Trial-row mappings are unavailable; refit the model before exporting predictions.")
  cells <- which(!is.na(as.vector(map)))
  cells <- cells[order(map[cells])]
  participant_index <- ((cells - 1L) %% nrow(map)) + 1L
  trial_index <- ((cells - 1L) %/% nrow(map)) + 1L
  ids <- attr(object$jags_data, "participants")
  rows <- map[cells]
  bounds <- response_bounds(object)
  observed <- bounds[1] + object$jags_data$y[cells] / 10 * diff(bounds)
  result <- data.frame(
    participant = ids[participant_index],
    trials = if (!is.null(object$data$trials)) object$data$trials[rows] else trial_index,
    observed = observed, estimate = estimate[cells], row.names = NULL)
  if (!is.null(lower)) result$lower <- lower[cells]
  if (!is.null(upper)) result$upper <- upper[cells]
  result
}
