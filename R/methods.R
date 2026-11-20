#' Print a fitted fear generalization model
#' @param x A fitted mfeargen model.
#' @param ... Additional arguments (ignored).
#' @return The input, invisibly.
#' @export
print.mfeargen <- function(x, ...) {
  cat("Fear generalization mixture model (", x$model_type, ")\n", sep = "")
  cat(x$jags_data$Nparticipants, "participants; up to",
      x$jags_data$Ntrials, "trials per participant\n")
  cat(x$n.chains, "chains;", x$n.iter, "iterations;", x$n.burnin, "burn-in\n")
  if (!is.null(x$convergence)) {
    print_convergence_status(x$convergence)
  }
  if (!is.null(x$group_membership)) print(table(x$group_membership$group_label))
  invisible(x)
}

#' Summarize a fitted model
#' @param object A fitted mfeargen model.
#' @param ... Additional arguments (ignored).
#' @return Parameter summaries, group probabilities, and convergence diagnostics.
#' @export
summary.mfeargen <- function(object, ...) {
  s <- object$jags_fit$summary
  structure(list(
    call = object$call,
    model_type = object$model_type,
    n_participants = object$jags_data$Nparticipants,
    n_trials = object$jags_data$Ntrials,
    parameters = stats::coef(object),
    group_membership = object$group_membership,
    group_probabilities = s[grepl("^pi\\[", rownames(s)), , drop = FALSE],
    convergence = object$convergence,
    dic = object$jags_fit$DIC
  ), class = "summary.mfeargen")
}

#' Print model summaries
#' @param x A summary.mfeargen object.
#' @param ... Additional arguments (ignored).
#' @return The input, invisibly.
#' @export
print.summary.mfeargen <- function(x, ...) {
  cat("Fear generalization model summary (", x$model_type, ")\n", sep = "")
  print(x$call)
  cat("\nGroup probabilities:\n")
  print(round(x$group_probabilities, 3))
  cat("\nParameter estimates:\n")
  print(x$parameters)
  cat("\n")
  print_convergence_status(x$convergence)
  if (!is.null(x$dic)) cat("DIC:", round(x$dic, 2), "\n")
  invisible(x)
}

#' Extract Posterior Coefficient Summaries
#' @param object A fitted mfeargen model.
#' @param scale "response" expresses observation noise in original rating units;
#'   "normalized" keeps it on the internal 0-10 scale. Other parameters retain
#'   their model units in both cases.
#' @param ... Additional arguments (ignored).
#' @return A data frame with parameter names, participant IDs where applicable,
#'   group components, posterior means, posterior standard deviations, and 95
#'   percent credible intervals. Posterior standard deviations describe parameter
#'   uncertainty rather than frequentist standard errors.
#' @export
coef.mfeargen <- function(object, scale = c("response", "normalized"), ...) {
  scale <- match.arg(scale)
  s <- object$jags_fit$summary
  keep <- grepl("^(alpha|lambda|w0|w1|sigma|pi)(\\[|$)", rownames(s))
  s <- s[keep, , drop = FALSE]
  result <- data.frame(parameter = rownames(s), estimate = s[, "mean"],
                       posterior_sd = s[, "sd"], ci_lower = s[, "2.5%"],
                       ci_upper = s[, "97.5%"], row.names = NULL)
  base <- sub("\\[.*", "", result$parameter)
  index <- rep(1L, nrow(result))
  indexed <- grepl("[", result$parameter, fixed = TRUE)
  index[indexed] <- as.integer(sub(".*\\[([0-9]+)\\]$", "\\1", result$parameter[indexed]))
  result$participant <- NA_character_
  ids <- attr(object$jags_data, "participants")
  individual <- base %in% c("alpha", "lambda", "w0", "w1")
  if (!is.null(ids)) result$participant[individual] <- as.character(ids[index[individual]])
  result$component <- NA_character_
  result$component[base == "pi"] <- get_group_labels(FALSE)[index[base == "pi"]]
  result$component[base == "sigma"] <- c("Learners", "Non-learners")[index[base == "sigma"]]
  if (scale == "response") {
    noise <- base == "sigma"
    multiplier <- diff(response_bounds(object)) / 10
    result[noise, c("estimate", "posterior_sd", "ci_lower", "ci_upper")] <-
      result[noise, c("estimate", "posterior_sd", "ci_lower", "ci_upper")] * multiplier
  }
  result
}

#' Predictions and Intervals for Fitted Trials
#' @param object A fitted mfeargen model.
#' @param newdata Reserved for future use; must be NULL.
#' @param type Either "groups" or "response".
#' @param interval "none", "credible" for uncertainty in the expected response,
#'   or "prediction" for repeated ratings including observation noise.
#' @param level Probability covered by the equal-tailed interval.
#' @param format "matrix" or "data.frame". Intervals in matrix format return
#'   a participant by trial by statistic array. Data frames preserve input row
#'   order and include participant IDs, trial labels, observations, and estimates.
#' @param scale "response" for original units, or "normalized" for 0-10.
#' @param seed Optional seed for generating predictive ratings.
#' @param ... Arguments passed to predict_groups for group predictions.
#' @return Group assignments, expected ratings, or predictions and intervals.
#' @export
predict.mfeargen <- function(object, newdata = NULL,
                            type = c("groups", "response"),
                            interval = c("none", "credible", "prediction"),
                            level = 0.95, format = c("matrix", "data.frame"),
                            scale = c("response", "normalized"), seed = NULL, ...) {
  if (!is.null(newdata)) {
    stop("Prediction for newdata is unavailable; omit newdata to predict fitted participants.")
  }
  type <- match.arg(type)
  interval <- match.arg(interval)
  format <- match.arg(format)
  scale <- match.arg(scale)
  if (type == "groups") {
    if (interval != "none" || format != "matrix" || scale != "response") {
      stop("Response intervals and output formats require type = 'response'.")
    }
    return(predict_groups(object, ...))
  }
  validate_scalar(level, "level", lower = 0, upper = 1, lower_open = TRUE, upper_open = TRUE)
  estimate <- stats::fitted(object, scale = scale)
  lower <- upper <- NULL
  if (interval != "none") {
    draws <- posterior_response_draws(object, predictive = interval == "prediction", scale = scale, seed = seed)
    tail <- (1 - level) / 2
    limits <- apply(draws, c(2, 3), stats::quantile, probs = c(tail, 1 - tail), na.rm = TRUE)
    lower <- matrix(limits[1, , ], nrow = nrow(estimate), dimnames = dimnames(estimate))
    upper <- matrix(limits[2, , ], nrow = nrow(estimate), dimnames = dimnames(estimate))
  }
  if (format == "data.frame") {
    result <- trial_prediction_table(object, estimate, lower, upper)
    if (scale == "normalized") {
      bounds <- response_bounds(object)
      result$observed <- 10 * (result$observed - bounds[1]) / diff(bounds)
    }
  } else if (interval == "none") {
    result <- estimate
  } else {
    result <- array(c(estimate, lower, upper), c(dim(estimate), 3L),
                    dimnames = c(dimnames(estimate), list(c("estimate", "lower", "upper"))))
  }
  attr(result, "interval") <- interval
  if (interval != "none") attr(result, "level") <- level
  result
}

#' Posterior Mean Expected Ratings
#' @param object A fitted mfeargen model.
#' @param scale "response" for the original rating scale, or "normalized" for 0-10.
#' @param ... Additional arguments (ignored).
#' @return A participant by trial matrix retaining participant IDs and NA for
#'   padded trials. For censored models, these are the expectations of bounded
#'   ratings, integrating over observation noise as well as parameter uncertainty.
#' @export
fitted.mfeargen <- function(object, scale = c("response", "normalized"), ...) {
  scale <- match.arg(scale)
  if (response_model(object) == "censored") {
    return(apply(posterior_response_draws(object, scale = scale), c(2, 3), mean))
  }
  theta <- object$jags_fit$mean$theta
  if (is.null(theta) && !is.null(object$jags_fit$sims.list$theta)) {
    theta <- apply(object$jags_fit$sims.list$theta, c(2, 3), mean)
  }
  if (is.null(theta)) stop("Expected responses were not monitored; refit with theta monitoring.")
  if (length(theta) != length(object$jags_data$y)) stop("Unexpected theta dimensions in model output.")
  theta <- matrix(theta, nrow = object$jags_data$Nparticipants,
                  ncol = object$jags_data$Ntrials)
  rownames(theta) <- as.character(attr(object$jags_data, "participants"))
  row_index <- attr(object$jags_data, "row_index")
  if (!is.null(row_index)) theta[is.na(row_index)] <- NA_real_
  if (scale == "response") {
    bounds <- response_bounds(object)
    theta <- theta / 10 * diff(bounds) + bounds[1]
  }
  theta
}

#' Model residuals
#' @param object A fitted mfeargen model.
#' @param type "raw" on the original rating scale, or "standardized" by the
#'   posterior mean conditional observation standard deviation.
#' @param ... Additional arguments (ignored).
#' @return A participant by trial matrix, with NA for missing observations.
#' @export
residuals.mfeargen <- function(object, type = c("raw", "standardized"), ...) {
  type <- match.arg(type)
  result <- object$jags_data$y - stats::fitted(object, scale = "normalized")
  if (type == "standardized") {
    n_draws <- length(object$jags_fit$sims.list$gp) / nrow(result)
    sigma <- observation_sigma_draws(object, n_draws)
    if (response_model(object) == "censored") {
      theta <- object$jags_fit$sims.list$theta
      if (is.null(theta)) stop("theta draws are required for standardized residuals.")
      conditional_sd <- array(NA_real_, dim(theta))
      for (i in seq_len(nrow(result))) {
        values <- matrix(theta[, i, ], nrow = n_draws)
        conditional_sd[, i, ] <- clipped_normal_sd(values, sigma[, i])
      }
      participant_sigma <- apply(conditional_sd, c(2, 3), mean)
    } else {
      participant_sigma <- colMeans(sigma)
    }
    result <- result / participant_sigma
  } else {
    bounds <- object$rating_scale
    if (is.null(bounds)) bounds <- attr(object$jags_data, "rating_scale")
    if (!is.null(bounds)) result <- result / 10 * diff(bounds)
  }
  result
}

#' Posterior mean conditional log likelihood
#' @param object A fitted mfeargen model.
#' @param ... Additional arguments (ignored).
#' @return A logLik object. The parameter count is unavailable for this mixture
#'   model, so df is NA; use calculate_waic for predictive comparison.
#'   The likelihood is evaluated on the internal 0-10 rating scale.
#' @export
logLik.mfeargen <- function(object, ...) {
  ll <- object$jags_fit$mean$loglik
  if (is.null(ll)) stop("Log likelihood was not monitored; refit with loglik monitoring.")
  observed <- !is.na(object$jags_data$y)
  value <- sum(ll[observed])
  structure(value, class = "logLik", df = NA_integer_, nobs = sum(observed))
}


print_convergence_status <- function(convergence) {
  status <- convergence$status %||% if (isTRUE(convergence$converged)) "passed" else "requires review"
  cat("Sampling diagnostics:", status, "\n")
  if (!is.null(convergence$max_rhat)) cat("Maximum rank Rhat:", round(convergence$max_rhat, 3), "\n")
  if (!is.null(convergence$min_neff)) cat("Minimum bulk/tail ESS:", round(convergence$min_neff), "\n")
  if (!isTRUE(convergence$converged)) cat("Review the parameter diagnostics and trace plots before interpreting estimates.\n")
}
