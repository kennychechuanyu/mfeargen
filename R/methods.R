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
    cat("Convergence:", if (isTRUE(x$convergence$converged)) "passed" else "requires review", "\n")
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
  cat("\nConvergence:", if (isTRUE(x$convergence$converged)) "passed" else "requires review", "\n")
  if (!is.null(x$dic)) cat("DIC:", round(x$dic, 2), "\n")
  invisible(x)
}

#' Extract posterior coefficient summaries
#' @param object A fitted mfeargen model.
#' @param ... Additional arguments (ignored).
#' @return A data frame with parameter names, posterior means, standard
#'   deviations, and 95 percent credible intervals. Standard deviations describe
#'   posterior uncertainty, rather than frequentist standard errors.
#' @export
coef.mfeargen <- function(object, ...) {
  s <- object$jags_fit$summary
  keep <- grepl("^(alpha|lambda|w0|w1|sigma|pi)(\\[|$)", rownames(s))
  s <- s[keep, , drop = FALSE]
  data.frame(parameter = rownames(s), estimate = s[, "mean"],
             posterior_sd = s[, "sd"], ci_lower = s[, "2.5%"],
             ci_upper = s[, "97.5%"], row.names = NULL)
}

#' Predictions for fitted participants
#' @param object A fitted mfeargen model.
#' @param newdata Reserved for future use; must be NULL.
#' @param type Either "groups" or "response".
#' @param ... Arguments passed to predict_groups or fitted.
#' @return Group assignments or posterior mean expected responses.
#' @export
predict.mfeargen <- function(object, newdata = NULL,
                            type = c("groups", "response"), ...) {
  if (!is.null(newdata)) {
    stop("Prediction for newdata is not implemented; omit newdata to predict fitted participants.")
  }
  type <- match.arg(type)
  if (type == "groups") predict_groups(object, ...) else stats::fitted(object, ...)
}

#' Posterior mean expected responses
#' @param object A fitted mfeargen model.
#' @param scale "response" for the original rating scale, or "normalized" for 0-10.
#' @param ... Additional arguments (ignored).
#' @return A participant by trial matrix. Rows retain participant IDs; trailing
#'   padded trials are NA. Values are expectations from theta, not noisy
#'   posterior predictive draws.
#' @export
fitted.mfeargen <- function(object, scale = c("response", "normalized"), ...) {
  scale <- match.arg(scale)
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
    bounds <- object$rating_scale
    if (is.null(bounds)) bounds <- attr(object$jags_data, "rating_scale")
    if (!is.null(bounds)) theta <- theta / 10 * diff(bounds) + bounds[1]
  }
  theta
}

#' Model residuals
#' @param object A fitted mfeargen model.
#' @param type "raw" on the original rating scale, or "standardized" by the
#'   participant's posterior mean observation standard deviation.
#' @param ... Additional arguments (ignored).
#' @return A participant by trial matrix, with NA for missing observations.
#' @export
residuals.mfeargen <- function(object, type = c("raw", "standardized"), ...) {
  type <- match.arg(type)
  result <- object$jags_data$y - stats::fitted(object, scale = "normalized")
  if (type == "standardized") {
    sigma <- object$jags_fit$sims.list$sigma
    gp <- object$jags_fit$sims.list$gp
    if (is.null(sigma) || is.null(gp)) stop("sigma and gp samples are required for standardized residuals.")
    sigma <- matrix(sigma, ncol = 2)
    gp <- matrix(gp, ncol = nrow(result))
    participant_sigma <- vapply(seq_len(nrow(result)), function(i) {
      mean(sigma[cbind(seq_len(nrow(gp)), ifelse(gp[, i] == 1, 2L, 1L))])
    }, numeric(1))
    result <- result / participant_sigma
  } else {
    bounds <- object$rating_scale
    if (is.null(bounds)) bounds <- attr(object$jags_data, "rating_scale")
    if (!is.null(bounds)) result <- result * diff(bounds) / 10
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
