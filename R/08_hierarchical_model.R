#' Repeated-Session Model Structure
#'
#' The repeated-session model separates population locations, stable participant
#' components, and session-specific parameters. Learning rates are bounded by
#' zero and one, response slopes are positive, and group membership constrains
#' generalization gradients. Sessions share a participant's latent group but
#' have separate conditioning histories.
#'
#' Variance parameters describe the underlying normal distributions before
#' truncation and group constraints. The ratio
#' \eqn{\tau^2 / (\tau^2 + \sigma^2)} summarizes between-person variance
#' relative to between-person plus session variance on that underlying scale.
#' It is not the empirical intraclass correlation of constrained parameters.
#'
#' @name hierarchical_structure
NULL

utils::globalVariables(c(
  "individual", "stable_value", "session_value", "session", "lower", "upper",
  "Stability", "Parameter", "CI_Lower", "CI_Upper", "Interpretation"
))

hierarchical_draws <- function(fit, parameter, suffix = "") {
  if (!inherits(fit, "mfeargen_hierarchical")) {
    stop("fit must be an mfeargen_hierarchical object")
  }
  parameter <- match.arg(parameter, c("alpha", "lambda", "w0", "w1"))
  name <- paste0(parameter, suffix)
  draws <- fit$jags_fit$sims.list[[name]]
  if (is.null(draws)) stop("Posterior draws for ", name, " are unavailable")
  draws
}

hierarchical_variance_draws <- function(fit, parameter) {
  if (!inherits(fit, "mfeargen_hierarchical")) {
    stop("fit must be an mfeargen_hierarchical object")
  }
  parameter <- match.arg(parameter, c("alpha", "lambda", "w0", "w1"))
  draws <- fit$jags_fit$sims.list
  between <- draws[[paste0("tau_", parameter, "_sq")]]
  within <- draws[[paste0("sigma_", parameter, "_sq")]]
  if (is.null(between) || is.null(within) || length(between) != length(within) ||
      !length(between) || any(!is.finite(between)) || any(!is.finite(within)) ||
      any(between <= 0) || any(within <= 0)) {
    stop("Positive posterior variance draws are unavailable for ", parameter)
  }
  list(between = as.numeric(between), within = as.numeric(within))
}

#' Summarize Between-Person and Session Variance
#'
#' @param fit An `mfeargen_hierarchical` object.
#' @param parameters Parameters to summarize.
#' @return An `mfeargen_variance_decomposition` list containing posterior means
#'   of between-person, within-person, and total variance, and posterior means
#'   of their draw-wise variance ratios. See [hierarchical_structure] for scale
#'   and interpretation.
#' @export
calculate_variance_decomposition <- function(fit,
                                            parameters = c("alpha", "lambda", "w0", "w1")) {
  if (!is.character(parameters) || !length(parameters) || anyDuplicated(parameters)) {
    stop("parameters must be a nonempty vector of distinct parameter names")
  }
  results <- list(between_variance = list(), within_variance = list(),
                  total_variance = list(), stability = list(), icc = list())
  for (parameter in parameters) {
    draws <- hierarchical_variance_draws(fit, parameter)
    results$between_variance[[parameter]] <- mean(draws$between)
    results$within_variance[[parameter]] <- mean(draws$within)
    results$total_variance[[parameter]] <- mean(draws$between + draws$within)
    results$stability[[parameter]] <- mean(draws$between / (draws$between + draws$within))
    results$icc[[parameter]] <- results$stability[[parameter]]
  }
  results$interpretation <- lapply(results$stability, interpret_stability_value)
  class(results) <- c("mfeargen_variance_decomposition", "list")
  results
}

#' Summarize a Parameter's Variance Ratio
#'
#' @param fit An `mfeargen_hierarchical` object.
#' @param parameter One of `"alpha"`, `"lambda"`, `"w0"`, or `"w1"`.
#' @param credible_interval Credible interval width between zero and one.
#' @return An `mfeargen_stability` list with posterior summaries and draws.
#'   The ratio refers to the underlying normal distributions, before parameter
#'   constraints. The labels Low, Moderate, and High use descriptive cutoffs of
#'   0.3 and 0.7; they are not evidence of a psychological trait.
#' @export
calculate_stability <- function(fit, parameter, credible_interval = 0.95) {
  validate_scalar(credible_interval, "credible_interval", lower = 0, upper = 1,
                  lower_open = TRUE, upper_open = TRUE)
  draws <- hierarchical_variance_draws(fit, parameter)
  ratio <- draws$between / (draws$between + draws$within)
  tail <- (1 - credible_interval) / 2
  results <- list(
    parameter = parameter,
    mean_stability = mean(ratio),
    median_stability = stats::median(ratio),
    sd_stability = stats::sd(ratio),
    ci_lower = unname(stats::quantile(ratio, tail)),
    ci_upper = unname(stats::quantile(ratio, 1 - tail)),
    credible_interval = credible_interval,
    between_variance = mean(draws$between),
    within_variance = mean(draws$within),
    variance_ratio = mean(draws$between) / mean(draws$within),
    samples = ratio,
    interpretation = interpret_stability_value(mean(ratio))
  )
  class(results) <- c("mfeargen_stability", "list")
  results
}

#' Interpret a Variance Ratio
#' @keywords internal
interpret_stability_value <- function(stability) {
  validate_scalar(stability, "stability", lower = 0, upper = 1)
  if (stability >= 0.7) {
    list(level = "High", description = "Most underlying variance is between people")
  } else if (stability >= 0.3) {
    list(level = "Moderate", description = "Both between-person and session variance contribute")
  } else {
    list(level = "Low", description = "Most underlying variance is between sessions")
  }
}

#' Interpret Variance Ratios
#' @keywords internal
interpret_stability <- function(stability_list) {
  lapply(stability_list, interpret_stability_value)
}

#' Extract Stable Participant Parameters
#'
#' @param fit An `mfeargen_hierarchical` object.
#' @param parameter Parameter to extract.
#' @return A matrix with participants in rows and posterior draws in columns.
#'   Row names retain participant IDs.
#' @export
extract_stable_components <- function(fit, parameter = "alpha") {
  draws <- hierarchical_draws(fit, parameter, "_stable")
  n_participants <- fit$data_info$n_participants
  draws <- matrix(draws, ncol = n_participants)
  result <- t(draws)
  rownames(result) <- as.character(fit$data_info$participants)
  result
}

#' Extract Session Deviations from Stable Parameters
#'
#' @param fit An `mfeargen_hierarchical` object.
#' @param parameter Parameter to extract.
#' @return An array with dimensions participant, session, posterior draw.
#'   Missing participant-session combinations contain `NA`. The `sessions`
#'   attribute records each participant's original session IDs.
#' @export
extract_session_deviations <- function(fit, parameter = "alpha") {
  stable <- extract_stable_components(fit, parameter)
  session_draws <- hierarchical_draws(fit, parameter)
  n_participants <- fit$data_info$n_participants
  n_sessions <- max(fit$jags_data$Nsessions)
  session_draws <- array(session_draws, c(ncol(stable), n_participants, n_sessions))
  result <- array(NA_real_, c(n_participants, n_sessions, ncol(stable)),
                  dimnames = list(rownames(stable), as.character(seq_len(n_sessions)), NULL))
  for (i in seq_len(n_participants)) {
    for (s in seq_len(fit$jags_data$Nsessions[i])) {
      result[i, s, ] <- session_draws[, i, s] - stable[i, ]
    }
  }
  attr(result, "sessions") <- fit$data_info$sessions
  result
}

#' Plot Stable and Session Parameters
#'
#' Plots posterior means and 95% intervals for actual participants and sessions.
#'
#' @param fit An `mfeargen_hierarchical` object.
#' @param parameter Parameter to visualize.
#' @param n_individuals Maximum number of participants to display.
#' @return A ggplot object.
#' @export
plot_hierarchical_structure <- function(fit, parameter = "alpha", n_individuals = 10) {
  validate_scalar(n_individuals, "n_individuals", lower = 1, integer = TRUE)
  stable <- extract_stable_components(fit, parameter)
  deviations <- extract_session_deviations(fit, parameter)
  n_show <- min(n_individuals, nrow(stable))
  summary_row <- function(draws) {
    c(mean = mean(draws), lower = unname(stats::quantile(draws, 0.025)),
      upper = unname(stats::quantile(draws, 0.975)))
  }
  stable_summary <- t(vapply(seq_len(n_show), function(i) summary_row(stable[i, ]), numeric(3)))
  stable_data <- data.frame(
    individual = rownames(stable)[seq_len(n_show)], stable_value = stable_summary[, "mean"],
    lower = stable_summary[, "lower"], upper = stable_summary[, "upper"]
  )
  rows <- list()
  for (i in seq_len(n_show)) {
    for (s in seq_len(fit$jags_data$Nsessions[i])) {
      values <- summary_row(deviations[i, s, ] + stable[i, ])
      rows[[length(rows) + 1L]] <- data.frame(
        individual = rownames(stable)[i], session = as.character(fit$data_info$sessions[[i]][s]),
        stable_value = stable_summary[i, "mean"], session_value = values["mean"],
        lower = values["lower"], upper = values["upper"]
      )
    }
  }
  session_data <- do.call(rbind, rows)
  population <- fit$population_parameters[[paste0("mu_", parameter)]]
  ggplot2::ggplot(session_data, ggplot2::aes(x = individual, y = session_value)) +
    ggplot2::geom_hline(yintercept = population, linetype = "dashed", color = "grey50") +
    ggplot2::geom_segment(ggplot2::aes(xend = individual, y = stable_value, yend = session_value),
                          alpha = 0.3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.1, alpha = 0.5) +
    ggplot2::geom_point(ggplot2::aes(color = session), size = 2) +
    ggplot2::geom_errorbar(data = stable_data,
      ggplot2::aes(x = individual, ymin = lower, ymax = upper),
      inherit.aes = FALSE, width = 0.3, color = "black") +
    ggplot2::geom_point(data = stable_data,
      ggplot2::aes(x = individual, y = stable_value), inherit.aes = FALSE,
      shape = 18, size = 3, color = "black") +
    ggplot2::labs(title = paste("Repeated-session parameters:", parameter),
      subtitle = "Posterior means and 95% intervals; black diamonds show stable components",
      x = "Participant", y = parameter, color = "Session") +
    ggplot2::theme_minimal()
}

#' Compare Variance Ratios Across Parameters
#' @param fit An `mfeargen_hierarchical` object.
#' @param parameters Parameters to compare.
#' @return A data frame with posterior means, 95% intervals, variance components,
#'   and descriptive labels.
#' @export
compare_parameter_stability <- function(fit, parameters = c("alpha", "lambda", "w0", "w1")) {
  if (!is.character(parameters) || !length(parameters) || anyDuplicated(parameters)) {
    stop("parameters must be a nonempty vector of distinct parameter names")
  }
  rows <- lapply(parameters, function(parameter) {
    result <- calculate_stability(fit, parameter)
    data.frame(
      Parameter = parameter, Stability = result$mean_stability,
      CI_Lower = result$ci_lower, CI_Upper = result$ci_upper,
      Between_Var = result$between_variance, Within_Var = result$within_variance,
      Interpretation = result$interpretation$level
    )
  })
  comparison <- do.call(rbind, rows)
  comparison[order(comparison$Stability, decreasing = TRUE), ]
}

#' Plot Variance Ratio Comparisons
#' @param stability_comparison Output from [compare_parameter_stability()].
#' @return A ggplot object.
#' @export
plot_stability_comparison <- function(stability_comparison) {
  required <- c("Parameter", "Stability", "CI_Lower", "CI_Upper", "Interpretation")
  if (!is.data.frame(stability_comparison) || !all(required %in% names(stability_comparison)) ||
      !nrow(stability_comparison)) stop("stability_comparison must contain variance ratio summaries")
  stability_comparison$Parameter <- factor(stability_comparison$Parameter,
    levels = stability_comparison$Parameter[order(stability_comparison$Stability)])
  ggplot2::ggplot(stability_comparison, ggplot2::aes(x = Parameter, y = Stability)) +
    ggplot2::geom_hline(yintercept = c(0.3, 0.7), linetype = "dashed", alpha = 0.3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
    ggplot2::geom_point(ggplot2::aes(color = Interpretation), size = 3) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Parameter", y = "Between-person / total underlying variance",
                  color = "Ratio") +
    ggplot2::theme_minimal()
}

#' Print Variance Components
#' @param x An `mfeargen_variance_decomposition` object.
#' @param ... Unused.
#' @export
print.mfeargen_variance_decomposition <- function(x, ...) {
  cat("Underlying parameter variance components\n")
  for (parameter in names(x$stability)) {
    cat(sprintf("%s: between %.4f; session %.4f; mean variance ratio %.3f\n",
      parameter, x$between_variance[[parameter]], x$within_variance[[parameter]],
      x$stability[[parameter]]))
  }
  invisible(x)
}

#' Print a Variance Ratio Summary
#' @param x An `mfeargen_stability` object.
#' @param ... Unused.
#' @export
print.mfeargen_stability <- function(x, ...) {
  cat(sprintf("%s: mean underlying variance ratio %.3f\n", x$parameter, x$mean_stability))
  cat(sprintf("%.0f%% interval: [%.3f, %.3f]\n", 100 * x$credible_interval,
                x$ci_lower, x$ci_upper))
  invisible(x)
}
