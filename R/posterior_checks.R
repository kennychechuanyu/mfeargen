#' Compare Replicated Ratings with Observations
#'
#' Summarize posterior predictive checks of response location, spread, and the
#' proportion of ratings at or beyond each scale boundary. Only observed ratings enter
#' each statistic; padded trials and missing ratings are excluded.
#'
#' @param fit A fitted mfeargen model.
#' @param n_draws Number of predictive draws to use. NULL uses up to 100 draws.
#' @param level Width of the equal-tailed interval for replicated statistics.
#' @param seed Optional seed for generating ratings; seeded calls preserve RNG state.
#' @return A data frame containing observed and replicated statistics, intervals,
#'   and posterior predictive tail probabilities. These probabilities describe
#'   the fraction of replicated statistics at least as large as the observed
#'   statistic; they are not frequentist p-values. Spread requires two ratings.
#' @export
posterior_predictive_check <- function(fit, n_draws = NULL, level = 0.95, seed = NULL) {
  validate_scalar(level, "level", lower = 0, upper = 1, lower_open = TRUE, upper_open = TRUE)
  if (is.null(n_draws)) {
    theta <- fit$jags_fit$sims.list$theta
    available <- if (is.null(dim(theta))) length(theta) else dim(theta)[1]
    n_draws <- min(100L, available)
  }
  draws <- posterior_predict(fit, n_draws = n_draws, scale = "normalized", seed = seed)
  observed <- as.vector(fit$jags_data$y)
  keep <- !is.na(observed)
  replicated <- matrix(draws, nrow = dim(draws)[1])[, keep, drop = FALSE]
  observed <- observed[keep]
  if (!length(observed) || any(!is.finite(replicated))) {
    stop("Posterior predictive checks require finite draws for observed ratings.")
  }
  statistics <- list(mean = mean, sd = stats::sd,
                     lower_boundary = function(x) mean(x <= 0),
                     upper_boundary = function(x) mean(x >= 10))
  tail <- (1 - level) / 2
  bounds <- response_bounds(fit)
  rows <- lapply(names(statistics), function(name) {
    statistic <- statistics[[name]]
    value <- statistic(observed)
    samples <- apply(replicated, 1, statistic)
    finite <- samples[is.finite(samples)]
    multiplier <- if (name %in% c("mean", "sd")) diff(bounds) / 10 else 1
    shift <- if (name == "mean") bounds[1] else 0
    interval <- if (length(finite)) stats::quantile(finite, c(tail, 1 - tail)) else c(NA_real_, NA_real_)
    data.frame(statistic = name, observed = shift + multiplier * value,
               estimate = if (length(finite)) shift + multiplier * mean(finite) else NA_real_,
               lower = shift + multiplier * interval[1], upper = shift + multiplier * interval[2],
               tail_probability = if (is.finite(value) && length(finite)) mean(finite >= value) else NA_real_,
               row.names = NULL)
  })
  result <- do.call(rbind, rows)
  attr(result, "level") <- level
  attr(result, "n_observations") <- length(observed)
  attr(result, "n_draws") <- nrow(replicated)
  result
}

#' Plot Posterior Predictive Checks
#' @param fit A fitted mfeargen model.
#' @param n_draws Number of predictive draws to use, or NULL for up to 100.
#' @param level Interval width for replicated statistics.
#' @param seed Optional seed for generated predictive ratings.
#' @return A ggplot object. Dots and intervals describe replicated statistics;
#'   crosses show the corresponding observed statistics.
#' @export
plot_posterior_predictive <- function(fit, n_draws = NULL, level = 0.95, seed = NULL) {
  checks <- posterior_predictive_check(fit, n_draws, level, seed)
  ggplot(checks, aes(x = statistic, y = estimate)) +
    geom_linerange(aes(ymin = lower, ymax = upper), linewidth = 0.8) +
    geom_point(size = 2.5) +
    geom_point(aes(y = observed), shape = 4, size = 3, color = "#D55E00") +
    facet_wrap(~ statistic, scales = "free") +
    labs(title = "Posterior predictive checks", x = NULL, y = "Statistic",
         subtitle = "Intervals: replicated ratings; crosses: observed ratings") +
    theme_minimal() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}
