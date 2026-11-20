# Source this file, then call analyze_trials() with your CSV path and rating bounds.
analyze_trials <- function(file, rating_scale, output_dir = "mfeargen-results",
                           model_type = "1cs", observation_model = "censored",
                           n.iter = 10000, n.burnin = 5000, n.chains = 3,
                           n.thin = 2, seed = 42, parallel = FALSE) {
  columns <- names(utils::read.csv(file, nrows = 0, check.names = FALSE))
  if (!"participant" %in% columns) stop("CSV is missing the participant column")
  trials <- utils::read.csv(file, colClasses = c(participant = "character"),
                            stringsAsFactors = FALSE, check.names = FALSE)
  mfeargen::validate_mfeargen_data(trials, model_type = model_type)
  fit <- mfeargen::fit_mfeargen(
    trials, model_type = model_type, rating_scale = rating_scale,
    observation_model = observation_model, n.iter = n.iter, n.burnin = n.burnin,
    n.chains = n.chains, n.thin = n.thin, seed = seed, parallel = parallel
  )
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
    stop("Could not create output directory: ", output_dir)
  }
  write_table <- function(data, name) {
    utils::write.csv(data, file.path(output_dir, name), row.names = FALSE)
  }
  write_table(stats::coef(fit), "parameters.csv")
  write_table(fit$group_membership, "groups.csv")
  write_table(stats::predict(fit, type = "response", interval = "credible", format = "data.frame"),
              "expected-ratings.csv")
  write_table(stats::predict(fit, type = "response", interval = "prediction", format = "data.frame", seed = seed),
              "predictive-intervals.csv")
  write_table(fit$convergence$diagnostics, "diagnostics.csv")
  saveRDS(fit, file.path(output_dir, "fit.rds"))
  ggplot2::ggsave(file.path(output_dir, "group-probabilities.png"),
    mfeargen::plot_group_probabilities(fit), width = 9, height = 5, dpi = 150)
  ggplot2::ggsave(file.path(output_dir, "posterior-predictive.png"),
    mfeargen::plot_posterior_predictive(fit, seed = seed), width = 9, height = 5, dpi = 150)
  if (!isTRUE(fit$convergence$converged)) {
    message("Sampling diagnostics require review. Inspect diagnostics.csv and trace plots before interpreting estimates.")
  }
  invisible(fit)
}
