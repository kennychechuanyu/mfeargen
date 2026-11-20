#' @title Parameter Recovery Functions for mfeargen
#' @description Functions for parameter recovery analysis
#' @name parameter_recovery
NULL

#' Run Parameter Recovery Analysis
#'
#' @description
#' Simulates data with known parameters and attempts to recover them
#' using the mfeargen model fitting procedure.
#'
#' @param n_participants Integer. Number of participants to simulate
#' @param model_type Character. "1cs" or "2cs"
#' @param n.iter Integer. MCMC iterations (default: 10000)
#' @param n.burnin Integer. MCMC burnin (default: 5000)
#' @param n.chains Integer. Number of chains (default: 3)
#' @param seed Integer. Random seed
#' @param verbose Logical. Print progress messages
#' @param simulation_args Named list of additional arguments to simulate_mfeargen.
#' @param ... Additional arguments to fit_mfeargen, such as parallel = FALSE.
#'
#' @return mfeargen_recovery object containing recovery results
#' @export
#' @examples
#' \dontrun{
#' recovery <- run_parameter_recovery(n_participants = 20)
#' plot_parameter_recovery(recovery)
#' }
run_parameter_recovery <- function(
  n_participants = 40,
  model_type = "1cs",
  n.iter = 10000,
  n.burnin = 5000,
  n.chains = 3,
  seed = NULL,
  verbose = TRUE,
  simulation_args = list(),
  ...
) {
  
  if (verbose) message("Starting parameter recovery analysis...")
  
  # Generate simulated data
  if (verbose) message("1. Generating simulated data...")
  if (!is.list(simulation_args) || (length(simulation_args) &&
      (is.null(names(simulation_args)) || any(!nzchar(names(simulation_args)))))) {
    stop("simulation_args must be a named list")
  }
  reserved <- c("n_participants", "model_type", "return_true_params", "seed")
  if (any(names(simulation_args) %in% reserved)) {
    stop("Set n_participants, model_type and seed directly, not in simulation_args")
  }
  sim_data <- do.call(simulate_mfeargen, c(list(
    n_participants = n_participants, model_type = model_type,
    return_true_params = TRUE, seed = seed), simulation_args))

  true_params <- attr(sim_data, "true_parameters")
  
  # Fit model
  if (verbose) message("2. Fitting model to recover parameters...")
  fit <- fit_mfeargen(
    data = sim_data,
    model_type = model_type,
    n.iter = n.iter,
    n.burnin = n.burnin,
    n.chains = n.chains,
    verbose = FALSE,
    seed = seed,
    ...
  )
  
  # Extract recovered parameters
  if (verbose) message("3. Extracting recovered parameters...")
  recovered_params <- extract_recovered_parameters(fit, true_params)
  
  # Calculate metrics
  if (verbose) message("4. Calculating recovery metrics...")
  metrics <- calculate_recovery_metrics(true_params, recovered_params)
  
  # Prepare results
  results <- list(
    true_parameters = true_params,
    recovered_parameters = recovered_params,
    metrics = metrics,
    fit = fit,
    data = sim_data,
    settings = list(
      n_participants = n_participants,
      model_type = model_type,
      n.iter = n.iter,
      n.burnin = n.burnin,
      n.chains = n.chains
    )
  )
  
  class(results) <- c("mfeargen_recovery", "list")
  
  if (verbose) {
    message("\nParameter recovery complete!")
    print_recovery_summary(metrics)
  }
  
  return(results)
}

#' Extract Recovered Parameters
#' @keywords internal
extract_recovered_parameters <- function(fit, true_params) {
  estimates <- fit$jags_fit$mean
  ids <- attr(fit$jags_data, "participants")
  if (is.null(ids)) ids <- true_params$participant
  recovered <- data.frame(participant = ids)
  for (parameter in c("alpha", "lambda", "w0", "w1")) {
    values <- estimates[[parameter]]
    if (is.null(values)) stop("Missing posterior parameter: ", parameter)
    recovered[[paste0(parameter, "_recovered")]] <- as.vector(values)
  }
  membership <- predict_groups(fit, threshold = 0)
  recovered$group_recovered <- membership$group[match(ids, membership$participant)]
  recovered[match(true_params$participant, recovered$participant), , drop = FALSE]
}

#' Calculate Recovery Metrics
#' @keywords internal
calculate_recovery_metrics <- function(true_params, recovered_params) {
  index <- match(true_params$participant, recovered_params$participant)
  if (anyNA(index)) stop("Recovered participant IDs do not match simulated IDs")
  recovered <- recovered_params[index, , drop = FALSE]
  combined <- cbind(true_params, recovered[setdiff(names(recovered), "participant")])
  metrics <- list()
  for (parameter in c("alpha", "lambda", "w0", "w1")) {
    truth <- combined[[paste0(parameter, "_true")]]
    estimate <- combined[[paste0(parameter, "_recovered")]]
    valid <- is.finite(truth) & is.finite(estimate)
    if (parameter == "alpha") valid <- valid & combined$group != 1L
    if (parameter == "lambda") valid <- valid & combined$group > 2L
    metrics[[paste0(parameter, "_correlation")]] <- if (sum(valid) > 2L &&
        stats::sd(truth[valid]) > 0 && stats::sd(estimate[valid]) > 0) {
      stats::cor(truth[valid], estimate[valid])
    } else NA_real_
    metrics[[paste0(parameter, "_rmse")]] <- if (any(valid)) {
      sqrt(mean((truth[valid] - estimate[valid])^2))
    } else NA_real_
  }
  valid <- !is.na(combined$group_recovered)
  metrics$group_accuracy <- if (any(valid)) {
    mean(combined$group[valid] == combined$group_recovered[valid])
  } else NA_real_
  metrics$confusion_matrix <- table(
    True = factor(combined$group, levels = 1:4),
    Recovered = factor(combined$group_recovered, levels = 1:4))
  metrics$combined_data <- combined
  metrics
}

#' Print Recovery Summary
#' @keywords internal
print_recovery_summary <- function(metrics) {
  cat("\nParameter Recovery Results:\n")
  cat("---------------------------\n")
  
  if (!is.null(metrics$alpha_correlation)) {
    cat(sprintf("Alpha:  r = %.3f, RMSE = %.3f\n", 
                metrics$alpha_correlation, metrics$alpha_rmse))
  }
  
  if (!is.null(metrics$lambda_correlation)) {
    cat(sprintf("Lambda: r = %.3f, RMSE = %.3f\n", 
                metrics$lambda_correlation, metrics$lambda_rmse))
  }
  
  if (!is.null(metrics$group_accuracy)) {
    cat(sprintf("Group Accuracy: %.1f%%\n", metrics$group_accuracy * 100))
  }
}