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
  verbose = TRUE
) {
  
  if (verbose) message("Starting parameter recovery analysis...")
  
  # Generate simulated data
  if (verbose) message("1. Generating simulated data...")
  sim_data <- simulate_mfeargen(
    n_participants = n_participants,
    model_type = model_type,
    return_true_params = TRUE,
    seed = seed
  )
  
  true_params <- attr(sim_data, "true_parameters")
  
  # Fit model
  if (verbose) message("2. Fitting model to recover parameters...")
  fit <- fit_mfeargen(
    data = sim_data,
    model_type = model_type,
    n.iter = n.iter,
    n.burnin = n.burnin,
    n.chains = n.chains,
    verbose = FALSE
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
  n_participants <- nrow(true_params)
  mcmc_output <- fit$jags_fit
  
  recovered <- data.frame(
    participant = 1:n_participants,
    alpha_recovered = if("alpha" %in% names(mcmc_output$mean)) mcmc_output$mean$alpha else NA,
    lambda_recovered = if("lambda" %in% names(mcmc_output$mean)) mcmc_output$mean$lambda else NA,
    w0_recovered = if("w0" %in% names(mcmc_output$mean)) mcmc_output$mean$w0 else NA,
    w1_recovered = if("w1" %in% names(mcmc_output$mean)) mcmc_output$mean$w1 else NA,
    group_recovered = if("gp" %in% names(mcmc_output$mean)) round(mcmc_output$mean$gp) else NA
  )
  
  return(recovered)
}

#' Calculate Recovery Metrics
#' @keywords internal
calculate_recovery_metrics <- function(true_params, recovered_params) {
  combined <- cbind(true_params, recovered_params)
  metrics <- list()
  
  # Alpha correlation
  valid <- !is.na(combined$alpha_recovered) & combined$alpha_true > 0
  if (sum(valid) > 2) {
    metrics$alpha_correlation <- cor(combined$alpha_true[valid], combined$alpha_recovered[valid])
    metrics$alpha_rmse <- sqrt(mean((combined$alpha_true[valid] - combined$alpha_recovered[valid])^2))
  }
  
  # Lambda correlation
  valid <- !is.na(combined$lambda_recovered) & combined$lambda_true > 0
  if (sum(valid) > 2) {
    metrics$lambda_correlation <- cor(combined$lambda_true[valid], combined$lambda_recovered[valid])
    metrics$lambda_rmse <- sqrt(mean((combined$lambda_true[valid] - combined$lambda_recovered[valid])^2))
  }
  
  # Group accuracy
  if (!is.na(combined$group_recovered[1])) {
    confusion_matrix <- table(True = combined$group, Recovered = combined$group_recovered)
    metrics$group_accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
    metrics$confusion_matrix <- confusion_matrix
  }
  
  metrics$combined_data <- combined
  return(metrics)
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