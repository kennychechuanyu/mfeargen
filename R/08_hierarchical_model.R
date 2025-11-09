#' @title Three-Level Hierarchical Model for mfeargen
#' @description Implementation of three-level hierarchical structure separating
#' population, individual-stable, and session-specific variation
#' @name hierarchical_model
#' @importFrom stats var sd quantile cor median rnorm
NULL

# Global variables for ggplot2 (to avoid R CMD check notes)
utils::globalVariables(c(
  "individual", "stable_value", "session_value", "session",
  "Stability", "Parameter", "CI_Lower", "CI_Upper", "Interpretation"
))

#' Three-Level Hierarchical Model Structure
#' 
#' @description
#' The mfeargen package implements a three-level hierarchical Bayesian model
#' that separates population-level distributions, individual stable components
#' (between-person differences), and session-specific realizations (within-person
#' variation).
#' 
#' @details
#' ## Model Hierarchy
#' 
#' ### Level 1: Population Distributions (Hyperparameters)
#' - \eqn{\mu_\alpha, \tau^2_\alpha}: Population mean and variance for learning rates
#' - \eqn{\mu_\lambda, \tau^2_\lambda}: Population mean and variance for generalization
#' - \eqn{\mu_{\beta_0}, \tau^2_{\beta_0}}: Population mean and variance for intercepts
#' - \eqn{\mu_{\beta_1}, \tau^2_{\beta_1}}: Population mean and variance for slopes
#' 
#' ### Level 2: Individual Stable Components (Between-Person)
#' - \eqn{\alpha^{stable}_i \sim \mathcal{N}(\mu_\alpha, \tau^2_\alpha)}
#' - \eqn{\lambda^{stable}_i \sim \mathcal{N}(\mu_\lambda, \tau^2_\lambda)}
#' - \eqn{\beta^{stable}_{0,i} \sim \mathcal{N}(\mu_{\beta_0}, \tau^2_{\beta_0})}
#' - \eqn{\beta^{stable}_{1,i} \sim \mathcal{N}(\mu_{\beta_1}, \tau^2_{\beta_1})}
#' 
#' ### Level 3: Session-Specific Realizations (Within-Person)
#' - \eqn{\alpha_{i,t} \sim \mathcal{N}(\alpha^{stable}_i, \sigma^2_\alpha)}
#' - \eqn{\lambda_{i,t} \sim \mathcal{N}(\lambda^{stable}_i, \sigma^2_\lambda)}
#' - \eqn{\beta_{0,i,t} \sim \mathcal{N}(\beta^{stable}_{0,i}, \sigma^2_{\beta_0})}
#' - \eqn{\beta_{1,i,t} \sim \mathcal{N}(\beta^{stable}_{1,i}, \sigma^2_{\beta_1})}
#' 
#' ## Variance Components
#' 
#' - **Between-person variance** (\eqn{\tau^2}): How much individuals differ in stable traits
#' - **Within-person variance** (\eqn{\sigma^2}): How much individuals vary across sessions
#' - **Individual means** (\eqn{\theta^{stable}_i}): Person i's characteristic value
#' - **Session deviations** (\eqn{\delta_{i,t} = \theta_{i,t} - \theta^{stable}_i})
#' 
#' ## Stability Metrics
#' 
#' Stability = \eqn{\tau^2 / (\tau^2 + \sigma^2)}
#' 
#' - High stability (>0.7): Most variance is between people
#' - Moderate stability (0.3-0.7): Mixed between/within variance
#' - Low stability (<0.3): Most variance is within people
#' 
#' @name hierarchical_structure
NULL

#' Calculate Variance Decomposition
#'
#' @description
#' Decomposes total variance into between-person (stable) and within-person
#' (session-specific) components for all model parameters.
#'
#' @param fit mfeargen_hierarchical object from fit_mfeargen_hierarchical()
#' @param parameters Character vector of parameters to decompose
#' @return List containing variance components and stability metrics
#' @export
#' @examples
#' \dontrun{
#' # Fit hierarchical model
#' fit <- fit_mfeargen_hierarchical(data)
#' 
#' # Calculate variance decomposition
#' var_decomp <- calculate_variance_decomposition(fit)
#' print(var_decomp$stability)
#' }
calculate_variance_decomposition <- function(fit, 
                                            parameters = c("alpha", "lambda", "w0", "w1")) {
  
  if (!inherits(fit, "mfeargen_hierarchical")) {
    stop("Input must be a mfeargen_hierarchical object")
  }
  
  # Extract MCMC samples
  mcmc_samples <- fit$jags_fit$samples
  
  # Initialize results
  results <- list(
    between_variance = list(),
    within_variance = list(),
    total_variance = list(),
    stability = list(),
    icc = list()
  )
  
  for (param in parameters) {
    # Extract variance components from MCMC
    tau_sq_name <- paste0("tau_", param, "_sq")
    sigma_sq_name <- paste0("sigma_", param, "_sq")
    
    if (tau_sq_name %in% names(mcmc_samples[[1]]) && 
        sigma_sq_name %in% names(mcmc_samples[[1]])) {
      
      # Get posterior means
      tau_sq <- mean(unlist(lapply(mcmc_samples, function(x) x[[tau_sq_name]])))
      sigma_sq <- mean(unlist(lapply(mcmc_samples, function(x) x[[sigma_sq_name]])))
      
      # Calculate variance components
      results$between_variance[[param]] <- tau_sq
      results$within_variance[[param]] <- sigma_sq
      results$total_variance[[param]] <- tau_sq + sigma_sq
      
      # Calculate stability (ICC)
      results$stability[[param]] <- tau_sq / (tau_sq + sigma_sq)
      results$icc[[param]] <- results$stability[[param]]
    }
  }
  
  # Add interpretation
  results$interpretation <- interpret_stability(results$stability)
  
  class(results) <- c("mfeargen_variance_decomposition", "list")
  return(results)
}

#' Calculate Parameter Stability
#'
#' @description
#' Calculates stability metrics for model parameters, quantifying the proportion
#' of variance attributable to stable individual differences vs. session variation.
#'
#' @param fit mfeargen_hierarchical object
#' @param parameter Character. Parameter name ("alpha", "lambda", "w0", "w1")
#' @param credible_interval Numeric. Credible interval width (default: 0.95)
#' @return List with stability metrics and credible intervals
#' @export
calculate_stability <- function(fit, parameter, credible_interval = 0.95) {
  
  if (!inherits(fit, "mfeargen_hierarchical")) {
    stop("Input must be a mfeargen_hierarchical object")
  }
  
  # Extract variance components
  tau_sq_name <- paste0("tau_", parameter, "_sq")
  sigma_sq_name <- paste0("sigma_", parameter, "_sq")
  
  mcmc_samples <- fit$jags_fit$samples
  
  # Combine chains
  tau_sq_samples <- unlist(lapply(mcmc_samples, function(x) x[[tau_sq_name]]))
  sigma_sq_samples <- unlist(lapply(mcmc_samples, function(x) x[[sigma_sq_name]]))
  
  # Calculate stability for each MCMC sample
  stability_samples <- tau_sq_samples / (tau_sq_samples + sigma_sq_samples)
  
  # Calculate summary statistics
  alpha <- (1 - credible_interval) / 2
  
  results <- list(
    parameter = parameter,
    mean_stability = mean(stability_samples),
    median_stability = median(stability_samples),
    sd_stability = sd(stability_samples),
    ci_lower = quantile(stability_samples, alpha),
    ci_upper = quantile(stability_samples, 1 - alpha),
    between_variance = mean(tau_sq_samples),
    within_variance = mean(sigma_sq_samples),
    variance_ratio = mean(tau_sq_samples) / mean(sigma_sq_samples),
    samples = stability_samples,
    interpretation = interpret_stability_value(mean(stability_samples))
  )
  
  class(results) <- c("mfeargen_stability", "list")
  return(results)
}

#' Interpret Stability Values
#' @keywords internal
interpret_stability_value <- function(stability) {
  if (stability >= 0.7) {
    interpretation <- list(
      level = "High",
      description = "Strong individual differences, minimal session variation",
      implication = "Parameter represents a stable trait-like characteristic"
    )
  } else if (stability >= 0.3) {
    interpretation <- list(
      level = "Moderate",
      description = "Balanced individual differences and session variation",
      implication = "Parameter shows both trait and state components"
    )
  } else {
    interpretation <- list(
      level = "Low",
      description = "Minimal individual differences, high session variation",
      implication = "Parameter is primarily state-dependent"
    )
  }
  return(interpretation)
}

#' Interpret Stability Results
#' @keywords internal
interpret_stability <- function(stability_list) {
  interpretations <- list()
  
  for (param in names(stability_list)) {
    value <- stability_list[[param]]
    interpretations[[param]] <- interpret_stability_value(value)
  }
  
  return(interpretations)
}

#' Extract Individual Stable Components
#'
#' @description
#' Extracts the stable (characteristic) parameter values for each individual,
#' representing their average tendency across sessions.
#'
#' @param fit mfeargen_hierarchical object
#' @param parameter Character. Parameter to extract
#' @return Matrix of individual stable values (participants x MCMC samples)
#' @export
extract_stable_components <- function(fit, parameter = "alpha") {
  
  if (!inherits(fit, "mfeargen_hierarchical")) {
    stop("Input must be a mfeargen_hierarchical object")
  }
  
  stable_name <- paste0(parameter, "_stable")
  
  # Extract from MCMC samples
  mcmc_samples <- fit$jags_fit$samples
  n_chains <- length(mcmc_samples)
  
  # Combine chains
  stable_samples <- NULL
  for (chain in 1:n_chains) {
    if (stable_name %in% names(mcmc_samples[[chain]])) {
      chain_samples <- mcmc_samples[[chain]][[stable_name]]
      if (is.null(stable_samples)) {
        stable_samples <- chain_samples
      } else {
        stable_samples <- rbind(stable_samples, chain_samples)
      }
    }
  }
  
  if (is.null(stable_samples)) {
    stop(paste("Stable component for", parameter, "not found in model output"))
  }
  
  return(stable_samples)
}

#' Extract Session Deviations
#'
#' @description
#' Extracts session-specific deviations from individual stable values,
#' quantifying within-person variation across sessions.
#'
#' @param fit mfeargen_hierarchical object
#' @param parameter Character. Parameter to analyze
#' @return Array of deviations (participants x sessions x MCMC samples)
#' @export
extract_session_deviations <- function(fit, parameter = "alpha") {
  
  if (!inherits(fit, "mfeargen_hierarchical")) {
    stop("Input must be a mfeargen_hierarchical object")
  }
  
  # Get stable components
  stable <- extract_stable_components(fit, parameter)
  
  # Get session-specific values
  session_name <- paste0(parameter, "_session")
  mcmc_samples <- fit$jags_fit$samples
  
  # Calculate deviations
  # This will depend on the specific model structure
  # Placeholder for actual implementation
  
  message("Extracting session deviations for ", parameter)
  
  # Return placeholder
  return(NULL)
}

#' Plot Hierarchical Structure
#'
#' @description
#' Creates a visual representation of the three-level hierarchical model
#' structure, showing population, individual, and session levels.
#'
#' @param fit mfeargen_hierarchical object
#' @param parameter Character. Parameter to visualize
#' @param n_individuals Integer. Number of individuals to show (default: 10)
#' @return ggplot object
#' @export
#' @import ggplot2
#' @import dplyr
plot_hierarchical_structure <- function(fit, 
                                       parameter = "alpha",
                                       n_individuals = 10) {
  
  if (!inherits(fit, "mfeargen_hierarchical")) {
    stop("Input must be a mfeargen_hierarchical object")
  }
  
  # Extract components
  var_decomp <- calculate_variance_decomposition(fit, parameters = parameter)
  
  # Create visualization data
  set.seed(123)  # For reproducibility of example
  
  # Population level
  pop_mean <- fit$population_parameters[[paste0("mu_", parameter)]]
  pop_sd <- sqrt(var_decomp$between_variance[[parameter]])
  
  # Individual stable values (sample for visualization)
  n_show <- min(n_individuals, fit$data_info$n_participants)
  ind_stable <- rnorm(n_show, pop_mean, pop_sd)
  
  # Session realizations (example with 3 sessions per person)
  session_sd <- sqrt(var_decomp$within_variance[[parameter]])
  session_data <- data.frame()
  
  for (i in 1:n_show) {
    for (s in 1:3) {
      session_data <- rbind(session_data, data.frame(
        individual = i,
        session = s,
        stable_value = ind_stable[i],
        session_value = rnorm(1, ind_stable[i], session_sd),
        parameter = parameter
      ))
    }
  }
  
  # Create plot
  p <- ggplot(session_data, aes(x = factor(individual))) +
    # Population level
    geom_hline(yintercept = pop_mean, color = "red", size = 1.5, alpha = 0.7) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, 
                  ymin = pop_mean - 2*pop_sd, 
                  ymax = pop_mean + 2*pop_sd),
              fill = "red", alpha = 0.1) +
    
    # Individual stable components
    geom_point(aes(y = stable_value), size = 4, color = "blue", shape = 18) +
    
    # Session realizations
    geom_point(aes(y = session_value, shape = factor(session)), 
               size = 2, alpha = 0.6) +
    geom_segment(aes(xend = factor(individual), 
                     y = stable_value, 
                     yend = session_value),
                 alpha = 0.3) +
    
    # Formatting
    labs(
      title = paste("Three-Level Hierarchical Structure:", parameter),
      subtitle = sprintf("Stability = %.2f (Between/Total Variance)", 
                        var_decomp$stability[[parameter]]),
      x = "Individual",
      y = paste(parameter, "value"),
      shape = "Session"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      legend.position = "right"
    ) +
    
    # Add annotations
    annotate("text", x = 0.5, y = pop_mean, 
             label = "Population Mean", 
             color = "red", hjust = 0) +
    annotate("text", x = n_show + 0.5, 
             y = max(session_data$session_value) * 0.95,
             label = paste("Between-person var:", round(pop_sd^2, 3)),
             size = 3, hjust = 1) +
    annotate("text", x = n_show + 0.5, 
             y = max(session_data$session_value) * 0.90,
             label = paste("Within-person var:", round(session_sd^2, 3)),
             size = 3, hjust = 1)
  
  return(p)
}

#' Compare Stability Across Parameters
#'
#' @description
#' Compares stability metrics across multiple parameters to identify
#' which aspects of behavior are more trait-like vs state-like.
#'
#' @param fit mfeargen_hierarchical object
#' @param parameters Character vector of parameters to compare
#' @return Data frame with stability comparisons
#' @export
compare_parameter_stability <- function(fit, 
                                       parameters = c("alpha", "lambda", "w0", "w1")) {
  
  if (!inherits(fit, "mfeargen_hierarchical")) {
    stop("Input must be a mfeargen_hierarchical object")
  }
  
  # Calculate stability for each parameter
  stability_results <- list()
  
  for (param in parameters) {
    tryCatch({
      stability_results[[param]] <- calculate_stability(fit, param)
    }, error = function(e) {
      message("Could not calculate stability for ", param, ": ", e$message)
    })
  }
  
  # Create comparison table
  comparison <- data.frame(
    Parameter = character(),
    Stability = numeric(),
    CI_Lower = numeric(),
    CI_Upper = numeric(),
    Between_Var = numeric(),
    Within_Var = numeric(),
    Interpretation = character(),
    stringsAsFactors = FALSE
  )
  
  for (param in names(stability_results)) {
    res <- stability_results[[param]]
    comparison <- rbind(comparison, data.frame(
      Parameter = param,
      Stability = res$mean_stability,
      CI_Lower = res$ci_lower,
      CI_Upper = res$ci_upper,
      Between_Var = res$between_variance,
      Within_Var = res$within_variance,
      Interpretation = res$interpretation$level,
      stringsAsFactors = FALSE
    ))
  }
  
  # Sort by stability
  comparison <- comparison[order(comparison$Stability, decreasing = TRUE), ]
  
  return(comparison)
}

#' Plot Stability Comparison
#'
#' @description
#' Creates a forest plot comparing stability metrics across parameters
#' with credible intervals.
#'
#' @param stability_comparison Output from compare_parameter_stability()
#' @return ggplot object
#' @export
#' @import ggplot2
plot_stability_comparison <- function(stability_comparison) {
  
  # Prepare data
  stability_comparison$Parameter <- factor(
    stability_comparison$Parameter,
    levels = stability_comparison$Parameter[order(stability_comparison$Stability)]
  )
  
  # Create forest plot
  p <- ggplot(stability_comparison, 
              aes(x = Stability, y = Parameter)) +
    
    # Reference lines
    geom_vline(xintercept = c(0.3, 0.7), 
               linetype = "dashed", alpha = 0.3) +
    
    # Credible intervals
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper),
                   height = 0.2) +
    
    # Point estimates
    geom_point(aes(color = Interpretation), size = 4) +
    
    # Formatting
    scale_color_manual(values = c("High" = "darkgreen",
                                  "Moderate" = "orange",
                                  "Low" = "red")) +
    scale_x_continuous(limits = c(0, 1), 
                       breaks = seq(0, 1, 0.2)) +
    
    labs(
      title = "Parameter Stability Comparison",
      subtitle = "Proportion of variance due to stable individual differences",
      x = "Stability (ICC)",
      y = "Parameter",
      color = "Stability Level"
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      legend.position = "right"
    ) +
    
    # Add region labels
    annotate("text", x = 0.15, y = 0.5, 
             label = "State-like", alpha = 0.5) +
    annotate("text", x = 0.5, y = 0.5, 
             label = "Mixed", alpha = 0.5) +
    annotate("text", x = 0.85, y = 0.5, 
             label = "Trait-like", alpha = 0.5)
  
  return(p)
}

#' S3 Print Method for Variance Decomposition
#' @param x A mfeargen_variance_decomposition object
#' @param ... Additional arguments (not used)
#' @export
print.mfeargen_variance_decomposition <- function(x, ...) {
  cat("Variance Decomposition Results\n")
  cat("==============================\n\n")
  
  for (param in names(x$stability)) {
    cat(sprintf("%s:\n", param))
    cat(sprintf("  Between-person variance: %.4f\n", x$between_variance[[param]]))
    cat(sprintf("  Within-person variance:  %.4f\n", x$within_variance[[param]]))
    cat(sprintf("  Stability (ICC):         %.3f\n", x$stability[[param]]))
    cat(sprintf("  Interpretation:          %s\n", x$interpretation[[param]]$level))
    cat("\n")
  }
}

#' S3 Print Method for Stability Results
#' @param x A mfeargen_stability object
#' @param ... Additional arguments (not used)
#' @export
print.mfeargen_stability <- function(x, ...) {
  cat(sprintf("Stability Analysis for %s\n", x$parameter))
  cat("================================\n\n")
  cat(sprintf("Mean Stability:     %.3f\n", x$mean_stability))
  cat(sprintf("95%% CI:            [%.3f, %.3f]\n", x$ci_lower, x$ci_upper))
  cat(sprintf("Between Variance:   %.4f\n", x$between_variance))
  cat(sprintf("Within Variance:    %.4f\n", x$within_variance))
  cat(sprintf("Variance Ratio:     %.2f\n", x$variance_ratio))
  cat(sprintf("\nInterpretation:     %s\n", x$interpretation$level))
  cat(sprintf("Description:        %s\n", x$interpretation$description))
  cat(sprintf("Implication:        %s\n", x$interpretation$implication))
}