#' @title Additional S3 Methods for mfeargen Objects
#' @description Additional generic methods for mfeargen class objects
#' @name s3_methods_additional
NULL

#' Extract Coefficients
#'
#' @param object mfeargen object
#' @param ... Additional arguments
#' @return Data frame of coefficients
#' @export
coef.mfeargen <- function(object, ...) {
  # Extract parameters more carefully
  if (!inherits(object, "mfeargen") || is.null(object$jags_fit)) {
    return(NULL)
  }
  
  jags_output <- object$jags_fit
  coef_table <- data.frame()
  
  # Extract key parameters
  param_names <- c("alpha", "lambda", "sigma", "pi")
  
  for (param_name in param_names) {
    if (!is.null(jags_output$mean[[param_name]])) {
      param_values <- jags_output$mean[[param_name]]
      param_sd <- jags_output$sd[[param_name]]
      param_lower <- if (!is.null(jags_output$q2.5)) jags_output$q2.5[[param_name]] else rep(NA, length(param_values))
      param_upper <- if (!is.null(jags_output$q97.5)) jags_output$q97.5[[param_name]] else rep(NA, length(param_values))
      
      n_values <- length(param_values)
      
      # Handle both participant-level and group-level parameters
      if (param_name == "pi") {
        # Group probabilities
        coef_table <- rbind(coef_table, data.frame(
          parameter = paste0("pi[", 1:n_values, "]"),
          estimate = param_values,
          std_error = param_sd,
          ci_lower = param_lower,
          ci_upper = param_upper,
          stringsAsFactors = FALSE
        ))
      } else {
        # Participant-level parameters
        coef_table <- rbind(coef_table, data.frame(
          parameter = paste0(param_name, "[", 1:n_values, "]"),
          estimate = param_values,
          std_error = param_sd,
          ci_lower = param_lower,
          ci_upper = param_upper,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  if (nrow(coef_table) == 0) {
    return(NULL)
  }
  
  return(coef_table)
}

#' Print Method for Recovery Results
#'
#' @param x mfeargen_recovery object
#' @param ... Additional arguments
#' @export
print.mfeargen_recovery <- function(x, ...) {
  cat("mfeargen Parameter Recovery Analysis\n")
  cat("====================================\n")
  cat("Participants:", x$settings$n_participants, "\n")
  cat("Model type:", x$settings$model_type, "\n")
  cat("MCMC iterations:", x$settings$n.iter, "\n\n")
  
  if (!is.null(x$metrics$alpha_correlation)) {
    cat("Recovery Correlations:\n")
    cat(sprintf("  Alpha:  r = %.3f\n", x$metrics$alpha_correlation))
    cat(sprintf("  Lambda: r = %.3f\n", x$metrics$lambda_correlation))
  }
  
  if (!is.null(x$metrics$group_accuracy)) {
    cat(sprintf("\nGroup Accuracy: %.1f%%\n", x$metrics$group_accuracy * 100))
  }
  
  invisible(NULL)
}