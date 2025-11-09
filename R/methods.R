# S3 Methods for mfeargen Objects

#' Print Method for mfeargen Objects
#'
#' @param x A mfeargen object
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the input object
#' @importFrom stats fitted
#' @export
print.mfeargen <- function(x, ...) {
  
  model_type <- if (!is.null(x$model_type)) x$model_type else NA
  n_participants <- if (!is.null(x$jags_data$Nparticipants)) x$jags_data$Nparticipants else NA
  n_trials <- if (!is.null(x$jags_data$Ntrials)) x$jags_data$Ntrials else NA
  lambda_threshold <- if (!is.null(x$lambda_threshold)) round(x$lambda_threshold, 4) else NA
  n_iter <- if (!is.null(x$n.iter)) x$n.iter else NA
  n_burn <- if (!is.null(x$n.burnin)) x$n.burnin else NA
  n_chains <- if (!is.null(x$n.chains)) x$n.chains else NA
  convergence_info <- x$convergence

  cat("mfeargen Bayesian Mixture Model\n")
  cat("===============================\n\n")
  
  # Basic model info
  cat("Model Type:", model_type, "\n")
  cat("Participants:", n_participants, "\n")
  cat("Trials:", n_trials, "\n")
  cat("Lambda Threshold:", lambda_threshold, "\n")
  cat("MCMC Iterations:", n_iter, "(burnin:", n_burn, ")\n")
  cat("Chains:", n_chains, "\n\n")
  
  # Convergence info
  if (!is.null(convergence_info)) {
    if (!is.null(convergence_info$converged) && convergence_info$converged) {
      cat("Convergence: GOOD\n")
    } else if (!is.null(convergence_info$converged) && !convergence_info$converged) {
      cat("Convergence: ISSUES DETECTED\n")
      if (!is.null(convergence_info$max_rhat)) {
        cat("  Max R-hat:", round(convergence_info$max_rhat, 3), "\n")
      }
      if (!is.null(convergence_info$min_neff)) {
        cat("  Min n.eff:", round(convergence_info$min_neff, 0), "\n")
      }
    } else {
      cat("Convergence: Not checked\n")
    }
  } else {
    cat("Convergence: Not available\n")
  }
  
  # Group membership summary
  if (!is.null(x$group_membership)) {
    cat("\nGroup Membership:\n")
    group_counts <- table(x$group_membership$group_label)
    for (group in names(group_counts)) {
      cat("  ", group, ":", group_counts[group], "\n")
    }
  }
  
  # Model fit
  if (!is.null(x$jags_fit) && "DIC" %in% names(x$jags_fit)) {
    cat("\nModel Fit:\n")
    cat("  DIC:", round(x$jags_fit$DIC, 2), "\n")
  }
  
  cat("\nUse summary() for detailed parameter estimates\n")
  cat("Use plot() for visualizations\n")
  
  invisible(x)
}

#' Summary Method for mfeargen Objects
#' 
#' @param object A mfeargen object
#' @param ... Additional arguments (ignored)
#' @return A summary object
#' @export
summary.mfeargen <- function(object, ...) {
  
  # Extract key parameters
  params <- extract_parameters(object)
  
  # Group probabilities
  pi_summary <- object$jags_fit$summary[grepl("^pi\\[", rownames(object$jags_fit$summary)), ]
  
  # Create summary object
  summary_obj <- list(
    call = object$call,
    model_type = object$model_type,
    n_participants = object$jags_data$Nparticipants,
    n_trials = object$jags_data$Ntrials,
    lambda_threshold = object$lambda_threshold,
    group_membership = object$group_membership,
    group_probabilities = pi_summary,
    convergence = object$convergence,
    dic = object$jags_fit$DIC
  )
  
  class(summary_obj) <- "summary.mfeargen"
  return(summary_obj)
}

#' Print Method for mfeargen Summary Objects
#' 
#' @param x A summary.mfeargen object
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the input object
#' @export
print.summary.mfeargen <- function(x, ...) {
  
  cat("mfeargen Model Summary\n")
  cat("=====================\n\n")
  
  # Call
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  # Model details
  cat("Model Type:", x$model_type, "\n")
  cat("Participants:", x$n_participants, "\n")
  cat("Trials:", x$n_trials, "\n")
  cat("Lambda Threshold:", round(x$lambda_threshold, 4), "\n\n")
  
  # Group probabilities
  cat("Group Probabilities:\n")
  print(round(x$group_probabilities, 3))
  cat("\n")
  
  # Group membership counts
  cat("Group Membership (threshold = 0.75):\n")
  group_counts <- table(x$group_membership$assigned_group)
  for (group in names(group_counts)) {
    pct <- round(100 * group_counts[group] / x$n_participants, 1)
    cat("  ", group, ":", group_counts[group], "(", pct, "%)\n")
  }
  
  # Convergence
  cat("\nConvergence:\n")
  if (x$convergence$convergence_ok) {
    cat("  Status: GOOD\n")
  } else {
    cat("  Status: ISSUES\n")
    cat("  Max R-hat:", round(x$convergence$max_rhat, 3), "\n")
    cat("  Min n.eff:", round(x$convergence$min_neff, 0), "\n")
  }
  
  # Model fit
  cat("\nModel Fit:\n")
  cat("  DIC:", round(x$dic, 2), "\n")
  
  invisible(x)
}


#' Coefficients for mfeargen Objects
#'
#' @param object A mfeargen object
#' @param ... Additional arguments (ignored)
#' @return Named list of parameter summaries
#' @export
coef.mfeargen <- function(object, ...) {
  if (!inherits(object, "mfeargen")) {
    stop("Input must be a mfeargen object")
  }
  summary_mat <- object$jags_fit$summary
  if (is.null(summary_mat)) {
    return(list())
  }
  extract_block <- function(pattern) {
    rows <- grepl(pattern, rownames(summary_mat))
    if (!any(rows)) return(NULL)
    summary_mat[rows, , drop = FALSE]
  }
  result <- list(
    alpha = extract_block("^alpha\\["),
    lambda = extract_block("^lambda\\["),
    sigma = extract_block("^sigma\\[")
  )
  Filter(Negate(is.null), result)
}


#' Predict Method for mfeargen Objects
#' 
#' @param object A mfeargen object
#' @param newdata Optional new data for prediction
#' @param type Character. Type of prediction ("groups", "response")
#' @param ... Additional arguments
#' @return Predictions based on type
#' @export
predict.mfeargen <- function(object, newdata = NULL, type = "groups", ...) {
  
  switch(type,
    "groups" = predict_groups(object, ...),
    "response" = {
      preds <- fitted(object)
      if (is.null(preds)) {
        stop("Posterior predictions are not available for this fit. ",
             "Re-run the model with posterior predictive monitoring enabled.")
      }
      preds
    },
    stop("Unknown prediction type: ", type, 
         ". Options: 'groups', 'response'")
  )
}

#' Residuals Method for mfeargen Objects
#' 
#' @param object A mfeargen object
#' @param type Character. Type of residuals ("raw", "standardized")
#' @param ... Additional arguments (ignored)
#' @return Matrix of residuals
#' @export
residuals.mfeargen <- function(object, type = "raw", ...) {
  
  # Extract observed data
  y_obs <- object$jags_data$y
  
  # Get fitted values using the fitted method
  theta_matrix <- fitted(object)
  
  if (is.null(theta_matrix)) {
    warning("Cannot compute residuals: no fitted values available")
    return(NULL)
  }
  
  # Ensure dimensions match
  if (!all(dim(y_obs) == dim(theta_matrix))) {
    warning("Dimension mismatch between observed and fitted values")
    return(NULL)
  }
  
  # Calculate residuals
  residuals <- y_obs - theta_matrix
  
  if (type == "standardized") {
    # Try to standardize by posterior SD
    if (!is.null(object$jags_fit$mean) && "sigma" %in% names(object$jags_fit$mean)) {
      sigma <- object$jags_fit$mean$sigma
      # Expand sigma to match residuals dimensions if needed
      if (length(sigma) == object$jags_data$Nparticipants) {
        sigma_matrix <- matrix(rep(sigma, each = object$jags_data$Ntrials), 
                              nrow = object$jags_data$Nparticipants, 
                              ncol = object$jags_data$Ntrials, 
                              byrow = TRUE)
        residuals <- residuals / sigma_matrix
      }
    } else if (!is.null(object$jags_fit$summary)) {
      sigma_rows <- grepl("^sigma\\[", rownames(object$jags_fit$summary))
      if (any(sigma_rows)) {
        sigma <- object$jags_fit$summary[sigma_rows, "mean"]
        if (length(sigma) == object$jags_data$Nparticipants) {
          sigma_matrix <- matrix(rep(sigma, each = object$jags_data$Ntrials), 
                                nrow = object$jags_data$Nparticipants, 
                                ncol = object$jags_data$Ntrials, 
                                byrow = TRUE)
          residuals <- residuals / sigma_matrix
        }
      }
    }
  }
  
  return(residuals)
}

#' Fitted Values Method for mfeargen Objects
#' 
#' @param object A mfeargen object
#' @param ... Additional arguments (ignored)
#' @return Matrix of fitted values
#' @export
fitted.mfeargen <- function(object, ...) {
  
  # Check if we have y_pre (posterior predictive samples)
  if (!is.null(object$jags_fit$sims.list) && "y_pre" %in% names(object$jags_fit$sims.list)) {
    # y_pre could be in different formats depending on JAGS model
    y_pre <- object$jags_fit$sims.list$y_pre
    
    # Check dimensions: could be (iterations, participants, trials) or (iterations, trials, participants)
    if (length(dim(y_pre)) == 3) {
      # If middle dimension matches Ntrials, it's likely (iter, trials, participants)
      if (dim(y_pre)[2] == object$jags_data$Ntrials && dim(y_pre)[3] == object$jags_data$Nparticipants) {
        # Need to transpose: (iter, trials, participants) -> (participants, trials)
        theta_matrix <- apply(y_pre, c(3, 2), mean)  # Average over iterations
      } else if (dim(y_pre)[2] == object$jags_data$Nparticipants && dim(y_pre)[3] == object$jags_data$Ntrials) {
        # Already in (iter, participants, trials) format
        theta_matrix <- apply(y_pre, c(2, 3), mean)  # Average over iterations
      } else {
        # Fallback: assume standard format
        theta_matrix <- apply(y_pre, c(2, 3), mean)
      }
      return(theta_matrix)
    }
  }
  
  # Try to extract theta from summary
  if (!is.null(object$jags_fit$summary)) {
    theta_rows <- grepl("^theta\\[", rownames(object$jags_fit$summary))
    if (any(theta_rows)) {
      theta <- object$jags_fit$summary[theta_rows, "mean"]
      # Create matrix (participants x trials)
      theta_matrix <- matrix(theta, 
                            nrow = object$jags_data$Nparticipants, 
                            ncol = object$jags_data$Ntrials, 
                            byrow = TRUE)
      return(theta_matrix)
    }
  }
  
  # If no theta, return the observed data as a fallback
  if (!is.null(object$jags_data$y)) {
    return(object$jags_data$y)
  }
  
  # No fitted values available
  return(NULL)
}

#' Log-likelihood Method for mfeargen Objects
#' 
#' @param object A mfeargen object
#' @param ... Additional arguments (ignored)
#' @return Log-likelihood value
#' @export
logLik.mfeargen <- function(object, ...) {
  
  # Extract log-likelihood from JAGS fit
  loglik <- object$jags_fit$summary["total_loglik", "mean"]
  
  # Create logLik object
  attr(loglik, "df") <- length(object$jags_fit$summary[, "mean"])  # Number of parameters
  attr(loglik, "nobs") <- object$jags_data$Nparticipants * object$jags_data$Ntrials
  
  class(loglik) <- "logLik"
  return(loglik)
}
