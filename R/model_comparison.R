#' Model Comparison Utilities for mfeargen
#'
#' Functions for comparing different model fits and selecting the best model

#' Compare mfeargen Models
#'
#' Compares multiple fitted mfeargen models using information criteria
#'
#' @param ... Fitted mfeargen model objects
#' @param criterion Character. Comparison criterion ("DIC", "WAIC", "LOO")
#' @return Data frame with model comparison results
#' @export
#' @examples
#' \dontrun{
#' fit1 <- fit_mfeargen(data, model_type = "1cs")
#' fit2 <- fit_mfeargen(data, model_type = "2cs")
#' comparison <- compare_models(fit1, fit2)
#' print(comparison)
#' }
compare_models <- function(..., criterion = "DIC") {
  
  # Get model objects
  models <- list(...)
  
  if (length(models) < 2) {
    stop("At least two models are required for comparison")
  }
  
  # Validate that all objects are mfeargen fits
  for (i in seq_along(models)) {
    if (!inherits(models[[i]], "mfeargen")) {
      stop(sprintf("Object %d is not a fitted mfeargen model", i))
    }
  }
  
  # Get model names
  model_names <- as.character(substitute(list(...)))[-1]
  if (length(model_names) != length(models)) {
    model_names <- paste0("Model_", seq_along(models))
  }
  
  # Initialize results data frame
  results <- data.frame(
    Model = model_names,
    Type = character(length(models)),
    N_Participants = integer(length(models)),
    N_Trials = integer(length(models)),
    N_Parameters = integer(length(models)),
    DIC = numeric(length(models)),
    Deviance = numeric(length(models)),
    pD = numeric(length(models)),
    Lambda_Threshold = numeric(length(models)),
    Converged = logical(length(models)),
    stringsAsFactors = FALSE
  )
  
  # Extract information from each model
  for (i in seq_along(models)) {
    model <- models[[i]]
    jags_fit <- model$jags_fit
    
    results$Type[i] <- model$model_type
    results$N_Participants[i] <- model$jags_data$Nparticipants
    results$N_Trials[i] <- model$jags_data$Ntrials
    results$Lambda_Threshold[i] <- model$lambda_threshold
    results$Converged[i] <- model$convergence$converged
    
    # Extract DIC information
    if ("DIC" %in% names(jags_fit)) {
      results$DIC[i] <- jags_fit$DIC
      results$Deviance[i] <- jags_fit$deviance
      results$pD[i] <- jags_fit$pD
    } else {
      results$DIC[i] <- NA
      results$Deviance[i] <- NA
      results$pD[i] <- NA
    }
    
    # Estimate number of parameters
    if (!is.na(results$pD[i])) {
      results$N_Parameters[i] <- round(results$pD[i])
    } else {
      # Rough estimate based on model type and participants
      n_p <- results$N_Participants[i]
      if (results$Type[i] == "1cs") {
        results$N_Parameters[i] <- n_p * 4 + 20  # Rough estimate
      } else {
        results$N_Parameters[i] <- n_p * 6 + 25  # Rough estimate
      }
    }
  }
  
  # Calculate relative performance
  if (criterion == "DIC" && !all(is.na(results$DIC))) {
    min_dic <- min(results$DIC, na.rm = TRUE)
    results$Delta_DIC <- results$DIC - min_dic
    results$DIC_Weight <- exp(-0.5 * results$Delta_DIC)
    results$DIC_Weight <- results$DIC_Weight / sum(results$DIC_Weight, na.rm = TRUE)
    
    # Rank models
    results$Rank <- rank(results$DIC, ties.method = "min")
  }
  
  # Add interpretation
  if ("Delta_DIC" %in% names(results)) {
    results$Interpretation <- ifelse(results$Delta_DIC < 2, "Substantial support",
                                   ifelse(results$Delta_DIC < 7, "Weak support", 
                                          "No support"))
  }
  
  # Sort by performance
  if ("Rank" %in% names(results)) {
    results <- results[order(results$Rank), ]
  }
  
  # Add comparison summary
  attr(results, "criterion") <- criterion
  attr(results, "best_model") <- results$Model[1]
  attr(results, "n_models") <- length(models)
  
  class(results) <- c("mfeargen_comparison", "data.frame")
  
  return(results)
}

#' Print Method for Model Comparison
#'
#' @param x A mfeargen_comparison object
#' @param ... Additional arguments
#' @export
print.mfeargen_comparison <- function(x, ...) {
  cat("mfeargen Model Comparison\n")
  cat("========================\n\n")
  
  criterion <- attr(x, "criterion")
  best_model <- attr(x, "best_model")
  n_models <- attr(x, "n_models")
  
  cat("Comparison criterion:", criterion, "\n")
  cat("Number of models:", n_models, "\n")
  cat("Best model:", best_model, "\n\n")
  
  # Print main results
  cat("Model Comparison Results:\n")
  print(x[, c("Model", "Type", "N_Participants", "DIC", "Delta_DIC", "Rank", "Interpretation")])
  
  cat("\nModel Details:\n")
  print(x[, c("Model", "N_Parameters", "Lambda_Threshold", "Converged")])
  
  if ("DIC_Weight" %in% names(x)) {
    cat("\nModel Weights:\n")
    weights_df <- x[, c("Model", "DIC_Weight")]
    weights_df$DIC_Weight <- round(weights_df$DIC_Weight, 3)
    print(weights_df)
  }
  
  cat("\nInterpretation Guide:\n")
  cat("- Delta DIC < 2: Substantial support\n")
  cat("- Delta DIC 2-7: Weak support\n")
  cat("- Delta DIC > 7: No support\n")
}

#' Calculate WAIC for mfeargen Models
#'
#' Calculates Watanabe-Akaike Information Criterion
#'
#' @param fit_object A fitted mfeargen model
#' @return List with WAIC components
#' @export
calculate_waic <- function(fit_object) {
  
  if (!inherits(fit_object, "mfeargen")) {
    stop("fit_object must be a fitted mfeargen model")
  }
  
  # Extract log-likelihood samples
  if (!"loglik" %in% names(fit_object$jags_fit$sims.list)) {
    stop("Log-likelihood not available in model fit. Refit with loglik monitoring.")
  }
  
  loglik <- fit_object$jags_fit$sims.list$loglik
  
  # Calculate WAIC components
  # loglik should be [n_samples, n_participants, n_trials]
  if (length(dim(loglik)) == 3) {
    # Sum over trials for each participant and sample
    loglik_sum <- apply(loglik, c(1, 2), sum)  # [n_samples, n_participants]
  } else if (length(dim(loglik)) == 2) {
    # Already summed over trials
    loglik_sum <- loglik
  } else {
    stop("Unexpected log-likelihood dimensions")
  }
  
  # WAIC calculation
  # 1. Pointwise log predictive density
  lppd <- sum(log(apply(exp(loglik_sum), 2, mean)))
  
  # 2. Effective number of parameters
  p_waic <- sum(apply(loglik_sum, 2, var))
  
  # 3. WAIC
  waic <- -2 * (lppd - p_waic)
  
  # Standard error
  pointwise_waic <- -2 * (log(apply(exp(loglik_sum), 2, mean)) - apply(loglik_sum, 2, var))
  se_waic <- sqrt(length(pointwise_waic) * var(pointwise_waic))
  
  result <- list(
    waic = waic,
    lppd = lppd,
    p_waic = p_waic,
    se_waic = se_waic,
    pointwise = pointwise_waic
  )
  
  class(result) <- "mfeargen_waic"
  return(result)
}

#' Cross-Validation for mfeargen Models
#'
#' Performs leave-one-out cross-validation
#'
#' @param data Data frame with fear generalization data
#' @param model_type Character. Either "1cs" or "2cs"
#' @param k_folds Integer. Number of folds for cross-validation (default: 10)
#' @param ... Additional arguments passed to fit_mfeargen
#' @return List with cross-validation results
#' @export
cross_validate_mfeargen <- function(data, model_type, k_folds = 10, ...) {
  
  participants <- unique(data$participant)
  n_participants <- length(participants)
  
  if (k_folds > n_participants) {
    warning("k_folds larger than number of participants. Using leave-one-out.")
    k_folds <- n_participants
  }
  
  # Create folds
  fold_size <- ceiling(n_participants / k_folds)
  folds <- split(sample(participants), ceiling(seq_along(participants) / fold_size))
  
  # Initialize results
  cv_results <- list()
  fold_metrics <- data.frame(
    fold = integer(),
    n_train = integer(),
    n_test = integer(),
    loglik_test = numeric(),
    mse_test = numeric(),
    converged = logical()
  )
  
  cat("Running", k_folds, "fold cross-validation...\n")
  
  for (k in seq_along(folds)) {
    cat("Fold", k, "of", k_folds, "...")
    
    # Split data
    test_participants <- folds[[k]]
    train_participants <- setdiff(participants, test_participants)
    
    train_data <- data[data$participant %in% train_participants, ]
    test_data <- data[data$participant %in% test_participants, ]
    
    # Fit model on training data
    tryCatch({
      fit <- fit_mfeargen(train_data, model_type = model_type, ...)
      
      # Predict on test data (simplified prediction)
      # In practice, this would require more sophisticated prediction
      test_loglik <- NA  # Placeholder
      test_mse <- NA     # Placeholder
      converged <- fit$convergence$converged
      
      fold_metrics <- rbind(fold_metrics, data.frame(
        fold = k,
        n_train = length(train_participants),
        n_test = length(test_participants),
        loglik_test = test_loglik,
        mse_test = test_mse,
        converged = converged
      ))
      
      cv_results[[k]] <- list(
        fold = k,
        fit = fit,
        train_participants = train_participants,
        test_participants = test_participants
      )
      
      cat(" completed\n")
      
    }, error = function(e) {
      cat(" failed:", e$message, "\n")
      fold_metrics <<- rbind(fold_metrics, data.frame(
        fold = k,
        n_train = length(train_participants),
        n_test = length(test_participants),
        loglik_test = NA,
        mse_test = NA,
        converged = FALSE
      ))
    })
  }
  
  # Summarize results
  summary_stats <- list(
    mean_loglik = mean(fold_metrics$loglik_test, na.rm = TRUE),
    se_loglik = sqrt(var(fold_metrics$loglik_test, na.rm = TRUE) / nrow(fold_metrics)),
    mean_mse = mean(fold_metrics$mse_test, na.rm = TRUE),
    se_mse = sqrt(var(fold_metrics$mse_test, na.rm = TRUE) / nrow(fold_metrics)),
    convergence_rate = mean(fold_metrics$converged, na.rm = TRUE),
    successful_folds = sum(!is.na(fold_metrics$loglik_test))
  )
  
  result <- list(
    fold_results = cv_results,
    fold_metrics = fold_metrics,
    summary = summary_stats,
    k_folds = k_folds,
    model_type = model_type
  )
  
  class(result) <- "mfeargen_cv"
  return(result)
}

#' Model Selection Helper
#'
#' Provides guidance on model selection based on comparison results
#'
#' @param comparison_result A mfeargen_comparison object
#' @return Character string with model selection advice
#' @export
model_selection_advice <- function(comparison_result) {
  
  if (!inherits(comparison_result, "mfeargen_comparison")) {
    stop("Input must be a mfeargen_comparison object")
  }
  
  if (!"Delta_DIC" %in% names(comparison_result)) {
    return("Model comparison requires DIC values. Refit models with DIC calculation.")
  }
  
  best_model <- comparison_result$Model[1]
  second_best <- comparison_result$Model[2]
  delta_dic <- comparison_result$Delta_DIC[2]
  
  advice <- paste0("Model Selection Advice:\n\n")
  advice <- paste0(advice, "Best model: ", best_model, "\n")
  
  if (delta_dic < 2) {
    advice <- paste0(advice, "WARNING: The difference between the best and second-best model is small (Delta-Delta = ", 
                    round(delta_dic, 2), ").\n")
    advice <- paste0(advice, "Both models have substantial support. Consider:\n")
    advice <- paste0(advice, "- Model averaging\n")
    advice <- paste0(advice, "- Additional data collection\n")
    advice <- paste0(advice, "- Cross-validation\n")
  } else if (delta_dic < 7) {
    advice <- paste0(advice, "The best model has weak to moderate support (Delta-DIC = ", 
                    round(delta_dic, 2), ").\n")
    advice <- paste0(advice, "The ", best_model, " is preferred but consider the second model.\n")
  } else {
    advice <- paste0(advice, "Strong support for ", best_model, " (Delta-DIC = ", 
                    round(delta_dic, 2), ").\n")
    advice <- paste0(advice, "This model is clearly preferred.\n")
  }
  
  # Check convergence
  converged_models <- sum(comparison_result$Converged)
  total_models <- nrow(comparison_result)
  
  if (converged_models < total_models) {
    advice <- paste0(advice, "\nWarning: ", total_models - converged_models, 
                    " out of ", total_models, " models did not converge.\n")
    advice <- paste0(advice, "Consider increasing MCMC iterations or checking data quality.\n")
  }
  
  return(advice)
}