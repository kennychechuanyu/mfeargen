validate_scalar <- function(x, name, lower = -Inf, upper = Inf, integer = FALSE,
                            lower_open = FALSE, upper_open = FALSE) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
  if (valid) {
    valid <- (if (lower_open) x > lower else x >= lower) &&
      (if (upper_open) x < upper else x <= upper) &&
      (!integer || x == floor(x))
  }
  if (!valid) stop(name, " must be a finite ", if (integer) "integer" else "number",
                   " in ", if (lower_open) "(" else "[", lower, ", ", upper,
                   if (upper_open) ")" else "]", call. = FALSE)
  invisible(x)
}

validate_rating_scale <- function(scale) {
  if (!is.numeric(scale) || length(scale) != 2L || anyNA(scale) ||
      any(!is.finite(scale)) || scale[1] >= scale[2]) {
    stop("rating_scale must contain two finite, increasing bounds", call. = FALSE)
  }
  invisible(scale)
}

#' @noRd

# Internal helper functions (not exported)
#' @noRd
calculate_physical_distance <- function(stimulus_values, cs_value) {
  abs(stimulus_values - cs_value)
}

calculate_perceptual_distance <- function(perceived_values, perceived_cs_mean) {
  abs(perceived_values - perceived_cs_mean)
}

calculate_similarity <- function(distances, lambda) {
  exp(-lambda * distances)
}

# Internal helper function (not exported)
format_mcmc_summary <- function(mcmc_summary) {
  
  if (is.null(mcmc_summary)) {
    return(NULL)
  }
  
  # Convert to data frame if necessary
  if (is.matrix(mcmc_summary)) {
    mcmc_summary <- as.data.frame(mcmc_summary)
  }
  
  # Round numeric columns
  numeric_cols <- sapply(mcmc_summary, is.numeric)
  mcmc_summary[numeric_cols] <- lapply(mcmc_summary[numeric_cols], round, 4)
  
  # Add parameter names as row names if not already present
  if (!"parameter" %in% names(mcmc_summary)) {
    mcmc_summary$parameter <- rownames(mcmc_summary)
  }
  
  # Reorder columns
  col_order <- c("parameter", "mean", "sd", "2.5%", "50%", "97.5%", "Rhat", "n.eff")
  available_cols <- intersect(col_order, names(mcmc_summary))
  mcmc_summary <- mcmc_summary[, available_cols, drop = FALSE]
  
  return(mcmc_summary)
}

# Internal helper function (not exported)
check_mcmc_convergence <- function(rhat_values, threshold = 1.1) {
  
  if (is.null(rhat_values) || length(rhat_values) == 0) {
    return(list(
      converged = FALSE,
      n_converged = 0,
      n_total = 0,
      convergence_rate = 0,
      problematic_params = character(0)
    ))
  }
  
  # Remove NAs
  valid_rhat <- rhat_values[!is.na(rhat_values)]
  
  if (length(valid_rhat) == 0) {
    return(list(
      converged = FALSE,
      n_converged = 0,
      n_total = 0,
      convergence_rate = 0,
      problematic_params = character(0)
    ))
  }
  
  # Check convergence
  converged_params <- valid_rhat < threshold
  n_converged <- sum(converged_params)
  n_total <- length(valid_rhat)
  convergence_rate <- n_converged / n_total
  
  # Identify problematic parameters
  problematic_params <- names(valid_rhat)[!converged_params]
  if (is.null(problematic_params)) {
    problematic_params <- which(!converged_params)
  }
  
  overall_converged <- convergence_rate >= 0.95  # 95% of parameters must converge
  
  result <- list(
    converged = overall_converged,
    n_converged = n_converged,
    n_total = n_total,
    convergence_rate = convergence_rate,
    problematic_params = problematic_params,
    threshold = threshold
  )
  
  return(result)
}

# Internal helper function (not exported)
get_group_labels <- function(include_numbers = TRUE) {
  if (include_numbers) {
    return(c("1: Non-Learners", "2: Overgeneralizers", 
             "3: Physical Generalizers", "4: Perceptual Generalizers"))
  } else {
    return(c("Non-Learners", "Overgeneralizers", 
             "Physical Generalizers", "Perceptual Generalizers"))
  }
}

# Internal helper function (not exported)
create_progress_bar <- function(current, total, width = 50) {
  
  if (total <= 0) return("")
  
  progress <- min(current / total, 1)
  filled <- floor(progress * width)
  empty <- width - filled
  
  bar <- paste0("[", 
                paste(rep("=", filled), collapse = ""), 
                paste(rep("-", empty), collapse = ""), 
                "]")
  
  percentage <- round(progress * 100, 1)
  
  return(sprintf("%s %s%% (%d/%d)", bar, percentage, current, total))
}

# Internal helper function (not exported)
safe_divide <- function(numerator, denominator, default_value = 0) {
  ifelse(denominator == 0, default_value, numerator / denominator)
}

# Internal helper function (not exported)
scale_to_range <- function(x, new_min = 0, new_max = 1) {
  if (length(x) == 0) return(x)
  
  old_min <- min(x, na.rm = TRUE)
  old_max <- max(x, na.rm = TRUE)
  
  if (old_min == old_max) {
    return(rep(mean(c(new_min, new_max)), length(x)))
  }
  
  (x - old_min) / (old_max - old_min) * (new_max - new_min) + new_min
}

#' Calculate Lambda Threshold
#' 
#' Calculates the lambda threshold for exponential similarity function based on 
#' maximum distance and criterion similarity level.
#' 
#' @param max_distance Numeric. Maximum distance between stimuli
#' @param criterion Numeric. Criterion similarity level (between 0 and 1, exclusive)
#' @param stimulus_range Numeric vector of length 2. Range of stimulus values
#' 
#' @return Numeric. The calculated lambda threshold
#' 
#' @details
#' The lambda threshold is calculated using the formula:
#' lambda = -log(criterion) / max_distance
#' 
#' This ensures that at the maximum distance, the similarity equals the criterion.
#' 
#' @examples
#' # Calculate threshold for max distance of 10 with 70% similarity criterion
#' threshold <- calculate_lambda_threshold(max_distance = 10, criterion = 0.7)
#' 
#' # Calculate threshold using stimulus range
#' threshold <- calculate_lambda_threshold(stimulus_range = c(1, 10), criterion = 0.5)
#' 
#' @export
calculate_lambda_threshold <- function(max_distance = NULL, criterion = 0.7, stimulus_range = NULL) {
  validate_scalar(criterion, "criterion", lower = 0, upper = 1,
                  lower_open = TRUE, upper_open = TRUE)
  if (!is.null(stimulus_range)) {
    validate_rating_scale(stimulus_range)
    max_distance <- diff(stimulus_range)
  }
  if (is.null(max_distance)) stop("Either max_distance or stimulus_range must be provided")
  validate_scalar(max_distance, "max_distance", lower = 0, lower_open = TRUE)
  -log(criterion) / max_distance
}

#' Predict Group Membership
#' 
#' Predicts group membership for participants based on model fit
#' 
#' @param object An mfeargen object from fit_mfeargen
#' @param threshold Numeric. Minimum posterior probability for assignment.
#'   Participants below this threshold receive NA group labels (default: 0.5).
#' @param ... Additional arguments (ignored)
#' 
#' @return Data frame with participant IDs and predicted group assignments
#' 
#' @examples
#' \dontrun{
#' # Fit model
#' fit <- fit_mfeargen(data)
#' 
#' # Predict groups
#' groups <- predict_groups(fit)
#' }
#' 
#' @export
predict_groups <- function(object, threshold = 0.5, ...) {
  if (!inherits(object, "mfeargen")) stop("object must be of class 'mfeargen'")
  validate_scalar(threshold, "threshold", lower = 0, upper = 1)
  membership <- object$group_membership
  if (is.null(membership)) {
    membership <- extract_group_membership(object$jags_fit)
    ids <- attr(object$jags_data, "participants")
    if (!is.null(ids)) membership$participant <- ids
  }
  probability_columns <- paste0("prob_group_", 1:4)
  if (!all(probability_columns %in% names(membership))) {
    stop("Posterior group probabilities are unavailable; refit the model")
  }
  probabilities <- as.matrix(membership[probability_columns])
  assignments <- max.col(probabilities, ties.method = "first")
  max_probability <- apply(probabilities, 1, max)
  labels <- c("Non-Learner", "Overgeneralizer", "Physical", "Perceptual")
  membership$group <- ifelse(max_probability >= threshold, assignments, NA_integer_)
  membership$group_label <- labels[membership$group]
  membership$assigned_group <- membership$group_label
  membership$max_probability <- max_probability
  membership
}
