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
      any(!is.finite(scale)) || scale[1] >= scale[2] || !is.finite(diff(scale))) {
    stop("rating_scale must contain two finite, increasing bounds", call. = FALSE)
  }
  invisible(scale)
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
  if (!is.null(stimulus_range) && !is.null(max_distance)) {
    stop("Supply only one of max_distance or stimulus_range")
  }
  if (!is.null(stimulus_range)) {
    validate_rating_scale(stimulus_range)
    max_distance <- diff(stimulus_range)
  }
  if (is.null(max_distance)) stop("Either max_distance or stimulus_range must be provided")
  validate_scalar(max_distance, "max_distance", lower = 0, lower_open = TRUE)
  threshold <- -log(criterion) / max_distance
  if (!is.finite(threshold) || threshold <= 0) {
    stop("The lambda threshold cannot be represented numerically; rescale stimulus distances")
  }
  threshold
}

#' Predict Group Membership
#' 
#' Predicts group membership for participants based on model fit
#' 
#' @param object A fitted mfeargen or mfeargen_hierarchical model.
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
  if (!inherits(object, c("mfeargen", "mfeargen_hierarchical"))) {
    stop("object must be a fitted mfeargen or mfeargen_hierarchical model")
  }
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
  if (!is.numeric(probabilities) || !nrow(probabilities) ||
      any(!is.finite(probabilities)) || any(probabilities < 0 | probabilities > 1) ||
      any(abs(rowSums(probabilities) - 1) > 1e-8)) {
    stop("Posterior group probabilities must be finite, non-negative, and sum to one")
  }
  assignments <- max.col(probabilities, ties.method = "first")
  max_probability <- apply(probabilities, 1, max)
  labels <- c("Non-Learner", "Overgeneralizer", "Physical", "Perceptual")
  membership$group <- ifelse(max_probability >= threshold, assignments, NA_integer_)
  membership$group_label <- labels[membership$group]
  membership$assigned_group <- membership$group_label
  membership$max_probability <- max_probability
  membership
}
