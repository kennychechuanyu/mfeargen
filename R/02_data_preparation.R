#' @title Data Preparation Functions for mfeargen
#' @description Functions for preparing and validating data for model fitting
#' @name data_preparation
#' @importFrom stats median rnorm aggregate
NULL

#' Validate mfeargen Data
#'
#' @description
#' Validates that a data frame has the required structure for mfeargen analysis.
#' Checks for required columns, data types, and value ranges.
#'
#' @param data Data frame to validate
#' @param model_type Character. "1cs" or "2cs" (optional)
#' @param verbose Logical. Print validation messages (default: TRUE)
#' @return TRUE if valid; invalid input raises an error.
#' @export
#' @examples
#' data <- simulate_mfeargen(n_participants = 10)
#' validate_mfeargen_data(data)
validate_mfeargen_data <- function(data, model_type = NULL, verbose = TRUE) {
  if (!is.data.frame(data) || nrow(data) == 0L) {
    stop("data must be a non-empty data frame")
  }
  required <- c("participant", "trials", "US_expect")
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))
  if (!(is.numeric(data$participant) || is.character(data$participant) ||
        is.factor(data$participant)) || anyNA(data$participant) ||
      any(!nzchar(as.character(data$participant)))) {
    stop("participant must contain non-missing numeric or character identifiers")
  }
  if (is.numeric(data$participant) && any(!is.finite(data$participant))) {
    stop("participant identifiers must be finite")
  }
  if (!is.numeric(data$trials) || anyNA(data$trials) ||
      any(!is.finite(data$trials) | data$trials < 1 | data$trials != floor(data$trials))) {
    stop("trials must contain positive integers")
  }
  if (anyDuplicated(data[c("participant", "trials")])) {
    stop("Duplicate trial numbers within participant")
  }
  if (!is.numeric(data$US_expect) || any(is.infinite(data$US_expect))) {
    stop("US_expect must be numeric, with finite values or NA")
  }
  if (any(vapply(split(data$US_expect, data$participant, drop = TRUE), function(x) all(is.na(x)), logical(1)))) {
    stop("Each participant must have at least one observed US_expect rating")
  }
  if (is.null(model_type)) model_type <- if (all(c("USp", "USm") %in% names(data))) "2cs" else "1cs"
  if (is.null(model_type)) {
    model_type <- if (all(c("USp", "USm") %in% names(data))) "2cs" else "1cs"
  }
  if (!is.null(model_type)) {
    model_type <- match.arg(model_type, c("1cs", "2cs"))
    outcomes <- if (model_type == "1cs") "US" else c("USp", "USm")
    missing <- setdiff(outcomes, names(data))
    if (length(missing)) stop("Missing required reinforcement columns: ", paste(missing, collapse = ", "))
    for (column in outcomes) {
      x <- data[[column]]
      allowed <- if (column == "USm") c(-1, 0, 1) else c(0, 1)
      if (!is.numeric(x) || anyNA(x) || any(!x %in% allowed)) {
        stop(column, " must contain observed reinforcement outcomes coded ", paste(allowed, collapse = ", "))
      }
    }
  }
  indicators <- if (model_type == "1cs") "k" else c("k_plus", "k_minus")
  distances <- if (model_type == "1cs") c("d_phy", "d_per_old") else
    c("d_p_phy", "d_m_phy", "d_p_per", "d_m_per")
  if (all(distances %in% names(data))) {
    for (column in distances) {
      if (!is.numeric(data[[column]]) || any(!is.finite(data[[column]]) | data[[column]] < 0)) {
        stop(column, " must contain finite, non-negative distances")
      }
    }
  } else {
    missing <- setdiff(c("Phy_size", "Per_size", "stimulus_type"), names(data))
    if (length(missing)) stop("Missing required columns for stimulus distances: ", paste(missing, collapse = ", "))
    for (column in c("Phy_size", "Per_size")) {
      if (!is.numeric(data[[column]]) || any(!is.finite(data[[column]]))) {
        stop(column, " must contain finite measurements; missing covariates cannot be inferred")
      }
    }
  }
  if (all(indicators %in% names(data))) {
    for (column in indicators) {
      if (!is.numeric(data[[column]]) || anyNA(data[[column]]) || any(!data[[column]] %in% c(0, 1))) {
        stop(column, " must contain learning indicators coded 0 or 1")
      }
    }
  } else if (!"stimulus_type" %in% names(data) || anyNA(data$stimulus_type)) {
    stop("Provide non-missing stimulus_type labels or explicit learning indicators")
  }
  if ("phase" %in% names(data) &&
      (anyNA(data$phase) || any(!data$phase %in% c("conditioning", "generalization")))) {
    stop("phase must contain 'conditioning' or 'generalization'")
  }
  if (verbose) message("Data validation passed")
  TRUE
}

#' Prepare Data for JAGS Model
#'
#' @description
#' Prepares data in the format required by JAGS for model fitting.
#' Handles both 1CS and 2CS models. Responses are rescaled internally to 0-10;
#' the original bounds and participant/trial mapping are stored as attributes.
#' Unequal trial counts are padded after the last trial. Missing ratings are
#' allowed when reinforcement and stimulus covariates are observed.
#'
#' @param data Data frame with participant, trials and US_expect, plus observed
#'   reinforcement outcomes (US for 1cs; USp and signed USm for 2cs). USm=-1
#'   denotes an inhibitory target, 0 no target, and 1 an excitatory target.
#'   Supply stimulus_type, Phy_size and Per_size, or explicit trial indicators
#'   and distance columns. When phase is supplied, learning updates occur
#'   during conditioning only.
#' @param model_type Character. "1cs" or "2cs"
#' @param lambda_criterion Numeric. Criterion for lambda threshold (default: 0.7)
#' @param rating_scale Numeric vector c(min, max). Uses simulation metadata
#'   when available, otherwise the observed range. Specify the experimental
#'   scale explicitly when participants did not use its full range.
#' @param priors Must be NULL; supply priors to fit_mfeargen instead.
#' @param use_flexible_priors Must be FALSE; supply priors to fit_mfeargen instead.
#' @return List formatted for JAGS
#' @export
#' @examples
#' data <- simulate_mfeargen(n_participants = 10)
#' jags_data <- prepare_mfeargen_data(data, model_type = "1cs")
prepare_mfeargen_data <- function(data,
                                 model_type = c("1cs", "2cs"),
                                 lambda_criterion = 0.7,
                                 rating_scale = NULL,
                                 priors = NULL,
                                 use_flexible_priors = FALSE) {
  model_type <- match.arg(model_type)
  if (!is.null(priors) || !identical(use_flexible_priors, FALSE)) {
    stop("Supply priors to fit_mfeargen(), not prepare_mfeargen_data()")
  }
  validate_mfeargen_data(data, model_type, verbose = FALSE)
  validate_scalar(lambda_criterion, "lambda_criterion", lower = 0, upper = 1,
                  lower_open = TRUE, upper_open = TRUE)

  participants <- unique(data$participant)
  participant_rows <- lapply(participants, function(id) {
    rows <- which(data$participant == id)
    rows[order(data$trials[rows])]
  })
  counts <- lengths(participant_rows)
  n_participants <- length(participants)
  n_trials <- max(counts)
  row_index <- matrix(NA_integer_, n_participants, n_trials)
  for (i in seq_along(participants)) row_index[i, seq_len(counts[i])] <- participant_rows[[i]]
  to_matrix <- function(values, padding = 0) {
    result <- matrix(padding, n_participants, n_trials)
    present <- !is.na(row_index)
    result[present] <- values[row_index[present]]
    result
  }

  if (is.null(rating_scale)) rating_scale <- attr(data, "rating_scale")
  if (is.null(rating_scale)) rating_scale <- detect_rating_scale(data$US_expect)
  validate_rating_scale(rating_scale)
  if (any(data$US_expect < rating_scale[1] | data$US_expect > rating_scale[2], na.rm = TRUE)) {
    stop("US_expect contains ratings outside rating_scale")
  }

  indicators <- if (model_type == "1cs") "k" else c("k_plus", "k_minus")
  distances <- if (model_type == "1cs") c("d_phy", "d_per_old") else {
    c("d_p_phy", "d_m_phy", "d_p_per", "d_m_per")
  }
  explicit_indicators <- all(indicators %in% names(data))
  explicit_distances <- all(distances %in% names(data))
  if (!explicit_indicators) {
    if (!"stimulus_type" %in% names(data) || anyNA(data$stimulus_type) ||
        any(!data$stimulus_type %in% c("CS+", "CS-", "GS"))) {
      stop("Provide stimulus_type ('CS+', 'CS-', 'GS') or explicit learning indicators: ",
           paste(indicators, collapse = ", "))
    }
    data[[indicators[1]]] <- as.integer(data$stimulus_type == "CS+")
    if (model_type == "2cs") data$k_minus <- as.integer(data$stimulus_type == "CS-")
  }
  for (column in indicators) {
    if (!is.numeric(data[[column]]) || anyNA(data[[column]]) || any(!data[[column]] %in% c(0, 1))) {
      stop(column, " must contain learning indicators coded 0 or 1")
    }
  }
  if (model_type == "2cs" && any(data$k_plus + data$k_minus > 1)) {
    stop("k_plus and k_minus cannot both be 1 on the same trial")
  }

  if (!explicit_distances) {
    required <- c("Phy_size", "Per_size", "stimulus_type")
    missing <- setdiff(required, names(data))
    if (length(missing)) stop("Missing required columns for stimulus distances: ", paste(missing, collapse = ", "))
    for (column in c("Phy_size", "Per_size")) {
      if (!is.numeric(data[[column]]) || any(!is.finite(data[[column]]))) {
        stop(column, " must contain finite measurements; missing covariates cannot be inferred")
      }
    }
    reference <- function(rows, label) {
      values <- unique(data$Phy_size[rows][data$stimulus_type[rows] == label])
      if (length(values) != 1L) {
        stop("Each participant needs exactly one physical reference for ", label,
             "; provide explicit distance columns for other designs")
      }
      values
    }
    for (column in distances) data[[column]] <- numeric(nrow(data))
    for (rows in participant_rows) {
      cs_plus <- reference(rows, "CS+")
      if (model_type == "1cs") {
        data$d_phy[rows] <- abs(data$Phy_size[rows] - cs_plus)
        data$d_per_old[rows] <- abs(data$Per_size[rows] - cs_plus)
      } else {
        cs_minus <- reference(rows, "CS-")
        if (cs_plus == cs_minus) stop("CS+ and CS- references must differ")
        data$d_p_phy[rows] <- abs(data$Phy_size[rows] - cs_plus)
        data$d_m_phy[rows] <- abs(data$Phy_size[rows] - cs_minus)
        data$d_p_per[rows] <- abs(data$Per_size[rows] - cs_plus)
        data$d_m_per[rows] <- abs(data$Per_size[rows] - cs_minus)
      }
    }
  }
  for (column in distances) {
    if (!is.numeric(data[[column]]) || any(!is.finite(data[[column]]) | data[[column]] < 0)) {
      stop(column, " must contain finite, non-negative distances")
    }
  }
  if ("phase" %in% names(data)) {
    for (column in indicators) data[[column]] <- data[[column]] * as.integer(data$phase == "conditioning")
  }

  physical_columns <- if (model_type == "1cs") "d_phy" else c("d_p_phy", "d_m_phy")
  max_distance <- max(unlist(data[physical_columns], use.names = FALSE))
  if (max_distance <= 0) stop("At least one nonzero physical distance is required to define lambda_threshold")
  jags_data <- list(
    y = to_matrix(10 * normalize_ratings(data$US_expect, rating_scale), NA_real_),
    Nparticipants = n_participants,
    Ntrials = n_trials,
    lambda_threshold = calculate_lambda_threshold(max_distance, lambda_criterion),
    rating_scale = c(0, 10)
  )
  for (column in c(indicators, distances)) jags_data[[column]] <- to_matrix(data[[column]])
  if (model_type == "1cs") {
    jags_data$r <- to_matrix(data$US)
  } else {
    jags_data$r_plus <- to_matrix(data$USp)
    jags_data$r_minus <- to_matrix(data$USm)
  }
  attr(jags_data, "participants") <- participants
  attr(jags_data, "trial_counts") <- counts
  attr(jags_data, "row_index") <- row_index
  attr(jags_data, "observed") <- !is.na(jags_data$y)
  attr(jags_data, "model_type") <- model_type
  attr(jags_data, "rating_scale") <- rating_scale
  jags_data
}

#' Detect Rating Scale Range
#'
#' @description
#' Automatically detects the rating scale range from response data.
#'
#' @param responses Numeric vector of responses
#' @param buffer Numeric. Buffer to add to detected range (default: 0)
#' @return Numeric vector c(min, max)
#' @export
#' @examples
#' responses <- c(1, 2, 5, 8, 9, 10)
#' detect_rating_scale(responses)
detect_rating_scale <- function(responses, buffer = 0) {
  if (!is.numeric(responses) || any(is.infinite(responses))) {
    stop("responses must be numeric, with finite values or NA")
  }
  validate_scalar(buffer, "buffer", lower = 0)
  responses <- responses[!is.na(responses)]
  if (!length(responses)) stop("Cannot detect a rating scale without observed responses")
  scale <- range(responses) + c(-buffer, buffer)
  if (scale[1] == scale[2]) stop("Constant ratings require an explicit rating_scale")
  scale
}

#' Normalize Ratings to 0-1 Scale
#'
#' @description
#' Normalizes ratings to a 0-1 scale for model fitting.
#'
#' @param ratings Numeric vector of ratings
#' @param scale Numeric vector c(min, max) of original scale
#' @return Numeric vector of normalized ratings
#' @export
#' @examples
#' ratings <- c(1, 5, 10)
#' normalize_ratings(ratings, scale = c(1, 10))
normalize_ratings <- function(ratings, scale = NULL) {
  if (!is.numeric(ratings) || any(is.infinite(ratings))) stop("ratings must be numeric and finite or NA")
  if (is.null(scale)) scale <- detect_rating_scale(ratings)
  validate_rating_scale(scale)
  if (any(ratings < scale[1] | ratings > scale[2], na.rm = TRUE)) stop("ratings lie outside scale")
  (ratings - scale[1]) / diff(scale)
}

#' Denormalize Ratings from 0-1 Scale
#'
#' @description
#' Converts normalized ratings back to original scale.
#'
#' @param normalized Numeric vector of normalized ratings (0-1)
#' @param scale Numeric vector c(min, max) of target scale
#' @return Numeric vector of denormalized ratings
#' @export
#' @examples
#' normalized <- c(0, 0.5, 1)
#' denormalize_ratings(normalized, scale = c(1, 10))
denormalize_ratings <- function(normalized, scale) {
  validate_rating_scale(scale)
  if (!is.numeric(normalized) || any(is.infinite(normalized)) ||
      any(normalized < 0 | normalized > 1, na.rm = TRUE)) {
    stop("normalized must contain values between 0 and 1, or NA")
  }
  normalized * diff(scale) + scale[1]
}
