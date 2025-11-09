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
#' @return Logical. TRUE if valid, FALSE otherwise
#' @export
#' @examples
#' data <- simulate_mfeargen(n_participants = 10)
#' validate_mfeargen_data(data)
validate_mfeargen_data <- function(data, model_type = NULL, verbose = TRUE) {
  
  if (!is.data.frame(data)) {
    stop("Input must be a data frame")
  }
  
  # Required columns
  required_cols <- c("participant", "trials", "US_expect")
  
  # Check for required columns
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check data types
  if (!is.numeric(data$participant)) {
    stop("'participant' must be positive integers")
  }
  
  if (!is.numeric(data$trials)) {
    stop("'trials' must be positive integers")
  }
  
  if (!is.numeric(data$US_expect)) {
    stop("'US_expect' must be numeric")
  }
  
  # Check for model-specific columns
  if (!is.null(model_type)) {
    if (model_type == "1cs") {
      if (!"US" %in% names(data)) {
        stop("1CS model requires 'US' column")
      }
    } else if (model_type == "2cs") {
      if (!all(c("USp", "USm") %in% names(data))) {
        stop("2CS model requires 'USp' and 'USm' columns")
      }
    }
  }
  
  # Check for distance columns
  distance_cols <- c("Per_size", "Phy_size")
  missing_distance <- setdiff(distance_cols, names(data))
  if (length(missing_distance) > 0) {
    stop(
      "Missing required columns: ", paste(missing_distance, collapse = ", "),
      ". Provide perceptual and physical distance measures before fitting."
    )
  }
  
  # Check for valid values
  if (any(is.na(data$US_expect))) {
    if (verbose) warning("Data contains missing values in US_expect")
  }
  
  # Check for extreme values
  if (any(data$US_expect < 0 | data$US_expect > 100, na.rm = TRUE)) {
    if (verbose) warning("Some US_expect values outside typical rating scale range")
  }
  
  if (verbose) message("Data validation passed")
  return(TRUE)
}


#' Prepare Data for JAGS Model
#'
#' @description
#' Prepares data in the format required by JAGS for model fitting.
#' Handles both 1CS and 2CS models.
#'
#' @param data Data frame with fear generalization data
#' @param model_type Character. "1cs" or "2cs"
#' @param lambda_criterion Numeric. Criterion for lambda threshold (default: 0.7)
#' @param rating_scale Numeric vector. Rating scale range (optional)
#' @param priors mfeargen_priors object from set_priors() (optional)
#' @param use_flexible_priors Logical. Use flexible prior model (default: FALSE)
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
  
  # Handle missing US columns before validation
  if (model_type == "1cs" && !"US" %in% names(data)) {
    # Create US column from stimulus_type if available
    if ("stimulus_type" %in% names(data)) {
      data$US <- as.integer(data$stimulus_type == "CS+")
    } else {
      # Default to zeros
      data$US <- 0
    }
  } else if (model_type == "2cs") {
    if (!"USp" %in% names(data) || !"USm" %in% names(data)) {
      # Create USp and USm from stimulus_type if available
      if ("stimulus_type" %in% names(data)) {
        data$USp <- as.integer(data$stimulus_type == "CS+")
        data$USm <- as.integer(data$stimulus_type == "CS-")
      } else {
        # Default to zeros
        data$USp <- 0
        data$USm <- 0
      }
    }
  }
  
  # Validate structure (ensures required distance columns exist)
  validate_mfeargen_data(data, model_type = model_type, verbose = FALSE)
  
  # Validate data
  if (!validate_mfeargen_data(data, model_type, verbose = FALSE)) {
    stop("Data validation failed. Run validate_mfeargen_data() for details.")
  }
  
  # Get unique participants and trials
  participants <- unique(data$participant)
  n_participants <- length(participants)
  n_trials <- max(data$trials)
  
  # Initialize matrices
  y <- matrix(NA, n_participants, n_trials)
  
  # Fill response matrix
  for (i in seq_along(participants)) {
    p_data <- subset(data, participant == participants[i])
    p_data <- p_data[order(p_data$trials), ]
    y[i, ] <- p_data$US_expect
  }
  
  # Detect rating scale if not provided
  if (is.null(rating_scale)) {
    rating_scale <- detect_rating_scale(data$US_expect)
    message("Detected rating scale: [", rating_scale[1], ", ", rating_scale[2], "]")
  }
  
  # Calculate distances
  d_phy <- matrix(NA, n_participants, n_trials)
  d_per <- matrix(NA, n_participants, n_trials)
  
  for (i in seq_along(participants)) {
    p_data <- subset(data, participant == participants[i])
    p_data <- p_data[order(p_data$trials), ]
    
    # Calculate physical distances from CS+
    if ("stimulus_type" %in% names(p_data) && "Phy_size" %in% names(p_data)) {
      cs_plus <- median(p_data$Phy_size[p_data$stimulus_type == "CS+"], na.rm = TRUE)
      if (is.na(cs_plus)) cs_plus <- median(p_data$Phy_size, na.rm = TRUE)
    } else {
      cs_plus <- median(p_data$Phy_size, na.rm = TRUE)
    }
    
    if (!is.finite(cs_plus)) {
      stop("Unable to determine reference stimulus for distance calculation; check 'Phy_size' values.")
    }
    
    phy_vals <- abs(p_data$Phy_size - cs_plus)
    if (anyNA(phy_vals)) {
      phy_vals[is.na(phy_vals)] <- 0
    }
    per_vals <- abs(p_data$Per_size - cs_plus)
    if (anyNA(per_vals)) {
      fallback <- abs(p_data$Phy_size - cs_plus)
      per_vals[is.na(per_vals) & !is.na(fallback)] <- fallback[is.na(per_vals) & !is.na(fallback)]
      per_vals[is.na(per_vals)] <- 0
    }
    d_phy[i, ] <- phy_vals
    d_per[i, ] <- per_vals
  }
  
  # Calculate lambda threshold
  max_distance <- max(d_phy, na.rm = TRUE)
  if (max_distance == 0) max_distance <- 10
  lambda_threshold <- -log(lambda_criterion) / max_distance
  
  # Prepare model-specific data
  jags_data <- list(
    y = y,
    Nparticipants = n_participants,
    Ntrials = n_trials,
    lambda_threshold = lambda_threshold,
    rating_scale = rating_scale
  )
  
  if (model_type == "1cs") {
    # 1CS model data
    k <- r <- matrix(0, n_participants, n_trials)
    
    for (i in seq_along(participants)) {
      p_data <- subset(data, participant == participants[i])
      p_data <- p_data[order(p_data$trials), ]
      
      if ("US" %in% names(p_data)) {
        r[i, ] <- p_data$US
        if ("stimulus_type" %in% names(p_data)) {
          k[i, ] <- as.integer(p_data$stimulus_type == "CS+")
        } else {
          # Default: assume CS+ trials are when US is presented
          k[i, ] <- as.integer(p_data$US > 0)
        }
      }
    }
    
    jags_data$k <- k
    jags_data$r <- r
    jags_data$d_phy <- d_phy
    jags_data$d_per_old <- d_per
    
  } else {
    # 2CS model data
    kplus <- kminus <- rplus <- rminus <- matrix(0, n_participants, n_trials)
    d_p_phy <- d_m_phy <- matrix(NA, n_participants, n_trials)
    
    for (i in seq_along(participants)) {
      p_data <- subset(data, participant == participants[i])
      p_data <- p_data[order(p_data$trials), ]
      
      # Handle stimulus type indicators
      if ("stimulus_type" %in% names(p_data)) {
        kplus[i, ] <- as.integer(p_data$stimulus_type == "CS+")
        kminus[i, ] <- as.integer(p_data$stimulus_type == "CS-")
      } else if (all(c("USp", "USm") %in% names(p_data))) {
        # Default: use US indicators for stimulus type
        kplus[i, ] <- as.integer(p_data$USp > 0)
        kminus[i, ] <- as.integer(p_data$USm > 0)
      }
      
      # Handle reinforcement - CS+ can be reinforced, CS- is never reinforced
      if ("USp" %in% names(p_data)) {
        rplus[i, ] <- p_data$USp
      }
      
      # CS- is never reinforced (r_minus should always be 0 for CS- trials)
      # Even if USm column exists, CS- trials should have r_minus = 0
      rminus[i, ] <- rep(0, n_trials)  # CS- is never reinforced
      
      # Calculate distances from CS+ and CS-
      if ("stimulus_type" %in% names(p_data) && "Phy_size" %in% names(p_data)) {
        cs_plus <- median(p_data$Phy_size[p_data$stimulus_type == "CS+"], na.rm = TRUE)
        cs_minus <- median(p_data$Phy_size[p_data$stimulus_type == "CS-"], na.rm = TRUE)
      } else {
        # Default values if columns missing
        cs_plus <- 5.5
        cs_minus <- 3.5
      }
      
      if (!is.na(cs_plus) && "Phy_size" %in% names(p_data)) {
        d_p_phy[i, ] <- abs(p_data$Phy_size - cs_plus)
      } else {
        # Default distances
        d_p_phy[i, ] <- seq(0, 10, length.out = nrow(p_data))
      }
      
      if (!is.na(cs_minus) && "Phy_size" %in% names(p_data)) {
        d_m_phy[i, ] <- abs(p_data$Phy_size - cs_minus)
      } else {
        # Default distances
        d_m_phy[i, ] <- seq(0, 10, length.out = nrow(p_data))
      }
    }
    
    jags_data$k_plus <- kplus
    jags_data$k_minus <- kminus
    jags_data$r_plus <- rplus
    jags_data$r_minus <- rminus
    jags_data$d_p_phy <- d_p_phy
    jags_data$d_m_phy <- d_m_phy
  }
  
  # Add prior parameters if using flexible priors
  if (use_flexible_priors) {
    # Use provided priors or defaults
    if (is.null(priors)) {
      if (!requireNamespace("mfeargen", quietly = TRUE)) {
        stop("mfeargen package required for default priors")
      }
      priors <- get_default_priors()
    }
    
    # Validate priors
    if (!inherits(priors, "mfeargen_priors")) {
      stop("Priors must be created with set_priors() function")
    }
    
    # Add prior parameters to JAGS data
    jags_data$prior_lambda_mu_mean <- priors$lambda_mu$mean
    jags_data$prior_lambda_mu_precision <- priors$lambda_mu$precision
    jags_data$prior_lambda_sigma_min <- priors$lambda_sigma$min
    jags_data$prior_lambda_sigma_max <- priors$lambda_sigma$max
    
    jags_data$prior_alpha_mu_a <- priors$alpha_mu$a
    jags_data$prior_alpha_mu_b <- priors$alpha_mu$b
    jags_data$prior_alpha_kappa_min <- priors$alpha_kappa$min
    jags_data$prior_alpha_kappa_max <- priors$alpha_kappa$max
    
    jags_data$prior_sigma1_min <- priors$sigma1$min
    jags_data$prior_sigma1_max <- priors$sigma1$max
    jags_data$prior_sigma2_min <- priors$sigma2$min
    jags_data$prior_sigma2_max <- priors$sigma2$max
    
    jags_data$prior_pi_alpha <- priors$pi$alpha
    
    jags_data$prior_w0_mu_mean <- priors$w0_mu$mean
    jags_data$prior_w0_mu_precision <- priors$w0_mu$precision
    jags_data$prior_w0_sigma_shape <- priors$w0_sigma$shape
    jags_data$prior_w0_sigma_rate <- priors$w0_sigma$rate
    
    jags_data$prior_w1_a_shape <- priors$w1_a$shape
    jags_data$prior_w1_a_rate <- priors$w1_a$rate
    jags_data$prior_w1_b_shape <- priors$w1_b$shape
    jags_data$prior_w1_b_rate <- priors$w1_b$rate
  }
  
  return(jags_data)
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
  
  responses <- responses[!is.na(responses)]
  
  if (length(responses) == 0) {
    warning("No valid responses to detect scale")
    return(c(1, 10))
  }
  
  min_val <- min(responses)
  max_val <- max(responses)
  
  # Return the actual observed range of the data
  # This is the correct approach - use what participants actually used
  return(c(min_val - buffer, max_val + buffer))
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
  
  if (is.null(scale)) {
    scale <- detect_rating_scale(ratings)
  }
  
  if (length(scale) != 2) {
    stop("Scale must be a vector of length 2: c(min, max)")
  }
  
  normalized <- (ratings - scale[1]) / (scale[2] - scale[1])
  
  # Clip to [0, 1]
  normalized <- pmax(0, pmin(1, normalized))
  
  return(normalized)
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
  
  if (length(scale) != 2) {
    stop("Scale must be a vector of length 2: c(min, max)")
  }
  
  denormalized <- normalized * (scale[2] - scale[1]) + scale[1]
  
  return(denormalized)
}
