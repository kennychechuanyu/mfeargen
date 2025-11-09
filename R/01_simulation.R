#' @title Simulation Functions for mfeargen
#' @description Functions for simulating fear generalization data
#' @name simulation
#' @importFrom stats rnorm runif rgamma rbeta
NULL

#' Simulate Fear Generalization Data
#'
#' @description
#' Generates simulated fear generalization data following the two-phase paradigm:
#' (1) Conditioning phase with reinforcement learning
#' (2) Generalization test phase without reinforcement
#'
#' @param n_participants Integer. Number of participants (default: 40)
#' @param stimuli Numeric vector. All stimulus values (default: 1:10)
#' @param cs_plus Numeric. CS+ stimulus value (default: 5)
#' @param cs_minus Numeric. CS- stimulus value for 2cs (default: 3)
#' @param n_conditioning_trials Integer. CS+ trials in conditioning (default: 24)
#' @param n_generalization_trials Integer. Trials per stimulus in test (default: 16)
#' @param reinforcement_rate Numeric. CS+ reinforcement rate (default: 0.75)
#' @param model_type Character. "1cs" or "2cs" (default: "1cs")
#' @param group_proportions Numeric vector. Proportions for 4 groups (default: equal)
#' @param rating_scale Numeric vector. c(min, max) of scale (default: c(1, 10))
#' @param perceptual_noise_sd Numeric. SD of perceptual noise (default: 0.5)
#' @param lambda_criterion Numeric. Criterion for lambda threshold (default: 0.7)
#' @param seed Integer. Random seed for reproducibility
#' @param return_true_params Logical. Return true parameters (default: FALSE)
#'
#' @return Data frame with simulated data
#' @export
#' @examples
#' # Basic simulation with default stimuli (1:10), CS+ at 5
#' data <- simulate_mfeargen(n_participants = 20)
#'
#' # Custom stimuli matching paper design
#' data <- simulate_mfeargen(
#'   n_participants = 40,
#'   stimuli = c(73.67, 81.3, 88.92, 96.54, 104.17, 111.79, 119.42),
#'   cs_plus = 96.54,
#'   n_conditioning_trials = 24,
#'   n_generalization_trials = 16,
#'   rating_scale = c(0, 100)
#' )
#'
#' # 2CS design
#' data <- simulate_mfeargen(
#'   model_type = "2cs",
#'   stimuli = 1:10,
#'   cs_plus = 7,
#'   cs_minus = 3
#' )
simulate_mfeargen <- function(
  n_participants = 40,
  stimuli = 1:10,
  cs_plus = 5,
  cs_minus = 3,
  n_conditioning_trials = 24,
  n_generalization_trials = 16,
  reinforcement_rate = 0.75,
  model_type = c("1cs", "2cs"),
  group_proportions = c(0.25, 0.25, 0.25, 0.25),
  rating_scale = c(1, 10),
  perceptual_noise_sd = 0.5,
  lambda_criterion = 0.7,
  seed = NULL,
  return_true_params = FALSE
) {

  model_type <- match.arg(model_type)

  # Validate inputs
  if (n_participants <= 0) stop("n_participants must be positive")
  if (n_conditioning_trials <= 0) stop("n_conditioning_trials must be positive")
  if (n_generalization_trials < 0) stop("n_generalization_trials must be non-negative")
  if (abs(sum(group_proportions) - 1) > 1e-10) stop("group_proportions must sum to 1")

  # Validate CS+ is in stimuli
  if (!cs_plus %in% stimuli) {
    stop("cs_plus (", cs_plus, ") must be one of the stimuli: ", paste(stimuli, collapse = ", "))
  }

  # Validate CS- for 2cs
  if (model_type == "2cs" && !cs_minus %in% stimuli) {
    stop("cs_minus (", cs_minus, ") must be one of the stimuli: ", paste(stimuli, collapse = ", "))
  }

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Calculate lambda threshold based on stimulus range
  stimulus_range <- range(stimuli)
  max_distance <- diff(stimulus_range)
  lambda_threshold <- -log(lambda_criterion) / max_distance

  # Initialize storage
  all_data <- list()
  true_params <- list()
  participant_params <- list()

  # Generate data for each participant
  for (p in 1:n_participants) {

    # Assign to latent group
    group <- sample(1:4, 1, prob = group_proportions)

    # Generate group-specific parameters
    params <- generate_group_parameters(group, lambda_threshold, perceptual_noise_sd)

    # Store parameters
    participant_params[[p]] <- list(
      participant = p,
      alpha = params$alpha,
      lambda = params$lambda,
      sigma = params$sigma
    )

    if (return_true_params) {
      true_params[[p]] <- data.frame(
        participant = p,
        group = group,
        alpha_true = params$alpha,
        lambda_true = params$lambda,
        w0_true = params$w0,
        w1_true = params$w1,
        sigma_true = params$sigma
      )
    }

    # Generate participant data
    participant_data <- simulate_participant(
      p = p,
      group = group,
      params = params,
      stimuli = stimuli,
      cs_plus = cs_plus,
      cs_minus = cs_minus,
      n_conditioning_trials = n_conditioning_trials,
      n_generalization_trials = n_generalization_trials,
      reinforcement_rate = reinforcement_rate,
      model_type = model_type,
      rating_scale = rating_scale,
      perceptual_noise_sd = perceptual_noise_sd
    )

    all_data[[p]] <- participant_data
  }

  # Combine all participants
  final_data <- do.call(rbind, all_data)

  # Add group labels
  group_labels <- c("Non-Learner", "Overgeneralizer", "Physical", "Perceptual")
  final_data$group_label <- factor(group_labels[final_data$group], levels = group_labels)

  # Sort and clean
  final_data <- final_data[order(final_data$participant, final_data$trials), ]
  rownames(final_data) <- NULL

  # Process participant parameters
  param_summary <- if (length(participant_params) > 0) {
    list(
      alpha = sapply(participant_params, function(x) x$alpha),
      lambda = sapply(participant_params, function(x) x$lambda),
      sigma = sapply(participant_params, function(x) x$sigma)
    )
  } else NULL

  # Add metadata
  attr(final_data, "simulation_info") <- list(
    model_type = model_type,
    n_participants = n_participants,
    stimuli = stimuli,
    cs_plus = cs_plus,
    cs_minus = if (model_type == "2cs") cs_minus else NULL,
    n_conditioning_trials = n_conditioning_trials,
    n_generalization_trials = n_generalization_trials,
    group_assignment = tapply(final_data$true_group, final_data$participant, function(x) x[1]),
    parameters = param_summary,
    reinforcement_rate = reinforcement_rate,
    group_proportions = group_proportions,
    lambda_threshold = lambda_threshold,
    seed = seed
  )

  if (return_true_params) {
    attr(final_data, "true_parameters") <- do.call(rbind, true_params)
  }

  return(final_data)
}


#' Simulate Individual Participant
#'
#' @description Internal function to simulate one participant's data
#' @keywords internal
simulate_participant <- function(p, group, params, stimuli, cs_plus, cs_minus,
                                 n_conditioning_trials, n_generalization_trials,
                                 reinforcement_rate, model_type, rating_scale,
                                 perceptual_noise_sd) {

  # Extract parameters
  alpha <- params$alpha
  lambda <- params$lambda
  w0 <- params$w0
  w1 <- params$w1
  sigma <- params$sigma

  # Initialize associative strengths
  V_cs_plus <- 0
  V_cs_minus <- 0
  overall_trial <- 0

  # ========== PHASE 1: CONDITIONING ==========
  conditioning_list <- vector("list", n_conditioning_trials)

  for (trial in 1:n_conditioning_trials) {
    overall_trial <- overall_trial + 1

    # Determine stimulus and reinforcement
    if (model_type == "1cs") {
      stimulus <- cs_plus
      stimulus_type <- "CS+"
      US_presented <- as.integer(runif(1) < reinforcement_rate)
      g_value <- V_cs_plus
    } else {
      # 2CS: alternate between CS+ and CS-
      if (trial %% 2 == 1) {
        stimulus <- cs_plus
        stimulus_type <- "CS+"
        US_presented <- as.integer(runif(1) < reinforcement_rate)
        g_value <- V_cs_plus
      } else {
        stimulus <- cs_minus
        stimulus_type <- "CS-"
        US_presented <- 0
        g_value <- V_cs_minus
      }
    }

    # Generate perceptual size with noise
    Per_size <- stimulus + rnorm(1, 0, perceptual_noise_sd)
    stimulus_range <- range(stimuli)
    Per_size <- pmax(stimulus_range[1], pmin(stimulus_range[2], Per_size))
    Phy_size <- stimulus

    # Response generation
    theta <- generate_response(g_value, w0, w1, rating_scale)
    US_expect <- rnorm(1, theta, sigma)
    US_expect <- pmax(rating_scale[1], pmin(rating_scale[2], US_expect))

    # Store trial
    trial_data <- data.frame(
      participant = p,
      trials = as.integer(overall_trial),
      phase = "conditioning",
      stimulus = stimulus,
      stimulus_type = stimulus_type,
      US_expect = US_expect,
      Per_size = Per_size,
      Phy_size = Phy_size,
      group = group,
      true_group = as.integer(group),
      true_alpha = alpha,
      true_lambda = lambda,
      stringsAsFactors = FALSE
    )

    # Add US columns based on model type
    if (model_type == "2cs") {
      trial_data$USp <- ifelse(stimulus_type == "CS+", US_presented, 0)
      trial_data$USm <- ifelse(stimulus_type == "CS-", US_presented, 0)
    } else {
      trial_data$US <- US_presented
    }

    conditioning_list[[trial]] <- trial_data

    # Learning update (only for CS+/CS-, not GS)
    if (group != 1) {
      if (stimulus_type == "CS+") {
        V_cs_plus <- V_cs_plus + alpha * (US_presented - V_cs_plus)
      } else if (stimulus_type == "CS-") {
        V_cs_minus <- V_cs_minus + alpha * (US_presented - V_cs_minus)
      }
    }
  }
  conditioning_data <- do.call(rbind, conditioning_list)

  # ========== PHASE 2: GENERALIZATION ==========
  generalization_list <- list()

  if (n_generalization_trials > 0) {
    # Create test sequence: repeat all stimuli n_generalization_trials times
    test_sequence <- rep(stimuli, each = n_generalization_trials)
    test_sequence <- sample(test_sequence)  # Randomize order

    for (trial in seq_along(test_sequence)) {
      overall_trial <- overall_trial + 1
      stimulus <- test_sequence[trial]

      # Determine stimulus type
      if (abs(stimulus - cs_plus) < 1e-9) {
        stimulus_type <- "CS+"
      } else if (model_type == "2cs" && abs(stimulus - cs_minus) < 1e-9) {
        stimulus_type <- "CS-"
      } else {
        stimulus_type <- "GS"
      }

      # No reinforcement during test
      US_presented <- 0

      # Generate perceptual size with noise
      Per_size <- stimulus + rnorm(1, 0, perceptual_noise_sd)
      stimulus_range <- range(stimuli)
      Per_size <- pmax(stimulus_range[1], pmin(stimulus_range[2], Per_size))
      Phy_size <- stimulus

      # Calculate generalization
      if (group == 1) {
        g_value <- 0  # Non-learners
      } else {
        # Distance calculation based on group
        if (group == 4) {
          # Perceptual generalizers use perceived distance
          distance_plus <- abs(Per_size - cs_plus)
        } else {
          # Physical generalizers use physical distance
          distance_plus <- abs(Phy_size - cs_plus)
        }

        similarity_plus <- exp(-lambda * distance_plus)
        g_value <- V_cs_plus * similarity_plus

        # Add CS- contribution for 2CS
        if (model_type == "2cs") {
          if (group == 4) {
            distance_minus <- abs(Per_size - cs_minus)
          } else {
            distance_minus <- abs(Phy_size - cs_minus)
          }
          similarity_minus <- exp(-lambda * distance_minus)
          g_value <- g_value + V_cs_minus * similarity_minus
        }
      }

      # Response generation
      theta <- generate_response(g_value, w0, w1, rating_scale)
      US_expect <- rnorm(1, theta, sigma)
      US_expect <- pmax(rating_scale[1], pmin(rating_scale[2], US_expect))

      # Store trial
      trial_data <- data.frame(
        participant = p,
        trials = as.integer(overall_trial),
        phase = "generalization",
        stimulus = stimulus,
        stimulus_type = stimulus_type,
        US_expect = US_expect,
        Per_size = Per_size,
        Phy_size = Phy_size,
        group = group,
        true_group = as.integer(group),
        true_alpha = alpha,
        true_lambda = lambda,
        stringsAsFactors = FALSE
      )

      # Add US columns (all zeros in test phase)
      if (model_type == "2cs") {
        trial_data$USp <- 0
        trial_data$USm <- 0
      } else {
        trial_data$US <- 0
      }
      generalization_list[[trial]] <- trial_data
    }
    generalization_data <- if(length(generalization_list) > 0) do.call(rbind, generalization_list) else data.frame()
  } else {
    generalization_data <- data.frame()
  }

  return(rbind(conditioning_data, generalization_data))
}


#' Generate Group-Specific Parameters
#'
#' @param group Integer (1-4) indicating latent group
#' @param lambda_threshold Numeric. Threshold for lambda values
#' @param perceptual_noise_sd Numeric. SD for perceptual noise
#' @return List with parameters
#' @keywords internal
generate_group_parameters <- function(group, lambda_threshold, perceptual_noise_sd = 0.5) {

  if (group == 1) {
    # Non-learners: No learning or generalization
    list(
      alpha = 0,
      lambda = 0,
      w0 = rnorm(1, 0, 2),
      w1 = 0,
      sigma = runif(1, 2, 3)
    )
  } else if (group == 2) {
    # Overgeneralizers: Good learning, broad generalization
    list(
      alpha = rbeta(1, 3, 2),
      lambda = runif(1, 0, lambda_threshold * 0.5),
      w0 = rnorm(1, -1, 0.5),
      w1 = rgamma(1, 3, 1) + 1,
      sigma = runif(1, 0.5, 1.5)
    )
  } else if (group == 3) {
    # Physical generalizers: Sharp generalization based on physical distance
    lambda_val <- rnorm(1, lambda_threshold * 2, lambda_threshold * 0.3)
    list(
      alpha = rbeta(1, 3, 2),
      lambda = pmax(lambda_threshold, lambda_val),
      w0 = rnorm(1, -1, 0.5),
      w1 = rgamma(1, 3, 1) + 1,
      sigma = runif(1, 0.4, 1.2)
    )
  } else {
    # Perceptual generalizers: Sharp generalization based on perceptual distance
    lambda_val <- rnorm(1, lambda_threshold * 2, lambda_threshold * 0.3)
    list(
      alpha = rbeta(1, 3, 2),
      lambda = pmax(lambda_threshold, lambda_val),
      w0 = rnorm(1, -1, 0.5),
      w1 = rgamma(1, 3, 1) + 1,
      sigma = runif(1, 0.4, 1.2)
    )
  }
}


#' Generate Response Using Sigmoid Transformation
#'
#' @param g_value Numeric. Generalized associative strength
#' @param w0 Numeric. Intercept parameter
#' @param w1 Numeric. Slope parameter
#' @param rating_scale Numeric vector. c(min, max) of scale
#' @return Numeric. Expected response value
#' @keywords internal
generate_response <- function(g_value, w0, w1, rating_scale) {
  A <- rating_scale[1]
  K <- rating_scale[2]
  theta <- A + (K - A) / (1 + exp(-(w0 + w1 * g_value)))
  return(theta)
}
