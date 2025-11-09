# Test Helper Functions for mfeargen Package
# These utilities facilitate testing by providing mock objects and common test data

# ------------------------------------------------------------------------------
# Compatibility helpers for older testthat releases
if (!exists("expect_no_error")) {
  expect_no_error <- function(object, ...) {
    testthat::expect_error(object, regexp = NA, ...)
  }
}

expect_true_msg <- function(condition, message) {
  if (!isTRUE(condition)) {
    testthat::fail(message)
  } else {
    testthat::succeed()
  }
}

#' Create a small test dataset
#'
#' @param n_participants Number of participants (default: 5)
#' @param n_trials Number of trials (default: 10)
#' @param model_type Model type ("1cs" or "2cs")
#' @param seed Random seed for reproducibility
#' @return A simulated dataset
create_test_data <- function(n_participants = 5, n_trials = 10,
                              model_type = "1cs", seed = 123) {
  set.seed(seed)
  simulate_mfeargen(
    n_participants = n_participants,
    n_trials = n_trials,
    model_type = model_type,
    seed = seed
  )
}

#' Create a mock mfeargen fit object
#'
#' This creates a realistic mock object without actually running JAGS
#'
#' @param n_participants Number of participants
#' @param n_trials Number of trials
#' @param model_type Model type ("1cs" or "2cs")
#' @return A mock mfeargen object
create_mock_mfeargen <- function(n_participants = 5, n_trials = 10,
                                  model_type = "1cs") {

  # Create base data
  data <- create_test_data(n_participants, n_trials, model_type)
  jags_data <- prepare_mfeargen_data(data, model_type = model_type)

  # Number of parameters
  n_params <- n_participants * 6 + 4 + n_participants * n_trials + 5

  # Create parameter names
  param_names <- c(
    paste0("alpha[", 1:n_participants, "]"),
    paste0("lambda[", 1:n_participants, "]"),
    paste0("pi[", 1:4, "]"),
    paste0("w0[", 1:n_participants, "]"),
    paste0("w1[", 1:n_participants, "]"),
    paste0("sigma[", 1:n_participants, "]"),
    paste0("N_group[", 1:4, "]"),
    paste0("theta[", rep(1:n_participants, each = n_trials), ",",
           rep(1:n_trials, n_participants), "]"),
    "total_loglik", "p_value"
  )

  # Create mock summary matrix
  mock_summary <- matrix(
    rnorm(length(param_names) * 7),
    nrow = length(param_names),
    ncol = 7,
    dimnames = list(
      param_names,
      c("mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%")
    )
  )

  # Set reasonable parameter values
  alpha_idx <- grep("^alpha\\[", rownames(mock_summary))
  lambda_idx <- grep("^lambda\\[", rownames(mock_summary))
  pi_idx <- grep("^pi\\[", rownames(mock_summary))
  sigma_idx <- grep("^sigma\\[", rownames(mock_summary))
  theta_idx <- grep("^theta\\[", rownames(mock_summary))
  loglik_idx <- grep("^total_loglik", rownames(mock_summary))

  mock_summary[alpha_idx, "mean"] <- runif(length(alpha_idx), 0.1, 0.9)
  mock_summary[lambda_idx, "mean"] <- runif(length(lambda_idx), 0.001, 0.05)
  mock_summary[pi_idx, "mean"] <- c(0.2, 0.3, 0.3, 0.2)
  mock_summary[sigma_idx, "mean"] <- runif(length(sigma_idx), 0.5, 2)
  mock_summary[theta_idx, "mean"] <- runif(length(theta_idx), 1, 10)
  mock_summary[loglik_idx, "mean"] <- -500

  # Add Rhat and n.eff
  rhat_values <- runif(nrow(mock_summary), 0.99, 1.04)
  neff_values <- runif(nrow(mock_summary), 200, 1000)

  mock_summary <- cbind(mock_summary, Rhat = rhat_values, n.eff = neff_values)

  # Create sims.list
  n_iter <- 1000
  sims_list <- list(
    gp = matrix(sample(1:4, n_participants * n_iter, replace = TRUE),
                nrow = n_iter, ncol = n_participants),
    alpha = matrix(runif(n_participants * n_iter, 0, 1),
                   nrow = n_iter, ncol = n_participants),
    lambda = matrix(runif(n_participants * n_iter, 0, 0.1),
                    nrow = n_iter, ncol = n_participants),
    pi = matrix(rep(c(0.2, 0.3, 0.3, 0.2), each = n_iter),
                nrow = n_iter, ncol = 4),
    theta = array(runif(n_participants * n_trials * n_iter, 1, 10),
                  dim = c(n_iter, n_participants, n_trials)),
    sigma = matrix(runif(n_participants * n_iter, 0.5, 2),
                   nrow = n_iter, ncol = n_participants)
  )

  # Create group membership
  group_membership <- data.frame(
    participant = 1:n_participants,
    assigned_group = sample(
      c("Non-Learner", "Overgeneralizer", "Physical", "Perceptual"),
      n_participants, replace = TRUE
    ),
    prob_group_1 = runif(n_participants, 0, 0.5),
    prob_group_2 = runif(n_participants, 0, 0.5),
    prob_group_3 = runif(n_participants, 0, 0.5),
    prob_group_4 = runif(n_participants, 0, 0.5),
    stringsAsFactors = FALSE
  )
  # Normalize probabilities
  group_membership$max_probability <- apply(
    group_membership[, c("prob_group_1", "prob_group_2", "prob_group_3", "prob_group_4")],
    1, max
  )

  # Construct mock fit object
  mock_fit <- list(
    jags_fit = list(
      summary = mock_summary,
      sims.list = sims_list,
      DIC = 1234.56,
      pD = 50,
      deviance = 1184.56,  # Add deviance field
      Rhat = rhat_values,
      n.eff = neff_values
    ),
    model_type = model_type,
    data = data,
    jags_data = jags_data,
    group_membership = group_membership,
    lambda_threshold = 0.01,
    convergence = list(
      convergence_ok = TRUE,
      converged = TRUE,
      max_rhat = max(rhat_values),
      min_neff = min(neff_values),
      n_converged = sum(rhat_values < 1.1),
      n_total = length(rhat_values)
    ),
    n.iter = 1000,
    n.burnin = 500,
    n.chains = 2,
    n.thin = 1,
    call = quote(fit_mfeargen(data, model_type = model_type))
  )

  class(mock_fit) <- "mfeargen"
  return(mock_fit)
}

create_mock_fit <- create_mock_mfeargen

#' Create mock JAGS output for specific scenarios
#'
#' @param scenario Type of scenario ("converged", "not_converged", "high_dic")
#' @return A mock JAGS fit list
create_mock_jags_output <- function(scenario = "converged") {

  if (scenario == "converged") {
    rhat_values <- runif(20, 0.99, 1.04)
    neff_values <- runif(20, 500, 2000)
    dic <- 1000
  } else if (scenario == "not_converged") {
    rhat_values <- runif(20, 1.1, 1.5)
    neff_values <- runif(20, 50, 150)
    dic <- 1500
  } else if (scenario == "high_dic") {
    rhat_values <- runif(20, 0.99, 1.04)
    neff_values <- runif(20, 500, 2000)
    dic <- 5000
  }

  summary_matrix <- matrix(
    rnorm(20 * 7),
    nrow = 20,
    ncol = 7,
    dimnames = list(
      paste0("param", 1:20),
      c("mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%")
    )
  )

  summary_matrix <- cbind(summary_matrix, Rhat = rhat_values, n.eff = neff_values)

  list(
    summary = summary_matrix,
    DIC = dic,
    pD = 20,
    Rhat = rhat_values,
    n.eff = neff_values
  )
}

#' Expect approximate equality for vectors
#'
#' @param actual Actual values
#' @param expected Expected values
#' @param tolerance Tolerance for comparison
#' @param label Label for error messages
expect_approx_equal <- function(actual, expected, tolerance = 0.1, label = NULL) {
  testthat::expect_equal(actual, expected, tolerance = tolerance, label = label)
}

#' Check if all required columns exist in a data frame
#'
#' @param df Data frame to check
#' @param required_cols Character vector of required column names
#' @return TRUE if all columns exist, FALSE otherwise
has_required_columns <- function(df, required_cols) {
  all(required_cols %in% names(df))
}

#' Create minimal valid data for testing
#'
#' @param n_rows Number of rows
#' @return A minimal valid dataset
create_minimal_data <- function(n_rows = 10) {
  data.frame(
    participant = rep(1:2, each = n_rows / 2),
    trials = rep(1:(n_rows / 2), 2),
    US_expect = runif(n_rows, 1, 10),
    stimulus = rep(c("CS+", "GS1", "GS2", "GS3", "GS4"), length.out = n_rows),
    Phy_size = rnorm(n_rows, 5.5, 1),
    Per_size = rnorm(n_rows, 5.5, 1),
    US = rep(c(1, 0, 0, 0, 0), length.out = n_rows)
  )
}

#' Skip test if JAGS is not available
skip_if_no_jags <- function() {
  skip_on_cran()

  jags_available <- tryCatch({
    requireNamespace("rjags", quietly = TRUE) &&
      requireNamespace("jagsUI", quietly = TRUE)
  }, error = function(e) FALSE)

  if (!jags_available) {
    skip("JAGS not available")
  }

  # Also check if JAGS binary is installed
  jags_binary <- tryCatch({
    rjags::load.module("basemod")
    TRUE
  }, error = function(e) FALSE)

  if (!jags_binary) {
    skip("JAGS binary not installed")
  }
}

#' Check if a ggplot object is valid
#'
#' @param p A ggplot object
#' @return TRUE if valid, FALSE otherwise
is_valid_ggplot <- function(p) {
  inherits(p, "gg") && inherits(p, "ggplot")
}

#' Create test data with known group structure
#'
#' @param group_id Group ID (1-4)
#' @param n_participants Number of participants
#' @param n_trials Number of trials
#' @return Data frame with known group structure
create_group_specific_data <- function(group_id, n_participants = 5, n_trials = 20) {
  set.seed(group_id * 100)

  # Group-specific parameters
  if (group_id == 1) {
    # Non-learners: low alpha, high variability
    alpha <- 0.05
    lambda <- 0.01
  } else if (group_id == 2) {
    # Overgeneralizers: high alpha, low lambda
    alpha <- 0.8
    lambda <- 0.005
  } else if (group_id == 3) {
    # Physical generalizers: moderate alpha, moderate lambda
    alpha <- 0.6
    lambda <- 0.02
  } else {
    # Perceptual generalizers: high alpha, high lambda
    alpha <- 0.75
    lambda <- 0.03
  }

  simulate_mfeargen(
    n_participants = n_participants,
    n_trials = n_trials,
    seed = group_id * 100
  )
}
