create_test_data <- function(n_participants = 4, n_trials = 8,
                             model_type = "1cs", seed = 123, ...) {
  simulate_mfeargen(
    n_participants = n_participants,
    n_conditioning_trials = n_trials,
    n_generalization_trials = 1,
    model_type = model_type,
    seed = seed,
    ...
  )
}

skip_if_no_jags <- function() {
  if (!requireNamespace("rjags", quietly = TRUE) ||
      !requireNamespace("jagsUI", quietly = TRUE)) {
    testthat::skip("JAGS is not available")
  }
}
