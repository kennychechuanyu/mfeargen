repeated_test_data <- function() {
  first <- simulate_mfeargen(n_participants = 2, n_conditioning_trials = 6,
                             n_generalization_trials = 1, seed = 1)
  second <- simulate_mfeargen(n_participants = 2, n_conditioning_trials = 6,
                              n_generalization_trials = 1, seed = 2)
  first$session <- "baseline"
  second$session <- "follow-up"
  data <- rbind(first, second)
  data$participant <- c("subject-b", "subject-a")[data$participant]
  data[!(data$participant == "subject-a" & data$session == "follow-up"), ]
}

test_that("hierarchical data retain sessions, learning masks, and rating units", {
  data <- repeated_test_data()
  prepared <- prepare_hierarchical_data(data, rating_scale = c(1, 10))
  expect_identical(attr(prepared, "participants"), c("subject-b", "subject-a"))
  expect_equal(prepared$Nsessions, c(2L, 1L))
  expect_equal(dim(prepared$y), c(2L, 2L, 16L))
  expect_equal(prepared$Ntrials, matrix(c(16L, 16L, 16L, 0L), nrow = 2))
  expect_true(all(is.na(prepared$y[2, 2, ])))
  expect_true(all(prepared$k[1, 1, 7:16] == 0))
  expect_equal(prepared$y[1, 1, ], (data$US_expect[1:16] - 1) / 9 * 10)
  expect_equal(attr(prepared, "rating_scale"), c(1, 10))
})

test_that("unsupported hierarchical requests fail before sampling", {
  data <- repeated_test_data()
  expect_error(fit_mfeargen_hierarchical(data, model_type = "2cs"), "only.*1cs")
  expect_error(fit_mfeargen_hierarchical(data, priors = list()), "Custom priors")
  expect_error(fit_mfeargen_hierarchical(data, n.iter = 10, n.burnin = 10), "retain")
  expect_error(prepare_hierarchical_data(subset(data, session == "baseline")), "repeated sessions")
})

test_that("variance summaries use draw-wise ratios", {
  fit <- structure(list(jags_fit = list(sims.list = list(
    tau_alpha_sq = c(1, 9), sigma_alpha_sq = c(1, 1)
  ))), class = "mfeargen_hierarchical")
  summary <- calculate_stability(fit, "alpha", credible_interval = 0.8)
  expect_equal(summary$samples, c(0.5, 0.9))
  expect_equal(summary$mean_stability, 0.7)
  expect_equal(summary$credible_interval, 0.8)
  decomposition <- calculate_variance_decomposition(fit, parameters = "alpha")
  expect_equal(decomposition$stability$alpha, 0.7)
  expect_equal(decomposition$between_variance$alpha, 5)
  expect_error(calculate_stability(fit, "alpha", credible_interval = 1), "credible_interval")
  expect_error(calculate_stability(fit, "unknown"), "arg")
})

test_that("session deviations subtract matching stable posterior draws", {
  fit <- structure(list(
    data_info = list(n_participants = 2L, participants = c("b", "a"),
                     sessions = list(c("x", "y"), "x")),
    jags_data = list(Nsessions = c(2L, 1L)),
    jags_fit = list(sims.list = list(
      alpha_stable = matrix(c(0.2, 0.4, 0.6, 0.8), nrow = 2),
      alpha = array(c(0.3, 0.5, 0.9, 1, 0.1, 0.2, NA, NA), c(2, 2, 2))
    ))
  ), class = "mfeargen_hierarchical")
  stable <- extract_stable_components(fit)
  expect_identical(rownames(stable), c("b", "a"))
  expect_equal(stable[1, ], c(0.2, 0.4))
  deviations <- extract_session_deviations(fit)
  expect_equal(as.numeric(deviations[1, 1, ]), c(0.1, 0.1))
  expect_equal(as.numeric(deviations[1, 2, ]), c(-0.1, -0.2))
  expect_true(all(is.na(deviations[2, 2, ])))
})

test_that("the repeated-session JAGS model runs and plots actual sessions", {
  skip_if_no_jags()
  data <- repeated_test_data()
  set.seed(27)
  rng_before <- .Random.seed
  fit <- fit_mfeargen_hierarchical(data, rating_scale = c(1, 10),
    n.iter = 80, n.burnin = 40, n.thin = 1, n.chains = 2,
    n.adapt = 100, parallel = FALSE, verbose = FALSE, seed = 20)
  expect_identical(.Random.seed, rng_before)
  expect_s3_class(fit, "mfeargen_hierarchical")
  expect_equal(fit$data_info$total_sessions, 3)
  expect_identical(fit$group_membership$participant, c("subject-b", "subject-a"))
  expect_equal(dim(extract_stable_components(fit)), c(2L, 80L))
  expect_equal(dim(extract_session_deviations(fit)), c(2L, 2L, 80L))
  alpha <- fit$jags_fit$sims.list$alpha
  expect_true(all(alpha[is.finite(alpha)] >= 0 & alpha[is.finite(alpha)] <= 1))
  expect_false(fit$convergence$converged)
  set.seed(19)
  old_seed <- .Random.seed
  plot <- plot_hierarchical_structure(fit)
  expect_identical(.Random.seed, old_seed)
  expect_equal(nrow(plot$data), 3)
  expect_setequal(plot$data$session, c("baseline", "follow-up"))
  expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
})
