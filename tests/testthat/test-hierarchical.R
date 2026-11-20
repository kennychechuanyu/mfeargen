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
  expect_error(fit_mfeargen_hierarchical(data, save_predictive = NA), "save_predictive")
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
  data$US_expect[c(1L, nrow(data))] <- c(1, 10)
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
  expect_equal(fit$observation_model, "censored")
  expect_equal(dim(fitted(fit)), c(2L, 2L, 16L))
  expect_true(all(is.na(fitted(fit)[2, 2, ])))
  expect_false("y_pre" %in% names(fit$jags_fit$sims.list))
  set.seed(15)
  state <- .Random.seed
  predictions <- predict(fit, type = "response", interval = "prediction", format = "data.frame", seed = 19)
  expect_identical(.Random.seed, state)
  expect_identical(predictions, predict(fit, type = "response", interval = "prediction", format = "data.frame", seed = 19))
  expect_equal(predictions$participant, data$participant)
  expect_equal(predictions$session, data$session)
  expect_equal(predictions$trials, data$trials)
  expect_true(all(predictions$lower >= 1 & predictions$upper <= 10))
  expect_true(all(is.finite(residuals(fit)[!is.na(fit$jags_data$y)])))
  expect_true(all(c("estimate", "posterior_sd", "participant", "session") %in% names(coef(fit))))
  raw_coefficients <- coef(fit)
  normalized_coefficients <- coef(fit, scale = "normalized")
  noise <- raw_coefficients$component == "observation_noise"
  expect_equal(raw_coefficients$estimate[noise], normalized_coefficients$estimate[noise] * 0.9)
  boundary_cells <- which(fit$jags_data$y %in% c(0, 10), arr.ind = FALSE)
  response <- fit$jags_data$y[boundary_cells]
  theta <- matrix(fit$jags_fit$sims.list$theta, nrow = 80)[1, boundary_cells]
  sd <- hierarchical_sigma_draws(fit, 80)[1, ((boundary_cells - 1L) %% 2L) + 1L]
  expected_loglik <- ifelse(response == 0, stats::pnorm(0, theta, sd, log.p = TRUE),
                            stats::pnorm(10, theta, sd, lower.tail = FALSE, log.p = TRUE))
  actual_loglik <- matrix(fit$jags_fit$sims.list$loglik, nrow = 80)[1, boundary_cells]
  expect_equal(actual_loglik, expected_loglik, tolerance = 1e-8)

  set.seed(19)
  old_seed <- .Random.seed
  plot <- plot_hierarchical_structure(fit)
  expect_identical(.Random.seed, old_seed)
  expect_equal(nrow(plot$data), 3)
  expect_setequal(plot$data$session, c("baseline", "follow-up"))
  expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
})

test_that("hierarchical predictions keep original rows, units, and uncertainty", {
  y <- array(c(3, 2, NA, 8), c(1, 2, 2))
  jags_data <- list(y = y, Nparticipants = 1, Nsessions = 2L,
                    Ntrials = matrix(2L, 1, 2))
  attr(jags_data, "row_index") <- array(c(3L, 1L, 4L, 2L), c(1, 2, 2))
  theta <- array(rep(1:4, 4), c(4, 1, 2, 2))
  future <- array(rep(c(-2, 0, 10, 12), 4), c(4, 1, 2, 2))
  fit <- structure(list(
    jags_data = jags_data, rating_scale = c(10, 20),
    data_info = list(n_participants = 1L, participants = "b",
                     sessions = list(c("baseline", "follow-up"))),
    data = data.frame(participant = "b", session = c("follow-up", "follow-up", "baseline", "baseline"),
                       trials = c(1L, 2L, 1L, 2L), US_expect = c(12, 18, 13, NA)),
    jags_fit = list(mean = list(theta = apply(theta, c(2, 3, 4), mean)),
                    sims.list = list(theta = theta, y_pre = future,
                      sigma_obs = matrix(rep(c(1, 2), each = 4), ncol = 2),
                      gp = c(1, 1, 3, 3)))
  ), class = "mfeargen_hierarchical")
  expect_equal(as.numeric(fitted(fit)), rep(12.5, 4))
  predictions <- predict(fit, type = "response", interval = "credible", format = "data.frame")
  expect_equal(predictions[c("participant", "session", "trials")], fit$data[c("participant", "session", "trials")])
  expect_equal(predictions$observed, fit$data$US_expect)
  expect_equal(predictions$estimate, rep(12.5, 4))
  expect_equal(predictions$lower, rep(11.075, 4))
  expect_equal(predictions$upper, rep(13.925, 4))
  future <- predict(fit, type = "response", interval = "prediction", format = "data.frame")
  expect_equal(future$lower, rep(8.15, 4))
  expect_equal(future$upper, rep(21.85, 4))
  expect_equal(as.numeric(residuals(fit, type = "standardized")), c(0.5, -0.5, NA, 5.5) / 1.5)
  expect_equal(dim(predict(fit, type = "response", interval = "credible")), c(1L, 2L, 2L, 3L))
  expect_error(predict(fit, newdata = data.frame(x = 1)), "newdata")
})

test_that("hierarchical preparation validates the original identifiers", {
  data <- repeated_test_data()
  data$session[1] <- ""
  expect_error(prepare_hierarchical_data(data), "session.*identifiers")
  data <- repeated_test_data()
  data$participant <- as.numeric(factor(data$participant))
  data$participant[1] <- Inf
  expect_error(prepare_hierarchical_data(data), "participant.*identifiers")
})

test_that("single-participant repeated sessions support normal response exports", {
  skip_if_no_jags()
  data <- subset(repeated_test_data(), participant == "subject-b")
  data <- data[rev(seq_len(nrow(data))), ]
  fit <- fit_mfeargen_hierarchical(data, observation_model = "normal", rating_scale = c(1, 10),
    n.iter = 40, n.burnin = 20, n.thin = 1, n.chains = 1,
    n.adapt = 100, parallel = FALSE, verbose = FALSE, seed = 24, save_predictive = TRUE)
  expect_true("y_pre" %in% names(fit$jags_fit$sims.list))
  expect_identical(fit$convergence$status, "unassessed")
  expect_true(is.na(fit$convergence$converged))
  expect_equal(fit$convergence$n_converged, 0)
  result <- predict(fit, type = "response", interval = "credible", format = "data.frame")
  expect_identical(result$participant, data$participant)
  expect_identical(result$session, data$session)
  expect_identical(result$trials, data$trials)
  expect_equal(result$observed, data$US_expect)
  expect_true(all(is.finite(result$estimate)))
  coefficients <- coef(fit)
  stable <- coefficients$parameter == "alpha_stable"
  expect_identical(coefficients$participant[stable], "subject-b")
  expect_equal(nrow(predict(fit, type = "groups")), 1L)
})
