result_fixture <- function() {
  y <- matrix(c(2, 4, 6, NA), nrow = 2)
  theta <- matrix(c(1, 3, 5, 7), nrow = 2)
  loglik <- array(rep(c(-1000, -1001, -999), 4), c(3, 2, 2))
  summary <- matrix(c(0.4, 0.1, 0.2, 0.6), nrow = 1,
                    dimnames = list("alpha[1]", c("mean", "sd", "2.5%", "97.5%")))
  jags_data <- list(y = y, Nparticipants = 2, Ntrials = 2)
  attr(jags_data, "participants") <- c("a", "b")
  attr(jags_data, "row_index") <- matrix(c(1, 2, 3, NA), nrow = 2)
  structure(list(jags_data = jags_data, model_type = "1cs", rating_scale = c(0, 100),
                 convergence = list(converged = FALSE),
                 jags_fit = list(mean = list(theta = theta, loglik = apply(loglik, c(2, 3), mean)),
                                 summary = summary, DIC = 12,
                                 sims.list = list(loglik = loglik,
                                                  gp = matrix(c(1, 1, 3, 3, 3, 3), ncol = 2),
                                                  sigma = matrix(rep(c(1, 2), each = 3), ncol = 2)))),
            class = "mfeargen")
}

test_that("fitted values are model expectations with scale and padding preserved", {
  fit <- result_fixture()
  expected <- matrix(c(10, 30, 50, NA), nrow = 2, dimnames = list(c("a", "b"), NULL))
  expect_equal(fitted(fit), expected)
  expect_equal(unname(residuals(fit)), matrix(c(10, 10, 10, NA), nrow = 2))
  expect_equal(unname(residuals(fit, type = "standardized")), matrix(c(0.6, 1, 0.6, NA), nrow = 2))
  expect_error(predict(fit, newdata = data.frame(x = 1)), "newdata")
  fit$jags_fit$mean$theta <- NULL
  expect_error(fitted(fit), "not monitored")
})

test_that("WAIC is stable and excludes missing observations", {
  fit <- result_fixture()
  result <- calculate_waic(fit)
  ll <- c(-1000, -1001, -999)
  expected_lppd <- -999 + log(mean(exp(ll + 999)))
  expect_equal(result$lppd, 3 * expected_lppd)
  expect_equal(result$p_waic, 3)
  expect_true(is.finite(result$waic))
  expect_length(result$pointwise, 3)
  expect_equal(attr(logLik(fit), "nobs"), 3L)
})

test_that("comparison requires identical observations and prints without recursion", {
  fit <- result_fixture()
  other <- fit
  other$jags_fit$DIC <- 14
  comparison <- compare_models(first = fit, second = other)
  expect_equal(comparison$Delta_DIC, c(0, 2))
  expect_output(print(comparison), "Model comparison")
  expect_equal(attr(compare_models(fit, other, criterion = "WAIC"), "criterion"), "WAIC")
  other$jags_data$y[1] <- 8
  expect_error(compare_models(fit, other), "same observations")
  expect_error(compare_models(fit, fit, criterion = "LOO"), "arg")
})

test_that("recovery accuracy uses matching labels when groups are absent", {
  truth <- data.frame(participant = c("a", "b"), group = c(1L, 3L),
                      alpha_true = c(0, 0.5), lambda_true = c(0, 0.2),
                      w0_true = c(0, 0), w1_true = c(0, 2))
  recovered <- data.frame(participant = c("b", "a"), group_recovered = c(4L, 2L),
                          alpha_recovered = c(0.5, 0), lambda_recovered = c(0.2, 0),
                          w0_recovered = c(0, 0), w1_recovered = c(2, 0))
  result <- mfeargen:::calculate_recovery_metrics(truth, recovered)
  expect_equal(result$group_accuracy, 0)
  expect_equal(dim(result$confusion_matrix), c(4L, 4L))
  expect_equal(result$alpha_rmse, 0)
  expect_true(is.na(result$alpha_correlation))
})

test_that("paper helpers use actual data and fail before unsupported replication", {
  expect_s3_class(load_paper_data("exp1"), "data.frame")
  expect_s3_class(get_paper_priors(), "mfeargen_priors")
  for (experiment in c("exp1", "exp2")) {
    data <- mfeargen:::prepare_paper_trials(experiment)
    prepared <- prepare_mfeargen_data(data, if (experiment == "exp1") "1cs" else "2cs")
    expect_equal(prepared$Nparticipants, 40)
    expect_equal(attr(prepared, "rating_scale"), c(1, 10))
  }
  expect_error(replicate_paper_analysis(model_type_paper = TRUE), "Exact paper replication")
  expect_error(replicate_paper_analysis("exp2", cs_reference = "csp"), "jointly")
})

test_that("held-out ratings receive finite predictions without training leakage", {
  skip_if_not_installed("rjags")
  data <- simulate_mfeargen(n_participants = 2, n_conditioning_trials = 4,
                           n_generalization_trials = 4, seed = 19)
  data$participant <- factor(data$participant, levels = c(1, 2, 3))
  set.seed(15)
  state <- .Random.seed
  cv <- cross_validate_mfeargen(data, k_folds = 2, seed = 22,
                               n.iter = 70, n.burnin = 30, n.thin = 1,
                               n.chains = 1, n.adapt = 100,
                               parallel = FALSE, verbose = FALSE)
  expect_identical(.Random.seed, state)
  expect_true(all(is.finite(cv$fold_metrics$loglik_test)))
  expect_true(all(is.finite(cv$fold_metrics$mse_test)))
  expect_equal(cv$summary$successful_folds, 2)
  for (fold in cv$fold_results) {
    expect_true(all(is.na(fold$fit$data$US_expect[fold$test_rows])))
    expect_true(all(is.finite(fold$predictions)))
  }
})


test_that("scalar participant coefficients are retained", {
  fit <- result_fixture()
  rownames(fit$jags_fit$summary) <- "alpha"
  expect_equal(coef(fit)$parameter, "alpha")
  expect_equal(coef(fit)$estimate, 0.4)
})


test_that("cross-validation rejects prepared data that would leak held-out ratings", {
  data <- simulate_mfeargen(n_participants = 1, n_conditioning_trials = 4,
                           n_generalization_trials = 1, seed = 12)
  expect_error(cross_validate_mfeargen(data, k_folds = "two"), "integer")
  expect_error(cross_validate_mfeargen(data, k_folds = 2,
                                      prepared_data = prepare_mfeargen_data(data)), "reused across folds")
})

test_that("group plots display posterior probabilities for each participant", {
  fit <- result_fixture()
  p <- plot_group_probabilities(fit)
  expect_equal(nrow(p$data), 8L)
  expect_equal(as.vector(tapply(p$data$probability, p$data$participant, sum)), c(1, 1))
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})
