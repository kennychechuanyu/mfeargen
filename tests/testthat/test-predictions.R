prediction_fixture <- function(observation_model = "normal") {
  n <- 200L
  rows <- matrix(c(4L, 3L, 1L, 2L, NA, 5L), nrow = 2)
  data <- data.frame(participant = c("B", "A", "A", "B", "A"),
                     trials = c(9L, 40L, 20L, 3L, 60L),
                     US_expect = c(40, 50, NA, 10, 70))
  y <- matrix(NA_real_, 2, 3)
  y[!is.na(rows)] <- data$US_expect[rows[!is.na(rows)]] / 10
  theta <- array(rep(c(3, 4, 5, 6, 7, 8), each = n), c(n, 2, 3))
  theta <- theta + seq(-0.2, 0.2, length.out = n)
  predictions <- theta + stats::qnorm((seq_len(n) - 0.5) / n)
  prepared <- list(y = y, Nparticipants = 2L, Ntrials = 3L)
  attr(prepared, "participants") <- c("B", "A")
  attr(prepared, "row_index") <- rows
  attr(prepared, "rating_scale") <- c(0, 100)
  structure(list(data = data, jags_data = prepared, rating_scale = c(0, 100),
                 observation_model = observation_model,
                 jags_fit = list(mean = list(theta = apply(theta, c(2, 3), mean)),
                                 sims.list = list(theta = theta, y_pre = predictions,
                                                  gp = matrix(3L, n, 2),
                                                  sigma = matrix(rep(c(1, 2), each = n), ncol = 2)))),
            class = "mfeargen")
}

test_that("prediction exports retain original row order, trial labels and missing observations", {
  fit <- prediction_fixture()
  result <- predict(fit, type = "response", interval = "credible", format = "data.frame")
  expect_equal(result$participant, fit$data$participant)
  expect_equal(result$trials, fit$data$trials)
  expect_equal(result$observed, fit$data$US_expect)
  expect_equal(result$estimate, c(50, 60, 40, 30, 80))
  expect_true(all(result$lower <= result$estimate & result$estimate <= result$upper))
  expect_equal(nrow(result), 5L)
  normalized <- predict(fit, type = "response", format = "data.frame", scale = "normalized")
  expect_equal(normalized$observed, fit$data$US_expect / 10)
  expect_equal(normalized$estimate, result$estimate / 10)
  intervals <- predict(fit, type = "response", interval = "credible")
  expect_equal(dim(intervals), c(2L, 3L, 3L))
  expect_true(all(is.na(intervals[1, 3, ])))
  expect_error(predict(fit, type = "response", level = 1), "level")
  expect_error(predict(fit, interval = "prediction"), "type")
})

test_that("predictive intervals include observation noise and preserve RNG state", {
  fit <- prediction_fixture()
  credible <- predict(fit, type = "response", interval = "credible", format = "data.frame")
  predictive <- predict(fit, type = "response", interval = "prediction", format = "data.frame")
  expect_true(all(predictive$upper - predictive$lower > credible$upper - credible$lower))
  set.seed(48)
  state <- .Random.seed
  draws <- posterior_predict(fit, n_draws = 15)
  expect_identical(.Random.seed, state)
  expect_equal(dim(draws), c(15L, 2L, 3L))
  expect_equal(dimnames(draws)[[2]], c("B", "A"))
  expect_true(all(is.na(draws[, 1, 3])))
  expect_error(posterior_predict(fit, n_draws = 201), "n_draws")
})

test_that("censored response moments agree with numerical integration", {
  for (location in c(0, 2, 5, 10)) {
    sigma <- 1.5
    mean_integral <- stats::integrate(function(x) x * stats::dnorm(x, location, sigma), 0, 10)$value +
      10 * stats::pnorm(10, location, sigma, lower.tail = FALSE)
    second_integral <- stats::integrate(function(x) x^2 * stats::dnorm(x, location, sigma), 0, 10)$value +
      100 * stats::pnorm(10, location, sigma, lower.tail = FALSE)
    expect_equal(mfeargen:::clipped_normal_mean(location, sigma), mean_integral)
    expect_equal(mfeargen:::clipped_normal_sd(location, sigma), sqrt(second_integral - mean_integral^2))
  }
  expect_equal(mfeargen:::clipped_normal_sd(c(2, 5), c(1e10, 1e15)), c(5, 5), tolerance = 1e-7)
  expect_equal(mfeargen:::clipped_normal_sd(10, 1e-9) / 1e-9,
               sqrt(0.5 - 1 / (2 * pi)), tolerance = 1e-7)
})

test_that("censored expectations integrate observation noise instead of clipping theta", {
  fit <- prediction_fixture("censored")
  fit$jags_fit$sims.list$theta[, 1, 1] <- 0
  expected <- stats::dnorm(0) * 10
  expect_equal(unname(fitted(fit)[1, 1]), expected, tolerance = 1e-6)
  expect_true(all(posterior_predict(fit) >= 0 & posterior_predict(fit) <= 100, na.rm = TRUE))
  expect_true(all(is.finite(residuals(fit, type = "standardized")[!is.na(fit$jags_data$y)])))
})

test_that("predictive checks exclude missing observations and padding", {
  fit <- prediction_fixture("censored")
  checks <- posterior_predictive_check(fit, n_draws = 20)
  expect_equal(attr(checks, "n_observations"), 4L)
  expect_equal(attr(checks, "n_draws"), 20L)
  expect_equal(checks$observed[checks$statistic == "mean"], mean(fit$data$US_expect, na.rm = TRUE))
  expect_true(all(checks$tail_probability >= 0 & checks$tail_probability <= 1))
  expect_s3_class(ggplot2::ggplot_build(plot_posterior_predictive(fit)), "ggplot_built")
})

test_that("endpoint scores are probability masses and interior scores are densities", {
  expect_equal(mfeargen:::response_log_density(0, 1, 2, "censored"),
               stats::pnorm(0, 1, 2, log.p = TRUE))
  expect_equal(mfeargen:::response_log_density(10, 9, 2, "censored"),
               stats::pnorm(10, 9, 2, lower.tail = FALSE, log.p = TRUE))
  expect_equal(mfeargen:::response_log_density(5, 6, 2, "censored"),
               stats::dnorm(5, 6, 2, log = TRUE))
  expect_equal(mfeargen:::response_log_density(0, 1, 2, "normal"),
               stats::dnorm(0, 1, 2, log = TRUE))
})


test_that("on-demand predictive draws are reproducible and use posterior observation noise", {
  fit <- prediction_fixture("censored")
  fit$jags_fit$sims.list$y_pre <- NULL
  set.seed(18)
  state <- .Random.seed
  samples <- posterior_predict(fit, n_draws = 100, seed = 43)
  expect_identical(.Random.seed, state)
  expect_identical(samples, posterior_predict(fit, n_draws = 100, seed = 43))
  expect_true(all(samples >= 0 & samples <= 100, na.rm = TRUE))
  intervals <- predict(fit, type = "response", interval = "prediction", seed = 43)
  expect_identical(intervals, predict(fit, type = "response", interval = "prediction", seed = 43))
  checks <- posterior_predictive_check(fit, seed = 43)
  expect_equal(attr(checks, "n_draws"), 100L)
  expect_identical(checks, posterior_predictive_check(fit, seed = 43))
})


test_that("gradients use physical coordinates and show model uncertainty", {
  fit <- prediction_fixture("censored")
  fit$data$Phy_size <- c(8, 4, 2, 6, 10)
  fit$data$stimulus <- c("GS", "CS+", "GS", "GS", "GS")
  plot <- plot_generalization_gradient(fit, group_separate = FALSE)
  expect_equal(sort(plot$data$stimulus), sort(fit$data$Phy_size))
  expect_true(all(plot$data$lower <= plot$data$predicted & plot$data$predicted <= plot$data$upper))
  expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
  fit$data$Phy_size <- rep(5, nrow(fit$data))
  single <- plot_generalization_gradient(fit, group_separate = FALSE)
  expect_s3_class(ggplot2::ggplot_build(single), "ggplot_built")
})
