test_that("seeded simulations are reproducible without changing the caller RNG", {
  set.seed(113)
  previous <- .Random.seed
  first <- simulate_mfeargen(3, n_conditioning_trials = 4, n_generalization_trials = 1,
                            seed = 19, return_true_params = TRUE)
  expect_identical(.Random.seed, previous)
  second <- simulate_mfeargen(3, n_conditioning_trials = 4, n_generalization_trials = 1,
                             seed = 19, return_true_params = TRUE)
  expect_identical(first, second)
})

test_that("simulation counts reflect conditioning totals and repetitions per stimulus", {
  x <- simulate_mfeargen(2, stimuli = c(1, 3, 5), cs_plus = 3,
                        n_conditioning_trials = 6, n_generalization_trials = 2,
                        reinforcement_rate = 1, seed = 11)
  expect_equal(nrow(x), 24)
  expect_equal(table(x$participant, x$phase)[, "conditioning"], c(`1` = 6L, `2` = 6L))
  expect_true(all(x$US[x$phase == "conditioning"] == 1))
  expect_true(all(x$US[x$phase == "generalization"] == 0))
  expect_true(all(table(x$participant[x$phase == "generalization"],
                        x$stimulus[x$phase == "generalization"]) == 2))
})

test_that("differential conditioning alternates CSs and leaves CS- unreinforced", {
  x <- simulate_mfeargen(2, model_type = "2cs", cs_plus = 7, cs_minus = 3,
                        n_conditioning_trials = 6, n_generalization_trials = 1,
                        reinforcement_rate = 1, seed = 21)
  conditioning <- x[x$phase == "conditioning", ]
  expect_true(all(conditioning$USp[conditioning$stimulus_type == "CS+"] == 1))
  expect_true(all(conditioning$USp[conditioning$stimulus_type == "CS-"] == 0))
  expect_true(all(conditioning$USm[conditioning$stimulus_type == "CS-"] == -1))
  expect_true(all(x$USm[x$phase == "generalization"] == 0))
  prepared <- prepare_mfeargen_data(x, "2cs")
  expect_equal(prepared$lambda_threshold, attr(x, "simulation_info")$lambda_threshold)
  expect_true(all(prepared$k_plus[, 7:ncol(prepared$k_plus)] == 0))
})

test_that("changing rating scale only changes response units", {
  args <- list(n_participants = 2, n_conditioning_trials = 4,
               n_generalization_trials = 1, seed = 33, return_true_params = TRUE)
  a <- do.call(simulate_mfeargen, c(args, list(rating_scale = c(0, 10))))
  b <- do.call(simulate_mfeargen, c(args, list(rating_scale = c(-20, 80))))
  expect_equal(b$US_expect, 10 * a$US_expect - 20)
  expect_equal(attr(b, "true_parameters")$sigma_true, 10 * attr(a, "true_parameters")$sigma_true)
  expect_equal(prepare_mfeargen_data(a)$y, prepare_mfeargen_data(b)$y)
})

test_that("simulation parameters reject invalid and ambiguous designs", {
  expect_error(simulate_mfeargen(n_participants = 1.2), "n_participants")
  expect_error(simulate_mfeargen(n_conditioning_trials = NA_real_), "n_conditioning_trials")
  expect_error(simulate_mfeargen(stimuli = c(5, 5)), "distinct")
  expect_error(simulate_mfeargen(group_proportions = c(0.5, 0.5)), "four")
  expect_error(simulate_mfeargen(group_proportions = c(-1, 1, 0.5, 0.5)), "non-negative")
  expect_error(simulate_mfeargen(rating_scale = c(10, 0)), "increasing")
  expect_error(simulate_mfeargen(reinforcement_rate = 2), "reinforcement_rate")
  expect_error(simulate_mfeargen(perceptual_noise_sd = -1), "perceptual_noise_sd")
  expect_error(simulate_mfeargen(lambda_criterion = 1), "criterion")
  expect_error(simulate_mfeargen(model_type = "2cs", cs_minus = 5), "distinct")
})

test_that("2cs simulation learns an inhibitory gradient around CS-", {
  x <- simulate_participant(
    p = 1, group = 3,
    params = list(alpha = 1, lambda = 2, w0 = 0, w1 = 4, sigma = 0),
    stimuli = 1:3, cs_plus = 3, cs_minus = 1,
    n_conditioning_trials = 4, n_generalization_trials = 1,
    reinforcement_rate = 1, model_type = "2cs", rating_scale = c(0, 10),
    perceptual_noise_sd = 0
  )
  test <- x[x$phase == "generalization", ]
  expect_equal(test$US_expect[test$stimulus == 1], 10 * plogis(4 * (exp(-4) - 1)))
  expect_equal(test$US_expect[test$stimulus == 3], 10 * plogis(4 * (1 - exp(-4))))
})

test_that("observation noise matches the fitted learner and non-learner categories", {
  x <- simulate_mfeargen(16, stimuli = 1:3, cs_plus = 2,
                        n_conditioning_trials = 4, n_generalization_trials = 1,
                        seed = 83, return_true_params = TRUE)
  truth <- attr(x, "true_parameters")
  expect_true(any(truth$group == 1))
  expect_true(any(truth$group != 1))
  expect_length(unique(truth$sigma_true[truth$group == 1]), 1L)
  expect_length(unique(truth$sigma_true[truth$group != 1]), 1L)
  expect_gt(truth$sigma_true[truth$group == 1][1], truth$sigma_true[truth$group != 1][1])
})

test_that("overgeneralizer responses use the same distances as the fitted model", {
  for (design in c("1cs", "2cs")) {
    set.seed(73)
    data <- simulate_participant(
      p = 1, group = 2,
      params = list(alpha = 1, lambda = 0.3, w0 = -0.5, w1 = 2, sigma = 0),
      stimuli = 1:3, cs_plus = 3, cs_minus = 1,
      n_conditioning_trials = 4, n_generalization_trials = 2,
      reinforcement_rate = 1, model_type = design, rating_scale = c(0, 10),
      perceptual_noise_sd = 1
    )
    test <- data[data$phase == "generalization", ]
    strength <- exp(-0.3 * abs(test$Per_size - 3))
    if (design == "2cs") strength <- strength - exp(-0.3 * abs(test$Per_size - 1))
    expect_equal(test$US_expect, 10 * plogis(-0.5 + 2 * strength))
  }
})

test_that("censored and normal simulation differ only at the instrument bounds", {
  args <- list(n_participants = 8, n_conditioning_trials = 8,
               n_generalization_trials = 3, seed = 94,
               group_proportions = c(1, 0, 0, 0))
  normal <- do.call(simulate_mfeargen, c(args, list(observation_model = "normal")))
  censored <- do.call(simulate_mfeargen, c(args, list(observation_model = "censored")))
  expect_true(any(normal$US_expect < 1 | normal$US_expect > 10))
  expect_equal(censored$US_expect, pmin(10, pmax(1, normal$US_expect)))
  expect_equal(normal$Per_size, censored$Per_size)
  expect_identical(attr(normal, "observation_model"), "normal")
  expect_identical(attr(censored, "observation_model"), "censored")
  expect_identical(attr(prepare_mfeargen_data(normal), "observation_model"), "normal")
  expect_equal(as.vector(t(prepare_mfeargen_data(normal)$y)),
               (normal$US_expect - 1) / 9 * 10)
})
