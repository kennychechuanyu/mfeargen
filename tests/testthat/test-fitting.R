test_that("group assignments use posterior probabilities and original IDs", {
  samples <- cbind(c(rep(1, 6), rep(4, 4)), rep(3, 10))
  membership <- extract_group_membership(list(sims.list = list(gp = samples)), c("A", "B"))
  expect_equal(membership$participant, c("A", "B"))
  expect_equal(membership$group, c(1L, 3L))
  expect_equal(membership$max_probability, c(0.6, 1))
  expect_equal(rowSums(membership[paste0("prob_group_", 1:4)]), c(1, 1))
  single <- extract_group_membership(list(sims.list = list(gp = samples[, 1])))
  expect_equal(single$group, 1L)
})

test_that("invalid sampling settings fail before starting JAGS", {
  data <- create_test_data(n_participants = 2, n_trials = 10)
  expect_error(fit_mfeargen(data, n.chains = 0), "n.chains")
  expect_error(fit_mfeargen(data, n.thin = NA_real_), "n.thin")
  expect_error(fit_mfeargen(data, n.iter = 20, n.burnin = 20), "greater")
  expect_error(fit_mfeargen(data, n.iter = 21, n.burnin = 20), "at least two")
  expect_error(fit_mfeargen(data, seed = -1), "seed")
  expect_error(fit_mfeargen(data, prepared_data = list()), "missing")
})

test_that("both designs fit with custom priors and reproducible posterior draws", {
  skip_if_no_jags()
  for (model_type in c("1cs", "2cs")) {
    data <- simulate_mfeargen(n_participants = 2, n_conditioning_trials = 8,
                             n_generalization_trials = 1, model_type = model_type, seed = 7)
    data$participant <- paste0("subject-", data$participant)
    priors <- set_priors_flexible(
      lambda_mu = prior_normal(0.1, 1, truncate = c(0, Inf)),
      w0_sigma = prior_half_cauchy(2), w1_a = prior_half_cauchy(2),
      w1_b = prior_half_cauchy(2), alpha_kappa = prior_uniform(20, 30)
    )
    fit_once <- function() suppressWarnings(fit_mfeargen(
      data, model_type = model_type, priors = priors,
      n.iter = 70, n.burnin = 20, n.chains = 2, n.thin = 1,
      seed = 42, n.adapt = 100, parallel = FALSE, verbose = FALSE
    ))
    set.seed(19)
    state <- .Random.seed
    fit <- fit_once()
    expect_identical(.Random.seed, state)
    repeat_fit <- fit_once()
    expect_identical(fit$jags_fit$sims.list, repeat_fit$jags_fit$sims.list)
    expect_s3_class(fit, "mfeargen")
    expect_equal(fit$group_membership$participant, unique(data$participant))
    expect_true(all(c("theta", "loglik", "y_pre") %in% names(fit$jags_fit$sims.list)))
    expect_true(all(is.finite(fit$jags_fit$sims.list$loglik)))
    expect_true(all(fit$jags_fit$sims.list$alpha_kappa >= 20 &
                    fit$jags_fit$sims.list$alpha_kappa <= 30))
    expect_equal(dim(fit$jags_fit$sims.list$theta)[1], 100)
    expect_true(all(rowSums(fit$group_membership[paste0("prob_group_", 1:4)]) == 1))
  }
})


test_that("manually prepared data cannot bypass response scale normalization", {
  data <- create_test_data()
  prepared <- prepare_mfeargen_data(data)
  prepared$rating_scale <- c(0, 100)
  expect_error(fit_mfeargen(data, prepared_data = prepared), "normalized 0-10")
})

test_that("the full study trial graph compiles and uses group-specific observation noise", {
  skip_if_no_jags()
  data <- prepare_paper_trials("exp1")
  prepared <- prepare_mfeargen_data(data, "1cs", rating_scale = c(1, 10))
  model <- rjags::jags.model(
    system.file("jags", "model_1cs.txt", package = "mfeargen"),
    data = prepared, inits = initialize_mcmc_parameters(prepared, 1),
    n.chains = 1, n.adapt = 0, quiet = TRUE
  )
  samples <- as.matrix(suppressWarnings(rjags::coda.samples(
    model, variable.names = c("gp", "sigma", "tau_obs"),
    n.iter = 2, progress.bar = "none"
  ))[[1]])
  expect_equal(prepared$Nparticipants, 40)
  expect_gt(prepared$Ntrials, 150)
  groups <- samples[, paste0("gp[", seq_len(prepared$Nparticipants), "]"), drop = FALSE]
  precision <- samples[, paste0("tau_obs[", seq_len(prepared$Nparticipants), "]"), drop = FALSE]
  expected <- 1 / ifelse(groups == 1, samples[, "sigma[2]"], samples[, "sigma[1]"])^2
  expect_equal(unname(precision), unname(expected))
})
