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
      seed = 42, n.adapt = 100, parallel = FALSE, verbose = FALSE, save_predictive = TRUE
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
    expect_equal(rowSums(fit$group_membership[paste0("prob_group_", 1:4)]), rep(1, nrow(fit$group_membership)))
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


test_that("prepared histories cannot be reused after the input data change", {
  data <- create_test_data(n_participants = 2)
  prepared <- prepare_mfeargen_data(data)
  changed <- data
  changed$US_expect[1] <- changed$US_expect[1] + 0.1
  expect_error(fit_mfeargen(changed, prepared_data = prepared), "does not match data")
  changed <- data[nrow(data):1, ]
  expect_error(fit_mfeargen(changed, prepared_data = prepared), "does not match data")
  invalid <- prepared
  invalid$k[1, 1] <- 2
  expect_error(fit_mfeargen(data, prepared_data = invalid), "coded 0 or 1")
  invalid <- prepared
  invalid$d_phy[1, 1] <- -1
  expect_error(fit_mfeargen(data, prepared_data = invalid), "non-negative")
  expect_error(fit_mfeargen(data, n_iter = 200), "Unknown JAGS arguments")
  expect_error(fit_mfeargen(data, DIC = NA), "DIC must")
})

test_that("censored likelihoods include boundary mass and correct observed deviance", {
  skip_if_no_jags()
  model_text <- c(
    "model {", "  for (i in 1:Nparticipants) {",
    "    tau_obs[i] <- pow(noise[i], -2)",
    "    for (j in 1:Ntrials) {",
    "      theta[i,j] <- location[i,j]",
    "      y[i,j] ~ dnorm(theta[i,j], tau_obs[i])",
    "      loglik[i,j] <- logdensity.norm(y[i,j], theta[i,j], tau_obs[i])",
    "      y_pre[i,j] ~ dnorm(theta[i,j], tau_obs[i])",
    "    }", "  }", "}"
  )
  observed <- c(0, 10, 5, NA, 0, 0, 10, 0, 0)
  location <- c(1, 9, 5, 5, 10, 10, 0, 1, 1)
  noise <- c(1, 1, 2, 1, 0.1, 1e-5, 1e-5, 0.1, 1 / 10.01)
  data <- list(y = matrix(observed, ncol = 1), location = matrix(location, ncol = 1),
               noise = noise, Nparticipants = 9L, Ntrials = 1L, rating_scale = c(0, 10))
  compiled <- compile_observation_model(model_text, data, "censored")
  rjags::load.module("dic", quiet = TRUE)
  compiled$data$rating_scale <- NULL
  model <- rjags::jags.model(textConnection(compiled$model_code), data = compiled$data,
                             n.chains = 1, n.adapt = 0, quiet = TRUE)
  samples <- as.matrix(rjags::coda.samples(model, c("loglik", "deviance", "y_pre"),
                                           n.iter = 2, progress.bar = "none")[[1]])
  expected <- c(stats::pnorm(-1, log.p = TRUE), stats::pnorm(-1, log.p = TRUE),
                stats::dnorm(5, 5, 2, log = TRUE), 0,
                stats::pnorm(-100, log.p = TRUE),
                stats::pnorm(-1e6, log.p = TRUE), stats::pnorm(-1e6, log.p = TRUE),
                stats::pnorm(-10, log.p = TRUE), stats::pnorm(-10.01, log.p = TRUE))
  likelihood <- as.vector(samples[1, paste0("loglik[", 1:9, ",1]")])
  expect_lt(max(abs(likelihood - expected) / pmax(1, abs(expected))), 1e-10)
  expect_equal(samples[, "deviance"], rep(-2 * sum(expected), 2))
  predictions <- samples[, paste0("y_pre[", 1:9, ",1]"), drop = FALSE]
  expect_true(all(predictions >= 0 & predictions <= 10))
  expect_true(all(is.na(compiled$data$y_cont[c(1, 2, 4, 5, 6, 7)])))
  expect_identical(data$y, matrix(observed, ncol = 1))
})

test_that("both designs honor censored simulation metadata", {
  skip_if_no_jags()
  for (model_type in c("1cs", "2cs")) {
    data <- simulate_mfeargen(n_participants = 2, n_conditioning_trials = 8,
                              n_generalization_trials = 1, model_type = model_type, seed = 7)
    data$US_expect[1:3] <- c(1, 10, NA)
    fit <- fit_mfeargen(data, model_type = model_type,
                        n.iter = 60, n.burnin = 20, n.chains = 1, n.thin = 1,
                        n.adapt = 100, parallel = FALSE, verbose = FALSE, seed = 27, save_predictive = TRUE)
    expect_identical(fit$observation_model, "censored")
    expect_true(is.finite(fit$jags_fit$DIC))
    expect_equal(fit$jags_fit$mean$loglik[1, 3], 0)
    expect_true(is.na(fit$jags_data$y[1, 3]))
    expect_true(all(fit$jags_fit$sims.list$y_pre >= 0 & fit$jags_fit$sims.list$y_pre <= 10))
    expect_true(is.na(fit$convergence$converged))
  }
})

test_that("convergence diagnostics require independent chains and effective samples", {
  set.seed(29)
  independent <- array(stats::rnorm(12000), c(1000, 4, 3),
                        dimnames = list(NULL, NULL, c("alpha[1]", "lambda[1]", "y_pre[1,1]")))
  make_fit <- function(draws) list(samples = lapply(seq_len(dim(draws)[2]), function(chain) {
    value <- draws[, chain, , drop = FALSE]
    matrix(value, nrow = dim(draws)[1], dimnames = list(NULL, dimnames(draws)[[3]]))
  }))
  fit <- make_fit(independent)
  diagnostics <- check_convergence(fit, verbose = FALSE)
  expect_true(diagnostics$converged)
  expect_equal(diagnostics$n_parameters, 2)
  expect_false(any(grepl("y_pre", diagnostics$diagnostics$parameter)))
  fit$samples <- fit$samples[1]
  expect_true(is.na(check_convergence(fit, verbose = FALSE)$converged))
  stuck <- independent
  for (chain in 1:4) stuck[, chain, 1] <- chain
  expect_false(check_convergence(make_fit(stuck), verbose = FALSE)$converged)
  expect_false(check_convergence(c(parameter = Inf), verbose = FALSE)$converged)
  expect_true(is.na(check_convergence(c(parameter = 1), verbose = FALSE)$converged))
  few <- independent[1:20, , , drop = FALSE]
  expect_false(check_convergence(make_fit(few), verbose = FALSE)$converged)
  constant <- independent
  constant[, , 1] <- 0
  answer <- check_convergence(make_fit(constant), verbose = FALSE)
  expect_equal(answer$structural_constants, "alpha[1]")
  expect_true(answer$converged)
})


test_that("singleton panels retain posterior dimensions without caching replicates", {
  skip_if_no_jags()
  data <- data.frame(participant = "one", trials = 1L, US_expect = 0,
                      US = 0, k = 0, d_phy = 1, d_per_old = 1)
  for (observation_model in c("normal", "censored")) {
    fit <- fit_mfeargen(data, rating_scale = c(0, 10), observation_model = observation_model,
                        n.iter = 50, n.burnin = 20, n.chains = 1, n.thin = 1,
                        n.adapt = 100, parallel = FALSE, verbose = FALSE, seed = 3)
    expect_equal(dim(fit$jags_fit$sims.list$theta), c(30L, 1L, 1L))
    expect_equal(dim(fit$jags_fit$sims.list$loglik), c(30L, 1L, 1L))
    expect_equal(dim(fit$jags_fit$mean$theta), c(1L, 1L))
    expect_false("y_pre" %in% names(fit$jags_fit$sims.list))
    expect_false(fit$save_predictive)
  }
})


test_that("parameter extraction handles summary-only outputs and named singleton estimates", {
  summary <- matrix(c(0.4, 0.1, 0.2, 0.4, 0.6, 0.3, 0.1, 0.1, 0.3, 0.5),
                     nrow = 2, byrow = TRUE,
                     dimnames = list(c("alpha", "lambda"), c("mean", "sd", "2.5%", "50%", "97.5%")))
  estimates <- extract_parameters(list(summary = summary), "alpha")
  expect_identical(rownames(estimates), "alpha")
  expect_identical(estimates$parameter, "alpha")
  expect_equal(estimates$mean, 0.4)
  expect_equal(estimates$q2.5, 0.2)
  expect_equal(estimates$q97.5, 0.6)
  expect_error(extract_parameters(list(summary = summary), summary_only = NA), "summary_only")
  expect_error(extract_parameters(list(summary = summary), "alpha", summary_only = FALSE), "Retained posterior draws")
  skip_if_no_jags()
  data <- create_test_data(n_participants = 1)
  fit <- fit_mfeargen(data, n.iter = 50, n.burnin = 20, n.chains = 1, n.thin = 1,
                      n.adapt = 100, parallel = FALSE, verbose = FALSE, seed = 3)
  estimates <- extract_parameters(fit, "alpha")
  expect_identical(estimates$parameter, "alpha")
  expect_equal(estimates$mean, as.vector(fit$jags_fit$mean$alpha))
  expect_equal(as.vector(extract_parameters(fit, "alpha", summary_only = FALSE)),
               as.vector(fit$jags_fit$sims.list$alpha))
})
