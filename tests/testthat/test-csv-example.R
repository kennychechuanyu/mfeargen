test_that("the installed CSV workflow preserves IDs and writes reviewable results", {
  skip_if_no_jags()
  script <- new.env(parent = globalenv())
  source(system.file("examples", "analyze_trials.R", package = "mfeargen"), local = script)
  data <- simulate_mfeargen(n_participants = 2, n_conditioning_trials = 4,
                            n_generalization_trials = 1, rating_scale = c(0, 100), seed = 18)
  data$participant <- c("001", "002")[data$participant]
  data$US_expect[5] <- NA_real_
  data <- data[rev(seq_len(nrow(data))), ]
  input <- tempfile(fileext = ".csv")
  output <- tempfile("mfeargen-csv-")
  on.exit(unlink(c(input, output), recursive = TRUE), add = TRUE)
  utils::write.csv(data, input, row.names = FALSE)
  fit <- suppressMessages(script$analyze_trials(
    input, rating_scale = c(0, 100), output_dir = output,
    n.iter = 60, n.burnin = 30, n.chains = 2, n.thin = 1, seed = 22
  ))
  expected <- c("parameters.csv", "groups.csv", "expected-ratings.csv",
    "predictive-intervals.csv", "diagnostics.csv", "fit.rds",
    "group-probabilities.png", "posterior-predictive.png")
  expect_true(all(file.exists(file.path(output, expected))))
  ratings <- utils::read.csv(file.path(output, "expected-ratings.csv"),
                             colClasses = c(participant = "character"))
  diagnostics <- utils::read.csv(file.path(output, "diagnostics.csv"))
  expect_identical(ratings$participant, data$participant)
  expect_identical(ratings$trials, data$trials)
  expect_equal(ratings$observed, data$US_expect)
  expect_true(all(is.finite(ratings$estimate)))
  expect_true(all(ratings$lower <= ratings$upper))
  expect_true(all(c("rhat", "ess_bulk", "ess_tail", "status") %in% names(diagnostics)))
  expect_gt(nrow(diagnostics), 0)
  expect_identical(readRDS(file.path(output, "fit.rds"))$observation_model, "censored")
  expect_identical(fit$group_membership$participant, unique(data$participant))
})
