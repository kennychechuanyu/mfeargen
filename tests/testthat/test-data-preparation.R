test_that("trial order and original participant IDs survive matrix preparation", {
  x <- data.frame(
    participant = c("B", "A", "B", "A", "B"),
    trials = c(9L, 4L, 2L, 8L, 5L),
    US_expect = c(80, 30, 20, NA, 50),
    US = c(0, 1, 1, 0, 0),
    stimulus_type = c("GS", "CS+", "CS+", "GS", "CS+"),
    Phy_size = c(8, 5, 5, 7, 5),
    Per_size = c(7.5, 5.1, 4.8, 6.8, 5.3),
    phase = c("generalization", "conditioning", "conditioning", "generalization", "generalization")
  )
  prepared <- prepare_mfeargen_data(x, rating_scale = c(0, 100))
  expect_identical(attr(prepared, "participants"), c("B", "A"))
  expect_identical(attr(prepared, "trial_counts"), c(3L, 2L))
  expect_equal(prepared$y, rbind(c(2, 5, 8), c(3, NA, NA)))
  expect_equal(prepared$k, rbind(c(1, 0, 0), c(1, 0, 0)))
  expect_equal(prepared$r, rbind(c(1, 0, 0), c(1, 0, 0)))
  expect_equal(prepared$d_phy, rbind(c(0, 0, 3), c(0, 2, 0)))
  expect_equal(prepared$d_per_old, rbind(c(0.2, 0.3, 2.5), c(0.1, 1.8, 0)))
  expect_equal(attr(prepared, "row_index"), rbind(c(3L, 5L, 1L), c(2L, 4L, NA)))
  expect_identical(attr(prepared, "observed"), !is.na(prepared$y))
  expect_equal(attr(prepared, "rating_scale"), c(0, 100))
  expect_equal(prepared$rating_scale, c(0, 10))
})

test_that("2cs retains measured perceptual distances and both outcome sequences", {
  x <- data.frame(
    participant = "p7", trials = 1:3, US_expect = c(2, 4, 6),
    USp = c(1, 0, 0), USm = c(0, 1, 0),
    k_plus = c(1, 0, 0), k_minus = c(0, 1, 0),
    d_p_phy = c(0, 4, 2), d_m_phy = c(4, 0, 2),
    d_p_per = c(0.2, 3.1, 1.4), d_m_per = c(3.7, 0.3, 2.6)
  )
  prepared <- prepare_mfeargen_data(x, "2cs", rating_scale = c(0, 10))
  expect_equal(as.vector(prepared$d_p_per), x$d_p_per)
  expect_equal(as.vector(prepared$d_m_per), x$d_m_per)
  expect_equal(as.vector(prepared$r_minus), x$USm)
  expect_equal(prepared$lambda_threshold, -log(0.7) / 4)
})

test_that("preparation rejects missing experimental information", {
  x <- simulate_mfeargen(1, n_conditioning_trials = 4, n_generalization_trials = 1, seed = 3)
  no_outcomes <- x
  no_outcomes$US <- NULL
  expect_error(prepare_mfeargen_data(no_outcomes), "reinforcement")
  missing_perception <- x
  missing_perception$Per_size[2] <- NA
  expect_error(prepare_mfeargen_data(missing_perception), "finite measurements")
  no_reference <- x
  no_reference$stimulus_type <- "GS"
  expect_error(prepare_mfeargen_data(no_reference), "reference")
  duplicate <- rbind(x, x[1, ])
  expect_error(prepare_mfeargen_data(duplicate), "Duplicate")
  x$US[1] <- NA
  expect_error(prepare_mfeargen_data(x), "reinforcement outcomes")
})

test_that("rating transforms validate bounds and retain missing ratings", {
  expect_equal(normalize_ratings(c(-5, 5, 15, NA), c(-5, 15)), c(0, 0.5, 1, NA))
  expect_equal(denormalize_ratings(c(0, 0.5, 1, NA), c(-5, 15)), c(-5, 5, 15, NA))
  expect_error(normalize_ratings(c(1, 2), c(2, 2)), "increasing")
  expect_error(normalize_ratings(c(1, 11), c(0, 10)), "outside")
  expect_error(denormalize_ratings(2, c(0, 10)), "between")
  expect_error(detect_rating_scale(rep(4, 3)), "explicit rating_scale")
  expect_error(detect_rating_scale(c(NA_real_, NA_real_)), "without observed")
  expect_error(detect_rating_scale(c(1, Inf)), "finite")
  expect_equal(detect_rating_scale(c(NA, 1, 9), buffer = 1), c(0, 10))
})

test_that("factor identifiers with unused levels remain valid", {
  x <- simulate_mfeargen(1, n_conditioning_trials = 4, n_generalization_trials = 1, seed = 7)
  x$participant <- factor("subject", levels = c("absent", "subject"))
  expect_true(validate_mfeargen_data(x, "1cs", verbose = FALSE))
  expect_identical(as.character(attr(prepare_mfeargen_data(x), "participants")), "subject")
})


test_that("prior settings cannot be silently discarded during preparation", {
  data <- create_test_data()
  expect_error(prepare_mfeargen_data(data, priors = set_priors()), "Supply priors to fit")
  expect_error(prepare_mfeargen_data(data, use_flexible_priors = TRUE), "Supply priors to fit")
})

test_that("partial explicit schemas are rejected rather than overwritten", {
  one <- create_test_data()
  one$d_phy <- 999
  expect_error(prepare_mfeargen_data(one), "Incomplete explicit trial columns")
  expect_error(validate_mfeargen_data(one, verbose = FALSE), "Incomplete explicit trial columns")
  two <- create_test_data(model_type = "2cs", cs_plus = 7)
  two$k_plus <- 1
  expect_error(prepare_mfeargen_data(two, "2cs"), "k_minus")
})

test_that("public validation checks the same trial covariates as preparation", {
  data <- create_test_data()
  invalid_labels <- data
  invalid_labels$stimulus_type[1] <- "CS?"
  expect_error(validate_mfeargen_data(invalid_labels, verbose = FALSE), "stimulus_type")
  duplicated_names <- data
  names(duplicated_names)[2] <- "participant"
  expect_error(validate_mfeargen_data(duplicated_names, verbose = FALSE), "distinct")
  empty_id <- data
  empty_id$participant <- "  "
  expect_error(validate_mfeargen_data(empty_id, verbose = FALSE), "identifiers")
  overlapping <- create_test_data(model_type = "2cs", cs_plus = 7)
  overlapping$k_plus <- overlapping$k_minus <- 1
  expect_error(validate_mfeargen_data(overlapping, "2cs", verbose = FALSE), "both be 1")
})

test_that("normal and censored data retain the selected observation model", {
  data <- create_test_data()
  data$US_expect[1] <- -3
  expect_error(prepare_mfeargen_data(data, rating_scale = c(0, 10)), "outside")
  normal <- prepare_mfeargen_data(data, rating_scale = c(0, 10), observation_model = "normal")
  expect_equal(normal$y[1, 1], -3)
  expect_identical(attr(normal, "observation_model"), "normal")
  data$US_expect[1] <- 0
  censored <- prepare_mfeargen_data(data, rating_scale = c(0, 10))
  expect_identical(attr(censored, "observation_model"), "censored")
})

test_that("logical outcome indicators work without coercing the original data", {
  data <- create_test_data()
  data$US <- as.logical(data$US)
  data$k <- data$stimulus_type == "CS+"
  expect_true(validate_mfeargen_data(data, verbose = FALSE))
  prepared <- prepare_mfeargen_data(data)
  expect_type(prepared$r, "double")
  expect_type(prepared$k, "double")
  expect_type(data$US, "logical")
})

test_that("rating scale arithmetic cannot overflow", {
  expect_error(normalize_ratings(0, c(-1e308, 1e308)), "finite")
  expect_error(detect_rating_scale(c(-1e308, 1e308)), "finite")
  expect_error(detect_rating_scale(c(0, 1e308), buffer = 1e308), "finite")
})

test_that("normalization respects extreme but representable response units", {
  data <- create_test_data(n_participants = 1)
  data$US_expect <- seq(-1e308, 0, length.out = nrow(data))
  prepared <- prepare_mfeargen_data(data, rating_scale = c(-1e308, 0), observation_model = "normal")
  expect_true(all(is.finite(prepared$y)))
  expect_equal(range(prepared$y), c(0, 10))
  expect_error(calculate_lambda_threshold(1e-320), "rescale stimulus")
  expect_error(calculate_lambda_threshold(1, stimulus_range = c(1, 2)), "only one")
})


test_that("instrument endpoints retain exact censoring codes after rescaling", {
  data <- create_test_data(n_participants = 1)
  data$US_expect <- rep(c(-0.3, 0.6), length.out = nrow(data))
  prepared <- prepare_mfeargen_data(data, rating_scale = c(-0.3, 0.6),
                                    observation_model = "censored")
  expect_identical(as.vector(prepared$y), rep(c(0, 10), length.out = nrow(data)))
})
