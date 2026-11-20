test_that("group predictions use posterior samples and retain IDs", {
  object <- structure(list(
    jags_fit = list(sims.list = list(gp = cbind(c(3, 3, 3, 1), c(1, 2, 3, 4)))),
    jags_data = structure(list(), participants = c("subject-9", "subject-2"))
  ), class = "mfeargen")
  groups <- predict_groups(object, threshold = 0.6)
  expect_identical(groups$participant, c("subject-9", "subject-2"))
  expect_equal(groups$max_probability, c(0.75, 0.25))
  expect_equal(groups$group, c(3L, NA_integer_))
  expect_identical(groups$group_label, c("Physical", NA_character_))
  expect_identical(groups$assigned_group, groups$group_label)
  expect_equal(predict_groups(object, threshold = 0)$group, c(3L, 1L))
  expect_error(predict_groups(object, threshold = 2), "threshold")
})
