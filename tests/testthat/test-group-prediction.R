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

test_that("hierarchical group prediction uses the same posterior probabilities", {
  object <- structure(list(
    jags_fit = list(sims.list = list(gp = cbind(c(4, 4, 4, 1), c(2, 2, 2, 2)))),
    jags_data = structure(list(), participants = c("A", "B"))
  ), class = "mfeargen_hierarchical")
  expect_identical(predict_groups(object)$group_label, c("Perceptual", "Overgeneralizer"))
  expect_identical(predict_groups(object)$participant, c("A", "B"))
  object$group_membership <- data.frame(prob_group_1 = 0.5, prob_group_2 = 0.5,
                                       prob_group_3 = 0.5, prob_group_4 = 0.5)
  expect_error(predict_groups(object), "sum to one")
  object$group_membership$prob_group_1 <- NA_real_
  expect_error(predict_groups(object), "finite")
})
