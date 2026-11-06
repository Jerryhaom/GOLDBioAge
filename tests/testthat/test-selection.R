test_that("Cox screening includes all named nonzero coefficients, including the first biomarker", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("ggplot2")
  set.seed(42)
  rng_before <- .Random.seed
  expect_message(fit <- gold_bioage(test_data, test_vars, feature_selection = TRUE,
                                    nfolds = 5, seed = 123), "Selected")
  expect_identical(.Random.seed, rng_before)
  weights <- as.matrix(stats::coef(fit$cv_fit, s = "lambda.min"))[, 1]
  expect_setequal(fit$selected_vars, c("age", names(weights)[weights != 0]))
  expect_true("marker_risk" %in% fit$selected_vars)
  expect_s3_class(fit$cv_plot, "ggplot")
  expect_true(all(is.finite(predict_gold_bioage(fit, test_data))))
  expect_message(repeated <- gold_bioage(test_data, test_vars, feature_selection = TRUE,
                                         nfolds = 5, seed = 123), "Selected")
  expect_equal(repeated$cv_fit$lambda.min, fit$cv_fit$lambda.min)
  expect_equal(repeated$coefficients, fit$coefficients)
})

test_that("elastic net and explicit no-selection are supported", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("ggplot2")
  expect_message(fit <- gold_bioage(test_data, test_vars, feature_selection = TRUE,
                                    selection_method = "elasticnet", alpha = 0.5,
                                    nfolds = 5), "Selected")
  expect_equal(fit$coefficients[["age"]], 1)
  none <- gold_bioage(test_data, test_vars, feature_selection = TRUE, selection_method = "none")
  expect_equal(none$selected_vars, test_vars)
  expect_null(none$cv_fit)
  expect_null(none$cv_plot)
})

test_that("invalid screening settings are rejected before cross-validation", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("ggplot2")
  expect_error(gold_bioage(test_data, test_vars, feature_selection = TRUE,
                          selection_method = "elasticnet", alpha = 2), "alpha")
  expect_error(gold_bioage(test_data, test_vars, feature_selection = TRUE,
                          nfolds = 2), "nfolds")
  expect_error(gold_bioage(transform(test_data, time = 0), test_vars,
                          feature_selection = TRUE), "strictly positive")
  d <- test_data
  d$status <- 0
  d$status[1:4] <- 1
  expect_error(gold_bioage(d, test_vars, feature_selection = TRUE,
                          nfolds = 5), "at least 'nfolds' events")
})
