test_that("Cox scoring centers training age and retains the legacy weight matrix", {
  fit <- cox_bioage(test_data, test_vars)
  expect_equal(fit$calibration, "centered")
  expect_equal(fit$coefficients[, 1], c(age = 1, marker_risk = 5.0637988452291,
                                       marker_protective = -3.57237221985725,
                                       noise = -0.616588835401311), tolerance = 1e-6)
  expect_equal(mean(fit$CoxBioage$Coxbioage - test_data$age), 0, tolerance = 1e-10)
  expect_equal(predict_cox_bioage(fit, test_data), fit$CoxBioage$Coxbioage)
  expect_equal(predict_cox_bioage(fit$score_coefficients, test_data), fit$CoxBioage$Coxbioage)
  expect_equal(fit$score_coefficients[["beta0"]], fit$intercept)
  expect_s3_class(fit$models$age_only, "coxph")
  expect_s3_class(fit$models$full, "coxph")
})

test_that("legacy calibration reproduces the previous intercept and scores", {
  fit <- cox_bioage(test_data, test_vars, calibration = "legacy")
  expect_equal(fit$intercept, -60.1767136885924, tolerance = 1e-10)
  old_weights <- c(age = 1, marker_risk = 5.0637988452291,
                   marker_protective = -3.57237221985725, noise = -0.616588835401311)
  old_scores <- as.vector(as.matrix(test_data[, test_vars]) %*% old_weights - 60.1767136885924)
  expect_equal(fit$CoxBioage$Coxbioage, old_scores, tolerance = 1e-8)
  # A saved 0.2.0 object had neither score_coefficients nor stored Cox models.
  old_object <- list(coefficients = fit$coefficients, intercept = fit$intercept)
  expect_equal(predict_cox_bioage(old_object, test_data), old_scores, tolerance = 1e-8)
})

test_that("centered Cox scores are invariant to translating biomarker values", {
  shifted <- test_data
  shifted$marker_risk <- shifted$marker_risk + 10
  a <- cox_bioage(test_data, test_vars)
  b <- cox_bioage(shifted, test_vars)
  expect_equal(a$coefficients, b$coefficients, tolerance = 1e-8)
  expect_equal(a$CoxBioage$Coxbioage, b$CoxBioage$Coxbioage, tolerance = 1e-8)
})

test_that("Cox prediction preserves missing rows and the frozen intercept", {
  d <- test_data
  d$marker_risk[2] <- NA_real_
  d$time[4] <- NA_real_
  expect_message(fit <- cox_bioage(d, test_vars), "Omitted 2 rows")
  scores <- predict_cox_bioage(fit, d)
  expect_length(scores, nrow(d))
  expect_true(is.na(scores[2]))
  expect_true(is.finite(scores[4]))
  expect_equal(scores[fit$rows_used], fit$CoxBioage$Coxbioage)
  expect_equal(predict_cox_bioage(fit, d[100:90, ]), scores[100:90])
  expect_error(predict_cox_bioage(list(coefficients = fit$coefficients), d), "finite intercept")
  expect_error(predict_cox_bioage(fit, d["age"]), "Missing predictor")
})
