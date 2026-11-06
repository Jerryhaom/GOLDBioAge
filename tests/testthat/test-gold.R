test_that("GOLD retains the previous no-selection numerical results", {
  fit <- gold_bioage(test_data, test_vars)
  # Frozen reference from the previous release, independently run before editing.
  expected <- c(age = 1, marker_risk = 5.0585766883775,
                marker_protective = -3.57927018071103,
                noise = -0.64797201006524, beta0 = -0.0838610023770074)
  expect_equal(fit$coefficients, expected, tolerance = 1e-6)
  expect_type(fit$GOLDBioage$GOLDbioage, "double")
  expect_null(dim(fit$GOLDBioage$GOLDbioage))
  expect_equal(predict_gold_bioage(fit, test_data), fit$GOLDBioage$GOLDbioage)
  expect_equal(mean(fit$GOLDBioage$GOLDbioage - test_data$age), 0, tolerance = 1e-10)
  expect_equal(fit[["model2 coef"]][["age"]], fit[["model1 coef"]][["age"]])
})

test_that("age can occur anywhere without changing the constrained parameter", {
  a <- gold_bioage(test_data, test_vars)
  reordered <- c("marker_risk", "age", "marker_protective", "noise")
  b <- gold_bioage(test_data, reordered)
  expect_equal(a$coefficients, b$coefficients, tolerance = 1e-8)
  expect_equal(b$original_vars, reordered)
  expect_equal(b$coefficients[["age"]], 1)
  expect_equal(b[["model2 coef"]][["age"]], b[["model1 coef"]][["age"]])
})

test_that("age-only GOLD returns chronological age", {
  fit <- gold_bioage(test_data, "age")
  expect_equal(fit$coefficients, c(age = 1, beta0 = 0))
  expect_equal(fit$GOLDBioage$GOLDbioage, test_data$age)
})

test_that("missing training rows are reported and prediction preserves alignment", {
  d <- test_data
  d$marker_risk[2] <- NA_real_
  d$time[4] <- NA_real_
  expect_message(fit <- gold_bioage(d, test_vars), "Omitted 2 rows")
  expect_equal(fit$rows_omitted, c(2L, 4L))
  expect_equal(fit$rows_used, setdiff(seq_len(nrow(d)), c(2L, 4L)))
  scores <- predict_gold_bioage(fit, d)
  expect_length(scores, nrow(d))
  expect_true(is.na(scores[2]))
  expect_true(is.finite(scores[4])) # Follow-up is not needed for scoring.
  expect_equal(scores[fit$rows_used], fit$GOLDBioage$GOLDbioage)
})

test_that("bundled NHANES refit reproduces the verified coefficients and counts", {
  data("NHANES4", package = "GOLDBioAge")
  variables <- c("age", "albumin", "alp", "creat", "glucose_mmol",
                 "lymph", "mcv", "rdw", "wbc", "ggt")
  expect_message(fit <- gold_bioage(NHANES4, variables), "Omitted 15597 rows")
  expect_equal(nrow(fit$GOLDBioage), 39348L)
  expect_equal(sum(fit$GOLDBioage$status), 4716)
  expected <- c(age = 1, albumin = -4.73575224, alp = 0.02604091,
                creat = 5.26909180, glucose_mmol = 0.57969768, lymph = -0.20318306,
                mcv = 0.33888497, rdw = 2.64452031, wbc = 0.44588830,
                ggt = 0.06075152, beta0 = -53.62873289)
  expect_equal(fit$coefficients, expected, tolerance = 1e-6)
})
