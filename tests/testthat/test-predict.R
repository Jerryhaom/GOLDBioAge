test_that("published NHANES coefficients match the calculator worked example", {
  coefficients <- c(age = 1, albumin = -4.73575224, alp = 0.02604091,
                    creat = 5.26909180, glucose_mmol = 0.57969768, lymph = -0.20318306,
                    mcv = 0.33888497, rdw = 2.64452031, wbc = 0.44588830,
                    ggt = 0.06075152, beta0 = -53.62873)
  person <- data.frame(age = 50, albumin = 4, alp = 80, creat = 1,
                       glucose_mmol = 5, lymph = 30, mcv = 90, rdw = 13, wbc = 7, ggt = 25)
  expect_equal(predict_gold_bioage(coefficients, person), 51.10203967, tolerance = 1e-10)
  a <- person
  a$lymph <- 40
  b <- person
  b$ggt <- 35
  expect_equal(predict_gold_bioage(coefficients, a) - predict_gold_bioage(coefficients, person),
               -2.0318306, tolerance = 1e-10)
  expect_equal(predict_gold_bioage(coefficients, b) - predict_gold_bioage(coefficients, person),
               0.6075152, tolerance = 1e-10)
})

test_that("prediction keeps the intercept fixed and uses coefficient names", {
  coefficients <- c(age = 1, marker = 2, beta0 = -4)
  d <- data.frame(marker = c(2, NA, 5), age = c(50, 60, 70), ignored = "unused")
  expect_equal(predict_gold_bioage(coefficients, d), c(50, NA_real_, 76))
  expect_equal(predict_gold_bioage(coefficients, d[3:1, ]), c(76, NA_real_, 50))
  expect_equal(predict_gold_bioage(coefficients, d[FALSE, ]), numeric())
  expect_error(predict_gold_bioage(coefficients, d["age"]), "Missing predictor")
  expect_error(predict_gold_bioage(c(age = 1), d), "beta0")
  expect_error(predict_gold_bioage(c(age = 1, beta0 = Inf), d), "finite")
  expect_error(predict_gold_bioage(c(age = 1, age = 2, beta0 = 0), d), "named coefficient")
  expect_error(predict_gold_bioage(coefficients, transform(d, marker = "bad")), "numeric")
  expect_error(predict_gold_bioage(coefficients, transform(d, marker = Inf)), "not Inf")
})
