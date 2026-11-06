load_coefficient_reference <- function() {
  read.csv(system.file("extdata", "published_coefficients.csv", package = "GOLDBioAge"),
           stringsAsFactors = FALSE)
}

test_that("the installed reference contains the five existing cohort models", {
  reference <- load_coefficient_reference()
  expect_equal(nrow(reference), 37L)
  expect_setequal(unique(paste(reference$dataset, reference$model)),
                  c("nhanes gold", "ukb gold", "nhanes light", "ukb light", "charls light"))
  expect_false(anyDuplicated(reference[, c("dataset", "model", "term")]) > 0)
  expect_true(all(is.finite(reference$coefficient)))
  for (key in unique(paste(reference$dataset, reference$model))) {
    selected <- reference[paste(reference$dataset, reference$model) == key, ]
    expect_equal(selected$coefficient[selected$term == "age"], 1)
    expect_equal(sum(selected$term == "beta0"), 1L)
  }
})

test_that("the original GOLD coefficients reproduce the worked individual scores", {
  reference <- load_coefficient_reference()
  person <- data.frame(age = 50, albumin = 4, alp = 80, creat = 1,
                       glucose_mmol = 5, lymph = 30, mcv = 90, rdw = 13, wbc = 7, ggt = 25)
  for (dataset in c("nhanes", "ukb")) {
    selected <- reference[reference$dataset == dataset & reference$model == "gold", ]
    beta <- setNames(selected$coefficient, selected$term)
    expected <- if (dataset == "nhanes") 51.10203967 else 50.7472
    expect_equal(predict_gold_bioage(beta, person), expected, tolerance = 1e-10)
  }
  nhanes <- reference[reference$dataset == "nhanes" & reference$model == "gold", ]
  expect_equal(nhanes$coefficient[nhanes$term == "lymph"], -0.20318306)
  expect_equal(nhanes$coefficient[nhanes$term == "ggt"], 0.06075152)
})

test_that("LIGHT units and transforms reproduce all three supplied cohort formulas", {
  reference <- load_coefficient_reference()
  expected <- c(nhanes = 50.439978417530952, ukb = 53.382622562072712,
                charls = 52.173009052401632)
  for (dataset in names(expected)) {
    selected <- reference[reference$dataset == dataset & reference$model == "light", ]
    crp <- selected[selected$term == "log_crp", ]
    expect_equal(crp$input_unit, if (dataset == "nhanes") "mg/dL" else "mg/L")
    expect_equal(crp$transformation, "log1p")
    raw_crp <- if (dataset == "nhanes") 0.3 else 3
    inputs <- data.frame(age = 50, creat = 1, glucose_mmol = 5, log_crp = log1p(raw_crp))
    expect_equal(predict_gold_bioage(setNames(selected$coefficient, selected$term), inputs),
                 unname(expected[dataset]), tolerance = 1e-10)
  }
})
