#' Fit a Cox regression biological age model
#'
#' @description
#' Fits an age-only Cox model and a biomarker model with the age coefficient
#' fixed by an offset. This is a Cox-based alternative to GOLD BioAge.
#' The default intercept centers the training mean score on chronological age.
#' @param d4 Data frame containing numeric model columns, non-negative `time`,
#'   and `status` coded 0 (censored) or 1 (event).
#' @param var Unique predictor column names including `age`, in any order.
#' @param calibration Intercept convention: `"centered"` (default) estimates
#'   `mean(age - linear_score)` on complete training rows. `"legacy"` reproduces
#'   the centered-prediction intercept used through version 0.2.0, whose scores
#'   do not generally have the same training mean as chronological age.
#' @details
#' The final Cox model uses Efron handling of tied events and fixes the age
#' coefficient to the age-only estimate using an offset. Biomarker coefficients
#' are divided by that reference age coefficient. The default intercept is
#' estimated once on training data and is retained for new observations.
#' This is a relative-risk-derived age score, not an absolute mortality
#' probability or proof that two fitted Cox baseline hazards are equal.
#' Use [predict_cox_bioage()] for frozen scoring. The two fitted Cox models
#' are also returned for inspecting their coefficients and fitted model objects.
#' @return A list with `"model1 coef"`, `"model2 coef"`, `coefficients` (the
#'   original coefficient matrix), `intercept`, `CoxBioage` (complete-case data
#'   with numeric `Coxbioage`), `rows_used`, and `rows_omitted`. Additional fields
#'   are `score_coefficients` (named weights including `beta0`), `calibration`,
#'   and `models` (age-only and constrained full Cox models).
#' @export
#' @examples
#' \donttest{
#' data(NHANES4)
#' variables <- c("age", "albumin", "creat", "glucose_mmol")
#' fit <- cox_bioage(NHANES4, variables)
#' fit$coefficients
#' fit$intercept
#' scores <- predict_cox_bioage(fit, NHANES4)
#' }
cox_bioage <- function(d4, var, calibration = c("centered", "legacy")) {
  calibration <- match.arg(calibration)
  prepared <- .prepare_bioage_data(d4, var)
  data <- prepared$data
  var <- prepared$var
  fitg1 <- survival::coxph(.bioage_formula("age"), data = data,
                           ties = "efron", model = TRUE, x = TRUE)
  coef1 <- stats::coef(fitg1)
  age_effect <- unname(coef1["age"])
  if (!is.finite(age_effect) || abs(age_effect) < sqrt(.Machine$double.eps)) {
    stop("The age-only age coefficient is zero or non-finite; cannot express effects in years.")
  }
  if (length(var) == 1L) {
    fit_fixed <- fitg1
    coef2 <- coef1
  } else {
    fit_fixed <- survival::coxph(.bioage_formula(setdiff(var, "age"), offset_age = age_effect),
                                 data = data, ties = "efron", model = TRUE, x = TRUE)
    coef2 <- c(age = age_effect, stats::coef(fit_fixed))
  }
  if (any(!is.finite(coef2))) {
    stop("The Cox model has non-finite coefficients; check collinearity and convergence.")
  }
  coefficients <- matrix(coef2[var] / age_effect, ncol = 1L, dimnames = list(var, NULL))
  linear_score <- as.vector(as.matrix(data[, var, drop = FALSE]) %*% coefficients)
  intercept <- if (calibration == "centered") {
    mean(data$age - linear_score)
  } else {
    mean(stats::predict(fitg1, type = "lp") -
           stats::predict(fit_fixed, type = "lp")) / age_effect
  }
  data$Coxbioage <- as.numeric(linear_score + intercept)
  list(`model1 coef` = coef1, `model2 coef` = coef2,
       coefficients = coefficients, intercept = intercept, CoxBioage = data,
       rows_used = prepared$rows_used, rows_omitted = prepared$rows_omitted,
       score_coefficients = c(stats::setNames(as.numeric(coefficients), rownames(coefficients)),
                              beta0 = intercept), calibration = calibration,
       models = list(age_only = fitg1, full = fit_fixed))
}
