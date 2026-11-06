#' Fit a Cox regression biological age model
#'
#' @description
#' Fits an age-only Cox model and a biomarker model with the age coefficient
#' fixed by an offset. This is a Cox-based alternative to GOLD BioAge.
#' @param d4 Data frame containing numeric model columns, non-negative `time`,
#'   and `status` coded 0 (censored) or 1 (event).
#' @param var Unique predictor column names including `age`, in any order.
#' @return A list with `"model1 coef"`, `"model2 coef"`, `coefficients` (the
#'   original coefficient matrix), `intercept`, `CoxBioage` (complete-case data
#'   with numeric `Coxbioage`), `rows_used`, and `rows_omitted`.
#' @export
#' @examples
#' \donttest{
#' data(NHANES4)
#' variables <- c("age", "albumin", "creat", "glucose_mmol")
#' fit <- cox_bioage(NHANES4, variables)
#' fit$coefficients
#' }
cox_bioage <- function(d4, var) {
  prepared <- .prepare_bioage_data(d4, var)
  data <- prepared$data
  var <- prepared$var
  fitg1 <- survival::coxph(.bioage_formula("age"), data = data)
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
                                 data = data)
    coef2 <- c(age = age_effect, stats::coef(fit_fixed))
  }
  coefficients <- matrix(coef2[var] / age_effect, ncol = 1L, dimnames = list(var, NULL))
  linear_score <- as.vector(as.matrix(data[, var, drop = FALSE]) %*% coefficients)
  intercept <- mean(stats::predict(fitg1, type = "lp") -
                    stats::predict(fit_fixed, type = "lp")) / age_effect
  data$Coxbioage <- as.numeric(linear_score + intercept)
  list(`model1 coef` = coef1, `model2 coef` = coef2,
       coefficients = coefficients, intercept = intercept, CoxBioage = data,
       rows_used = prepared$rows_used, rows_omitted = prepared$rows_omitted)
}
