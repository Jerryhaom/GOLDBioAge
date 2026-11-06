#' Calculate GOLD BioAge using frozen coefficients
#'
#' @param object A result returned by [gold_bioage()], or a finite named numeric
#'   coefficient vector containing predictor names and an intercept named `beta0`.
#' @param newdata Data frame with the named numeric predictor columns, in the
#'   same units and transformations used to derive the coefficients. Survival
#'   follow-up and event status are not required.
#' @return A numeric vector, one value per row of `newdata`, in the original row
#'   order. Rows with missing predictors receive `NA`. No model is refitted and
#'   the training intercept is kept unchanged.
#' @export
#' @examples
#' coefficients <- c(age = 1, marker = 2, beta0 = -4)
#' newdata <- data.frame(age = c(50, 60), marker = c(2, NA))
#' predict_gold_bioage(coefficients, newdata)
predict_gold_bioage <- function(object, newdata) {
  coefficients <- if (is.list(object)) object$coefficients else object
  if (!is.numeric(coefficients) || !is.null(dim(coefficients)) ||
      is.null(names(coefficients)) || anyNA(names(coefficients)) ||
      any(!nzchar(names(coefficients))) || anyDuplicated(names(coefficients)) ||
      !all(is.finite(coefficients)) || !"beta0" %in% names(coefficients) ||
      length(coefficients) < 2L) {
    stop("Supply a fitted model or a finite named coefficient vector including 'beta0'.")
  }
  if (!is.data.frame(newdata)) stop("'newdata' must be a data frame.")
  if (anyDuplicated(names(newdata))) stop("'newdata' must have unique column names.")
  predictors <- setdiff(names(coefficients), "beta0")
  missing <- setdiff(predictors, names(newdata))
  if (length(missing)) stop("Missing predictor columns: ", paste(missing, collapse = ", "))
  x <- newdata[, predictors, drop = FALSE]
  if (!all(vapply(x, is.numeric, logical(1)))) stop("Prediction columns must be numeric.")
  if (any(vapply(x, function(value) any(is.infinite(value)), logical(1)))) {
    stop("Prediction columns must contain finite values or NA, not Inf.")
  }
  scores <- as.vector(as.matrix(x) %*% coefficients[predictors] + coefficients[["beta0"]])
  scores[is.nan(scores)] <- NA_real_
  scores
}

#' Calculate Cox biological age using frozen coefficients
#'
#' @param object A result from [cox_bioage()], or a finite named numeric vector
#'   of score coefficients including an intercept named `beta0`. A fitted object
#'   saved with version 0.2.0 is also accepted, retaining its legacy intercept.
#' @param newdata Data frame with numeric predictor columns in the training units.
#'   Survival follow-up and event status are not needed.
#' @return A numeric vector in the input row order, with `NA` for rows missing
#'   predictors. The fitted intercept is retained without recalibration.
#' @details
#' The legacy matrix `object$coefficients` contains weights only. This function
#' also uses `object$intercept`; never omit that intercept when applying the score.
#' The equivalent frozen vector is available as `object$score_coefficients` in
#' version 0.2.1 and later. A Cox score is not a mortality probability.
#' @export
#' @examples
#' frozen <- c(age = 1, marker = 2, beta0 = -4)
#' predict_cox_bioage(frozen, data.frame(age = c(50, 60), marker = c(2, NA)))
predict_cox_bioage <- function(object, newdata) {
  if (is.list(object)) {
    weights <- object$coefficients
    intercept <- object$intercept
    if (!is.matrix(weights) || !is.numeric(weights) || ncol(weights) != 1L || is.null(rownames(weights)) ||
        !is.numeric(intercept) || length(intercept) != 1L || !is.finite(intercept)) {
      stop("Supply a fitted Cox model with a coefficient matrix and a finite intercept.")
    }
    coefficients <- c(stats::setNames(as.numeric(weights), rownames(weights)), beta0 = intercept)
  } else {
    coefficients <- object
  }
  predict_gold_bioage(coefficients, newdata)
}
