#' Fit a GOLD BioAge model using Gompertz regression
#'
#' @description
#' Fits an age-only Gompertz survival model, then an age-plus-biomarker model
#' with its age coefficient fixed to the age-only estimate. Biomarker effects
#' are expressed in years, and the intercept centers mean GOLD BioAge on mean
#' chronological age in the complete-case training sample.
#'
#' @param d4 Data frame with numeric predictors, non-negative follow-up `time`,
#'   and `status` coded 0 (censored) or 1 (event). Age must be measured in years.
#' @param var Unique predictor column names including `age`, in any order.
#'   Names refer to columns, not formula expressions. Rows missing any model
#'   variable are omitted and their original positions are returned.
#' @param feature_selection Apply optional Cox penalized regression to the
#'   candidate biomarkers before the final Gompertz fit. Age is always retained.
#' @param selection_method `"lasso"`, `"elasticnet"`, or `"none"`. The `"none"`
#'   option disables selection even if `feature_selection = TRUE`.
#' @param alpha Elastic-net mixing parameter between 0 and 1; ignored for LASSO.
#' @param nfolds Number of cross-validation folds, at least 3.
#' @param family Must be `"cox"`; retained for compatibility.
#' @param seed Seed for cross-validation (default 123). The caller's random
#'   number state is restored afterwards. Use `NULL` to use the caller's state.
#'
#' @details
#' Selection uses the biomarker-only Cox model at `lambda.min`; age is added
#' to the final Gompertz model. This screening step does not constitute external
#' validation. It requires the optional packages `glmnet` and `ggplot2`, strictly
#' positive follow-up times, and at least `nfolds` events. Training and selection
#' must be restricted to the development sample.
#' Breslow handling of tied events is retained across supported `glmnet` versions.
#'
#' This function refits a model on `d4`; it does not automatically apply the
#' published NHANES or UKB coefficients. Use [predict_gold_bioage()] with frozen
#' fitted coefficients or a named published coefficient vector to score new data.
#'
#' @return A list retaining the existing fields `"model1 coef"`, `"model2 coef"`,
#'   `coefficients` (a named numeric vector including `beta0`), and `GOLDBioage`
#'   (complete-case training data with a numeric `GOLDbioage` column). Additional
#'   fields are `selected_vars`, `original_vars`, `rows_used`, `rows_omitted`,
#'   `models` (the two fitted Gompertz models), and `cv_plot` (`NULL` without
#'   selection). The selected model also returns `cv_fit` for diagnostics.
#' @export
#' @examples
#' \donttest{
#' data(NHANES4)
#' variables <- c("age", "albumin", "alp", "creat", "glucose_mmol",
#'                "lymph", "mcv", "rdw", "wbc", "ggt")
#' fit <- gold_bioage(NHANES4, variables)
#' fit$coefficients
#' scores <- predict_gold_bioage(fit, NHANES4)
#' }
gold_bioage <- function(d4, var, feature_selection = FALSE,
                       selection_method = "lasso", alpha = 1,
                       nfolds = 10, family = "cox", seed = 123) {
  original_vars <- var
  prepared <- .prepare_bioage_data(d4, var)
  data <- prepared$data
  var <- prepared$var
  if (!is.logical(feature_selection) || length(feature_selection) != 1L ||
      is.na(feature_selection)) stop("'feature_selection' must be TRUE or FALSE.")
  selection_method <- match.arg(selection_method, c("lasso", "elasticnet", "none"))
  if (!identical(family, "cox")) stop("'family' must be 'cox'.")
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
                        !is.finite(seed) || seed < 0 || seed > .Machine$integer.max ||
                        seed != as.integer(seed))) {
    stop("'seed' must be NULL or a non-negative integer within the R integer range.")
  }

  biomarker_vars <- setdiff(var, "age")
  selected_vars <- var
  cv_fit <- NULL
  cv_plot <- NULL
  select <- feature_selection && selection_method != "none" && length(biomarker_vars) > 1L
  if (select) {
    for (pkg in c("glmnet", "ggplot2")) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Feature selection requires the optional package '", pkg, "'.")
      }
    }
    if (selection_method == "elasticnet" &&
        (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha > 1)) {
      stop("'alpha' must be between 0 and 1.")
    }
    if (!is.numeric(nfolds) || length(nfolds) != 1L || !is.finite(nfolds) ||
        nfolds < 3L || nfolds > nrow(data) || nfolds != as.integer(nfolds)) {
      stop("'nfolds' must be an integer between 3 and the number of complete rows.")
    }
    if (any(data$time <= 0)) stop("Cox feature selection requires strictly positive 'time'.")
    if (sum(data$status == 1) < nfolds) stop("Feature selection requires at least 'nfolds' events.")
    cv_fit <- .with_bioage_seed(seed, function() {
      args <- list(x = as.matrix(data[, biomarker_vars, drop = FALSE]),
                   y = survival::Surv(data$time, data$status), family = "cox",
                   alpha = if (selection_method == "lasso") 1 else alpha, nfolds = nfolds)
      # Preserve the historical Cox tie method as newer glmnet defaults change.
      if (utils::packageVersion("glmnet") >= "5.0") args$cox.ties <- "breslow"
      do.call(glmnet::cv.glmnet, args)
    })
    # Cox glmnet has no intercept row: retain every named biomarker coefficient.
    selected <- as.matrix(stats::coef(cv_fit, s = "lambda.min"))[, 1]
    selected_vars <- c("age", names(selected)[selected != 0])
    cv_plot <- ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(x = log(cv_fit$lambda), y = cv_fit$cvm)) +
      ggplot2::geom_vline(xintercept = log(cv_fit$lambda.min), linetype = "dashed", color = "red") +
      ggplot2::labs(title = paste("Cross-validation for", selection_method),
                    x = "Log(Lambda)", y = "Partial likelihood deviance") +
      ggplot2::theme_minimal()
    message("Selected ", length(selected_vars) - 1L, " of ", length(biomarker_vars), " biomarkers.")
  } else if (feature_selection && selection_method != "none" && length(biomarker_vars) <= 1L) {
    warning("Not enough biomarkers for feature selection. Using all variables.")
  }

  fitg1 <- flexsurv::flexsurvreg(.bioage_formula("age"), data = data, dist = "gompertz")
  coef1 <- fitg1$res[, "est"]
  age_effect <- unname(coef1["age"])
  if (!is.finite(age_effect) || abs(age_effect) < sqrt(.Machine$double.eps)) {
    stop("The age-only age coefficient is zero or non-finite; cannot express effects in years.")
  }
  if (length(selected_vars) == 1L) {
    fitg2 <- fitg1
  } else {
    formula <- .bioage_formula(selected_vars)
    initial_fit <- flexsurv::flexsurvreg(formula, data = data, dist = "gompertz")
    initial <- initial_fit$res[, "est"]
    age_index <- match("age", names(initial))
    # Fix chronological age, not the Gompertz shape parameter.
    initial[age_index] <- age_effect
    fitg2 <- flexsurv::flexsurvreg(formula, data = data, dist = "gompertz",
                                  fixedpars = age_index, inits = initial)
  }
  coef2 <- fitg2$res[, "est"]
  weights <- coef2[selected_vars] / age_effect
  linear_score <- as.vector(as.matrix(data[, selected_vars, drop = FALSE]) %*% weights)
  # Algebraically identical to the original mean log-hazard-ratio calibration,
  # computed on the linear scale to avoid exponentiation and overflow.
  intercept <- mean(data$age - linear_score)
  coefficients <- c(weights, beta0 = intercept)
  data$GOLDbioage <- as.numeric(linear_score + intercept)

  list(`model1 coef` = coef1, `model2 coef` = coef2,
       coefficients = coefficients, GOLDBioage = data,
       selected_vars = selected_vars, original_vars = original_vars,
       rows_used = prepared$rows_used, rows_omitted = prepared$rows_omitted,
       models = list(age_only = fitg1, full = fitg2), cv_plot = cv_plot, cv_fit = cv_fit)
}
