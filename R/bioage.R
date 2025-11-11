#' Calculate biological age using Gompertz model with feature selection
#'
#' @description
#' Computes biological age based on two Gompertz regression models with optional
#' Cox-LASSO feature selection for biomarkers
#'
#' @param d4 Dataframe containing survival information (time, status)
#' @param var Character vector of variable names (must include 'age')
#' @param feature_selection Logical indicating whether to perform feature selection (default: FALSE)
#' @param selection_method Method for feature selection: "lasso", "elasticnet", or "none"
#' @param alpha Elasticnet mixing parameter (0 = ridge, 1 = lasso)
#' @param nfolds Number of folds for cross-validation
#' @param family Survival family for glmnet (default: "cox")
#'
#' @return A list containing:
#' \itemize{
#'   \item model1_coef: Coefficients from Model 1
#'   \item model2_coef: Coefficients from Model 2
#'   \item coefficients: Biological age coefficients
#'   \item Bioage: Dataframe with calculated biological age
#'   \item selected_vars: Variables selected by feature selection (if performed)
#'   \item cv_plot: Cross-validation plot (if feature selection performed)
#' }
#'
#' @importFrom stats as.formula na.omit predict
#' @importFrom flexsurv flexsurvreg
#' @importFrom survival Surv
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom ggplot2 ggplot aes geom_point geom_vline labs theme_minimal
#' @export
#'
#' @examples
#' # Load dataset
#' data(NHANES4)
#' 
#' # Define variables
#' var <- c("age", "albumin", "alp", "creat", "glucose_mmol", "lymph", "mcv", "rdw", "wbc", "ggt")
#'
#' # Without feature selection
#' result1 <- gold_bioage(NHANES4, var)
#'
#' # With LASSO feature selection
#' result2 <- gold_bioage(NHANES4, var, 
#'                       feature_selection = TRUE,
#'                       selection_method = "lasso")
#'
#' # With Elasticnet feature selection
#' result3 <- gold_bioage(NHANES4, var, 
#'                       feature_selection = TRUE,
#'                       selection_method = "elasticnet",
#'                       alpha = 0.5)

gold_bioage <- function(d4, var, feature_selection = FALSE, 
                       selection_method = "lasso", alpha = 1, 
                       nfolds = 10, family = "cox") {
  
  # Check if required packages are available
  if (feature_selection && !requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required for feature selection but not installed.")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }
  
  # Validate inputs
  if (!"age" %in% var) {
    stop("'age' must be included in the variable list")
  }
  
  if (!all(c("time", "status") %in% names(d4))) {
    stop("Data must contain 'time' and 'status' columns")
  }
  
  # Prepare data
  data <- data.frame(
    time = d4$time,
    status = d4$status,
    d4[, var]
  )
  data <- na.omit(data)
  
  # Store original variables
  original_vars <- var
  biomarker_vars <- setdiff(var, "age")
  
  # Feature selection if requested
  selected_vars <- NULL
  cv_plot <- NULL
  
  if (feature_selection && length(biomarker_vars) > 1) {
    
    # Prepare data for glmnet
    x_data <- as.matrix(data[, biomarker_vars])
    y_data <- survival::Surv(data$time, data$status)
    
    # Set alpha based on selection method
    if (selection_method == "lasso") {
      alpha_val <- 1
    } else if (selection_method == "elasticnet") {
      alpha_val <- alpha
    } else {
      alpha_val <- 1  # default to lasso
    }
    
    # Perform cross-validated LASSO
    set.seed(123)  # for reproducibility
    cv_fit <- glmnet::cv.glmnet(x = x_data, y = y_data, 
                               family = family, 
                               alpha = alpha_val,
                               nfolds = nfolds)
    
    # Get selected variables (non-zero coefficients at lambda.1se)
    coefs <- as.matrix(coef(cv_fit, s = "lambda.1se"))
    selected_biomarkers <- biomarker_vars[coefs[-1, 1] != 0]  # exclude intercept
    
    # Always include age and selected biomarkers
    selected_vars <- c("age", selected_biomarkers)
    
    # Create cross-validation plot
    cv_plot <- ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(x = log(cv_fit$lambda), y = cv_fit$cvm)) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = log(cv_fit$lambda.1se)), 
                         linetype = "dashed", color = "red") +
      ggplot2::labs(title = paste("Cross-Validation for", selection_method),
                   x = "Log(Lambda)", 
                   y = "Partial Likelihood Deviance") +
      ggplot2::theme_minimal()
    
    message(paste("Selected", length(selected_biomarkers), 
                  "biomarkers out of", length(biomarker_vars), "using", selection_method))
    
  } else {
    # No feature selection or not enough biomarkers
    selected_vars <- var
    if (feature_selection && length(biomarker_vars) <= 1) {
      warning("Not enough biomarkers for feature selection. Using all variables.")
    }
  }
  
  data1 <- data
  
  # Model 1: Age only
  fitg1 <- flexsurv::flexsurvreg(
    formula = survival::Surv(time, status) ~ age,
    data = data1, dist = "gompertz"
  )
  coef1 <- fitg1$res[, 1]
  
  # Model 2: Age + (selected) biomarkers
  if (length(selected_vars) == 1) {
    # Only age selected, use Model 1 coefficients
    fitg2 <- fitg1
    coef2 <- coef1
  } else {
    fitg2 <- flexsurv::flexsurvreg(
      formula = as.formula(paste(
        "survival::Surv(time, status) ~",
        paste(selected_vars, collapse = "+")
      )),
      data = data1, dist = "gompertz"
    )
    
    # Fix shape parameter as in original function
    ce <- fitg2$res[, 1]
    ce[3] <- coef1[3]  # Use shape parameter from Model 1
    
    fitg2 <- flexsurv::flexsurvreg(
      formula = as.formula(paste(
        "survival::Surv(time, status) ~",
        paste(selected_vars, collapse = "+")
      )),
      fixedpars = 3,
      inits = ce,
      data = data1, dist = "gompertz"
    )
    
    coef2 <- fitg2$res[, 1]
  }
  
  # Calculate biological age coefficients
  coef21 <- matrix(coef2[3:length(coef2)]) / coef1[3]
  rownames(coef21) <- selected_vars
  
  # Calculate biological age
  x <- as.matrix(data[, selected_vars]) %*% as.matrix(coef21)
  pr <- predict(fitg1, type = "hazard", times = c(0))
  pr2 <- predict(fitg2, type = "hazard", times = c(0))
  risk <- pr$.pred_hazard
  risk2 <- pr2$.pred_hazard
  la <- risk / risk2
  x0 <- (log(coef2[2] / coef1[2])) / coef1[3] + mean(log(la) / coef1[3])
  bioage <- x0 + x
  
  data$GOLDbioage <- bioage
  names(coef21) <- selected_vars
  coef21 <- c(coef21, beta0 = as.numeric(x0))
  
  # Prepare results
  result <- list(
    `model1 coef` = coef1,
    `model2 coef` = coef2,
    `coefficients` = coef21,
    `GOLDBioage` = data
  )
  
  # Add feature selection results if performed
  if (feature_selection) {
    result$selected_vars <- selected_vars
    result$cv_plot <- cv_plot
    result$original_vars <- original_vars
  }
  
  return(result)
}
