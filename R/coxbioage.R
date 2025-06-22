#' Calculate biological age using Cox regression model
#'
#' @description
#' Computes biological age based on two Cox regression models:
#' \itemize{
#'   \item Cox Model 1: Surv(time, status) ~ age
#'   \item Cox Model 2: Surv(time, status) ~ age + Biomarkers
#' }
#'
#' @param d4 Dataframe containing survival information (time, status)
#' @param var Character vector of variable names (must include 'age')
#'
#' @return A list containing:
#' \itemize{
#'   \item model1_coef: Coefficients from Model 1
#'   \item model2_coef: Coefficients from Model 2
#'   \item coefficients: Biological age coefficients
#'   \item Bioage: Dataframe with calculated biological age
#' }
#'
#' @importFrom stats as.formula na.omit
#' @importFrom survival Surv
#' @export
#'
#' @examples
#' head(NHANES4)
#' var <- c(
#'   "age", "albumin", "alp", "creat",
#'   "glucose_mmol", "lymph",
#'   "mcv", "rdw", "wbc", "ggt"
#' )
#'
#' Cox.bioage <- cox_bioage(NHANES4, var)

cox_bioage <- function(d4, var) {
  data <- data.frame(
    time = d4$time,
    status = d4$status,
    d4[, var]
  )
  data <- na.omit(data)

  data1 <- data
  fitg1 <- survival::coxph(survival::Surv(time, status) ~ age, data = data1)
  coef1 <- coef(fitg1)

  ############
  coef_fixed <- coef1[1]
  formula_str <- paste0("survival::Surv(time, status) ~ offset(", coef_fixed, " * age) + ",
                        paste(var[-1], collapse = " + "))
  fit_fixed <- survival::coxph(as.formula(formula_str), data = data1)
  coef2 = c(coef_fixed, coef(fit_fixed))

  risk <- exp(predict(fitg1, type = "lp"))
  risk2 <- exp(predict(fit_fixed, type = "lp"))

  coef21=matrix(coef2[1:length(coef2)])/coef1[1]
  rownames(coef21)=var
  coef21

  x <- as.matrix(data[, var]) %*% as.matrix(coef21)
  la<-risk/risk2

  x0=mean(log(la)/coef1[1]) #(coef2[1]-coef1[1])
  bioage <- x0 + x
  data$Coxbioage <- bioage

  list(
    `model1 coef` = coef1,
    `model2 coef` = coef2,
    `coefficients` = coef21,
    `CoxBioage` = data
  )
}
