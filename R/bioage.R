#' gold_bioage
#' @description Develop biological age through Gompertz hazard function
#' \itemize{
#'  \item Model 1: Gompertz regression with Surv(time, status) ~ age
#'  \item Model 2: Gompertz regression with Surv(time, status) ~ age + Biomarkers
#' }
#' @param d4 Dataframe,containing survival information (time, status)
#' @param var Variables for developing biological age (Age and Biomarkers)
#'
#' @return a list
#' \itemize{
#' \item coef1, Model 1 cofficients
#' \item coef2, Model 2 cofficients
#' \item coefficients, Biological age cofficients
#' \item Bioage, Biological age data.frame
#' }

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
#' bioage <- gold_bioage(NHANES4, var)
gold_bioage <- function(d4, var) {
  data <- data.frame(
    time = d4$time,
    status = d4$status,
    d4[, var]
  )
  data <- na.omit(data)

  data1 <- data
  fitg1 <- flexsurv::flexsurvreg(
    formula = survival::Surv(time, status) ~ age, # age
    data = data1, dist = "gompertz"
  )
  coef1 <- fitg1$res[, 1]

  ############
  fitg2 <- flexsurv::flexsurvreg(
    formula = as.formula(paste(
      "survival::Surv(time, status) ~",
      paste(var, collapse = "+")
    )),
    data = data1, dist = "gompertz"
  )
  fitg2

  coef2 <- fitg2$res[, 1]

  ####
  coef21 <- matrix(coef2[3:length(coef2)]) / coef1[3]
  rownames(coef21) <- var
  coef21

  x <- as.matrix(data[, var]) %*% coef21
  pr <- predict(fitg1, type = "hazard", times = c(0))
  pr2 <- predict(fitg2, type = "hazard", times = c(0))
  risk <- pr$.pred_hazard
  risk2 <- pr2$.pred_hazard
  la <- risk / risk2
  x0 <- (log(coef2[2] / coef1[2])) / coef1[3] + mean(log(la) / coef1[3])
  bioage <- x0 + x
  data$bioage <- bioage
  names(coef21)=var
  coef21 <- c(coef21, beta0 = as.numeric(x0))
  list(
    `model1 coef` = coef1,
    `model2 coef` = coef2,
    `coefficients` = coef21,
    `Bioage` = data
  )
}
