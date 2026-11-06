#' NHANES example data for biological age modeling
#'
#' A bundled, processed data frame with 54,945 observations and 30 columns.
#' This is example analysis data, not a complete raw NHANES download or a
#' survey-design object. Re-fitting it illustrates the algorithm; scoring it
#' with the published NHANES coefficients is not independent validation.
#'
#' @format A data frame with 54,945 rows and 30 columns. Key columns are:
#' \describe{
#'   \item{time}{Non-negative follow-up time. Use one consistent follow-up unit
#'     within a fit; this field is not chronological age.}
#'   \item{status}{Event indicator: 0 = censored, 1 = death.}
#'   \item{sampleID}{Participant label; never include this as a predictor.}
#'   \item{age}{Baseline chronological age, years.}
#'   \item{albumin}{Albumin, g/dL.}
#'   \item{alp}{Alkaline phosphatase, U/L.}
#'   \item{creat}{Creatinine, mg/dL.}
#'   \item{glucose_mmol}{Glucose, mmol/L.}
#'   \item{lymph}{Lymphocyte percentage (30 means 30 percent).}
#'   \item{mcv}{Mean cell volume, fL.}
#'   \item{rdw}{Red cell distribution width, percent.}
#'   \item{wbc}{White blood cell count, 1000 cells/uL.}
#'   \item{ggt}{Gamma-glutamyl transferase, U/L.}
#'   \item{lncrp}{Pre-transformed CRP field. Do not assume its units or transform
#'     match a published Light BioAge formula without checking the data provenance.}
#'   \item{bmi}{Additional baseline body mass index field.}
#'   \item{bun}{Additional blood urea nitrogen field.}
#'   \item{ttbl}{Additional total bilirubin field.}
#'   \item{uap}{Additional distributed biomarker field with incomplete upstream metadata.}
#'   \item{basopa}{Additional basophil field.}
#'   \item{eosnpa}{Additional eosinophil field.}
#'   \item{monopa}{Additional monocyte field.}
#'   \item{neut}{Additional neutrophil field.}
#'   \item{rbc}{Additional red blood cell field.}
#'   \item{dbp}{Additional diastolic blood pressure field.}
#'   \item{sbp}{Additional systolic blood pressure field.}
#'   \item{meanbp}{Additional mean blood pressure field.}
#'   \item{hba1c}{Additional glycated hemoglobin field.}
#'   \item{hdl}{Additional HDL cholesterol field.}
#'   \item{ldl}{Additional LDL cholesterol field.}
#'   \item{totchol}{Additional total cholesterol field.}
#' }
#' The additional fields are not used in the GOLD tutorial. Their units and
#' processing should be verified from source records before fitting new models.
#' Missing biomarker values are present. No imputation is performed by the package.
#' @source Bundled NHANES-derived data distributed with GOLDBioAge.
#' @examples
#' data(NHANES4)
#' dim(NHANES4)
#' table(NHANES4$status)
"NHANES4"
