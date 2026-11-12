# Reproduce NHANES GOLD and Light BioAge with GOLDBioAge >= 0.2.2.
# See Usage.md for expected coefficients, units, and the original-model comparison.

# 1. Load the complete bundled dataset.
library(GOLDBioAge)
data("NHANES4", package = "GOLDBioAge")
stopifnot(utils::packageVersion("GOLDBioAge") >= "0.2.2")
stopifnot(nrow(NHANES4) == 54945L)

# 2. Refit GOLD using age and all nine biomarkers, without feature selection.
gold_variables <- c("age", "albumin", "alp", "creat", "glucose_mmol",
                    "lymph", "mcv", "rdw", "wbc", "ggt")
gold_fit <- gold_bioage(NHANES4, gold_variables, feature_selection = FALSE)
print(round(gold_fit$coefficients, 8))
gold_counts <- c(n = length(gold_fit$rows_used),
                 deaths = sum(NHANES4$status[gold_fit$rows_used]))
print(gold_counts)
stopifnot(gold_counts[["n"]] == 39348L, gold_counts[["deaths"]] == 4716L)

# 3. Refit Light from the full dataset, using its own complete cases.
# lncrp is already ln(CRP + 1), with raw CRP in mg/dL. Do not log it again.
light_variables <- c("age", "creat", "glucose_mmol", "lncrp")
light_fit <- gold_bioage(NHANES4, light_variables, feature_selection = FALSE)
print(round(light_fit$coefficients, 8))
light_counts <- c(n = length(light_fit$rows_used),
                  deaths = sum(NHANES4$status[light_fit$rows_used]))
print(light_counts)
stopifnot(light_counts[["n"]] == 38001L, light_counts[["deaths"]] == 6373L)

# 4. Compare newly fitted coefficients with the original fixed reference.
# This file is read only after fitting and does not supply the fitted coefficients.
reference <- read.csv(system.file("extdata", "published_coefficients.csv",
                                  package = "GOLDBioAge"))
compare_reference <- function(fit, model_name) {
  original <- reference[reference$dataset == "nhanes" &
                          reference$model == model_name, ]
  # The reference calls the transformed CRP term log_crp; NHANES4 calls it lncrp.
  fitted_terms <- ifelse(original$term == "log_crp", "lncrp", original$term)
  refitted <- unname(fit$coefficients[fitted_terms])
  data.frame(term = original$term, original = original$coefficient,
             refitted = refitted, difference = refitted - original$coefficient)
}
gold_comparison <- compare_reference(gold_fit, "gold")
light_comparison <- compare_reference(light_fit, "light")
print(gold_comparison, digits = 10, row.names = FALSE)
print(light_comparison, digits = 10, row.names = FALSE)

# Numerical checks allow the supplied intercept precision and minor fit differences.
stopifnot(all(abs(gold_comparison$difference) < 1e-5))
stopifnot(all(abs(light_comparison$difference) < 1e-4))

# 5. Calculate scores from the newly fitted models.
# Units: age years; albumin g/dL; creat mg/dL; glucose_mmol mmol/L;
# alp and ggt U/L; mcv fL; lymph and rdw percent; wbc 1000 cells/uL.
person <- data.frame(age = 50, albumin = 4, alp = 80, creat = 1,
                     glucose_mmol = 5, lymph = 30, mcv = 90, rdw = 13, wbc = 7, ggt = 25,
                     lncrp = log1p(0.3))  # A new raw CRP measurement: 0.3 mg/dL.
scores <- c(GOLD = predict_gold_bioage(gold_fit, person),
            Light = predict_gold_bioage(light_fit, person))
print(round(scores, 4))

# Prediction retains every original row; missing predictors yield NA scores.
nhanes_scores <- data.frame(
  age = NHANES4$age,
  GOLD = predict_gold_bioage(gold_fit, NHANES4),
  Light = predict_gold_bioage(light_fit, NHANES4)
)
stopifnot(nrow(nhanes_scores) == nrow(NHANES4))

# 6. Record the software environment used for your reproduction.
print(sessionInfo())
