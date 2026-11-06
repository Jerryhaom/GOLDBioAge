# Apply all five existing GOLD/Light cohort models without fitting survival models.
library(GOLDBioAge)
coefficient_file <- system.file("extdata", "published_coefficients.csv", package = "GOLDBioAge")
if (!nzchar(coefficient_file)) stop("Install GOLDBioAge 0.2.2 or later for the coefficient reference.")
coefficient_reference <- read.csv(coefficient_file, stringsAsFactors = FALSE)

# This is a local example helper, not a separate exported package function.
get_beta <- function(dataset, model) {
  selected <- coefficient_reference[
    coefficient_reference$dataset == dataset & coefficient_reference$model == model, ]
  if (!nrow(selected)) stop("No supplied coefficient set for this dataset/model combination.")
  setNames(selected$coefficient, selected$term)
}

# GOLD inputs: age years; albumin g/dL; creat mg/dL; glucose mmol/L;
# mcv fL; lymph and rdw percent; alp and ggt U/L; wbc 1000 cells/uL.
person <- data.frame(age = 50, albumin = 4, alp = 80, creat = 1,
                     glucose_mmol = 5, lymph = 30, mcv = 90, rdw = 13, wbc = 7, ggt = 25)
nhanes_gold <- predict_gold_bioage(get_beta("nhanes", "gold"), person)
ukb_gold <- predict_gold_bioage(get_beta("ukb", "gold"), person)

# LIGHT inputs: same age/creatinine/glucose plus CRP in the cohort-specific unit.
# The same raw concentration is 0.3 mg/dL = 3 mg/L; transform after conversion.
light_nhanes <- data.frame(age = 50, creat = 1, glucose_mmol = 5, log_crp = log1p(0.3))
light_mg_l <- data.frame(age = 50, creat = 1, glucose_mmol = 5, log_crp = log1p(3))
nhanes_light <- predict_gold_bioage(get_beta("nhanes", "light"), light_nhanes)
ukb_light <- predict_gold_bioage(get_beta("ukb", "light"), light_mg_l)
charls_light <- predict_gold_bioage(get_beta("charls", "light"), light_mg_l)

example_scores <- c(nhanes_gold = nhanes_gold, ukb_gold = ukb_gold,
                    nhanes_light = nhanes_light, ukb_light = ukb_light,
                    charls_light = charls_light)
print(round(example_scores, 8))
stopifnot(isTRUE(all.equal(unname(example_scores),
  c(51.10203967, 50.7472, 50.439978417530952, 53.382622562072712, 52.173009052401632),
  tolerance = 1e-10)))

# Example of manual linear scoring with the same named coefficients.
beta <- get_beta("nhanes", "gold")
predictors <- setdiff(names(beta), "beta0")
manual_score <- as.vector(as.matrix(person[, predictors, drop = FALSE]) %*%
                           beta[predictors] + beta[["beta0"]])
stopifnot(isTRUE(all.equal(manual_score, nhanes_gold)))
