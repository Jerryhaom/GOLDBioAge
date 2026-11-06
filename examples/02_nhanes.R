# Run after installing GOLDBioAge: source("examples/02_nhanes.R")
library(GOLDBioAge)
data("NHANES4", package = "GOLDBioAge")

# 1. Inspect sample size, event coding, and missingness before fitting.
variables <- c("age", "albumin", "alp", "creat", "glucose_mmol",
               "lymph", "mcv", "rdw", "wbc", "ggt")
required <- c("time", "status", variables)
print(dim(NHANES4))
print(table(NHANES4$status, useNA = "ifany"))
print(colSums(is.na(NHANES4[, required])))
complete_rows <- complete.cases(NHANES4[, required])
print(c(total = nrow(NHANES4), complete = sum(complete_rows),
        excluded = sum(!complete_rows), events_complete = sum(NHANES4$status[complete_rows])))

# 2. Refit a model using the bundled cohort (unweighted complete-case analysis).
nhanes_fit <- gold_bioage(NHANES4, variables)
print(round(nhanes_fit$coefficients, 8))
print(summary(nhanes_fit$GOLDBioage$GOLDbioage))

# 3. Restore fitted scores to original row positions without a join on row names.
nhanes_scored <- NHANES4
nhanes_scored$GOLDbioage_refit <- NA_real_
nhanes_scored$GOLDbioage_refit[nhanes_fit$rows_used] <- nhanes_fit$GOLDBioage$GOLDbioage
nhanes_scored$GOLDbioage_frozen <- predict_gold_bioage(nhanes_fit, nhanes_scored)
stopifnot(isTRUE(all.equal(nhanes_scored$GOLDbioage_refit,
                           nhanes_scored$GOLDbioage_frozen, tolerance = 1e-8)))

# 4. Alternatively apply the published NHANES GOLD coefficients without refitting.
# Units: age years; albumin g/dL; alp and ggt U/L; creat mg/dL;
# glucose_mmol mmol/L; lymph and rdw percent; mcv fL; wbc 1000 cells/uL.
published_nhanes <- c(
  age = 1, albumin = -4.73575224, alp = 0.02604091, creat = 5.26909180,
  glucose_mmol = 0.57969768, lymph = -0.20318306, mcv = 0.33888497,
  rdw = 2.64452031, wbc = 0.44588830, ggt = 0.06075152, beta0 = -53.62873
)
nhanes_scored$GOLDbioage_published <- predict_gold_bioage(
  published_nhanes, nhanes_scored[, variables, drop = FALSE]
)
nhanes_scored$age_difference <- nhanes_scored$GOLDbioage_published - nhanes_scored$age
print(summary(nhanes_scored$age_difference))

# A worked individual example with exactly the same biomarker units.
person <- data.frame(age = 50, albumin = 4, alp = 80, creat = 1,
                     glucose_mmol = 5, lymph = 30, mcv = 90, rdw = 13, wbc = 7, ggt = 25)
person_score <- predict_gold_bioage(published_nhanes, person)
print(person_score) # 51.10203967 years (51.1 at one decimal place)
stopifnot(abs(person_score - 51.10203967) < 1e-8)

# 5. Fit Cox on the same variables and complete cases, then compare weights/scores.
nhanes_cox_fit <- cox_bioage(NHANES4, variables)
stopifnot(identical(nhanes_cox_fit$rows_used, nhanes_fit$rows_used))
nhanes_scored$Coxbioage <- predict_cox_bioage(nhanes_cox_fit, nhanes_scored)
cox_weights <- nhanes_cox_fit$score_coefficients
print(round(cbind(Gompertz = nhanes_fit$coefficients[names(cox_weights)],
                  Cox = cox_weights), 8))
print(c(mean_age = mean(nhanes_cox_fit$CoxBioage$age),
        mean_Cox_age = mean(nhanes_cox_fit$CoxBioage$Coxbioage),
        score_correlation = cor(nhanes_scored$GOLDbioage_refit,
                                nhanes_scored$Coxbioage, use = "complete.obs")))
stopifnot(isTRUE(all.equal(nhanes_scored$Coxbioage[nhanes_cox_fit$rows_used],
                           nhanes_cox_fit$CoxBioage$Coxbioage, tolerance = 1e-8)))

# Optional local export (uncomment to write files in your current directory):
# write.csv(nhanes_scored, "NHANES_GOLDBioAge.csv", row.names = FALSE)
# saveRDS(nhanes_fit$coefficients, "NHANES_refitted_coefficients.rds")
# saveRDS(nhanes_cox_fit$score_coefficients, "NHANES_Cox_coefficients.rds")
# Scoring the development dataset is not independent validation.
