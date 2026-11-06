# Tutorial: GOLD BioAge with the bundled NHANES data

[Back to README](../README.md) · [中文说明](../README.zh-CN.md) · [Complete executable script](../examples/02_nhanes.R)

There are two distinct tasks: **refit a model using the bundled NHANES cohort**, or **apply the published NHANES coefficients**. The first requires survival follow-up and estimates coefficients. The second needs only baseline age and biomarkers with matching units.

## 1. Inspect the bundled data

```r
library(GOLDBioAge)
data("NHANES4", package = "GOLDBioAge")
dim(NHANES4)
table(NHANES4$status, useNA = "ifany")

variables <- c("age", "albumin", "alp", "creat", "glucose_mmol",
               "lymph", "mcv", "rdw", "wbc", "ggt")
required <- c("time", "status", variables)
colSums(is.na(NHANES4[, required]))
complete_rows <- complete.cases(NHANES4[, required])
c(total = nrow(NHANES4), complete = sum(complete_rows),
  excluded = sum(!complete_rows), events_complete = sum(NHANES4$status[complete_rows]))
```

For the bundled file, these counts are:

| Check | Value |
|---|---:|
| All observations | 54,945 |
| Columns | 30 |
| All observed deaths (`status = 1`) | 9,193 |
| Complete cases for GOLD fitting | 39,348 |
| Excluded for missing model values | 15,597 |
| Deaths among complete cases | 4,716 |

These counts describe the distributed processed file and this predictor set. They do not imply that every NHANES participant, survey wave, or paper analysis used the same sample.

## 2. Match columns and units

| Column | Definition | Unit for published NHANES GOLD |
|---|---|---|
| `age` | Baseline chronological age | years |
| `albumin` | Albumin | g/dL |
| `alp` | Alkaline phosphatase | U/L |
| `creat` | Creatinine | mg/dL |
| `glucose_mmol` | Glucose | mmol/L |
| `lymph` | Lymphocytes | percent, e.g. 30 for 30% |
| `mcv` | Mean cell volume | fL |
| `rdw` | Red cell distribution width | percent |
| `wbc` | White blood cell count | 1000 cells/uL |
| `ggt` | Gamma-glutamyl transferase | U/L |
| `time` | Follow-up elapsed from baseline | required only for fitting |
| `status` | Death/event indicator | 0 = censored, 1 = death |

The repository does not provide a complete raw-to-analysis preprocessing dictionary or explicit metadata for the stored follow-up unit. This tutorial retains `time` as distributed; custom cohorts must document and use a consistent follow-up unit, while `age` remains in years. Do not infer a unit from the column name alone.

`sampleID` is an identifier and is excluded from the predictors. `lncrp` is already transformed in the bundled file; do not use it as raw CRP or assume its transform/units match the published Light BioAge model. The example below uses the GOLD biomarker set and requires no CRP transformation.

## 3. Refit GOLD BioAge

```r
fit <- gold_bioage(NHANES4, variables)
round(fit$coefficients, 8)
summary(fit$GOLDBioage$GOLDbioage)
length(fit$rows_used)
length(fit$rows_omitted)
```

A verified full-data fit gives the following coefficients (displayed to eight decimal places):

| Term | Fitted coefficient |
|---|---:|
| `age` | 1.00000000 |
| `albumin` | -4.73575224 |
| `alp` | 0.02604091 |
| `creat` | 5.26909180 |
| `glucose_mmol` | 0.57969768 |
| `lymph` | -0.20318306 |
| `mcv` | 0.33888497 |
| `rdw` | 2.64452031 |
| `wbc` | 0.44588830 |
| `ggt` | 0.06075152 |
| `beta0` | -53.62873289 |

Small last-digit differences can occur across optimizer and dependency versions. The age coefficient is constrained to 1 after conversion to years. Missing rows are removed only for required model columns. The package does not impute missingness, fit survey weights, or incorporate strata and clusters; this example is an unweighted complete-case fit.

The displayed reference was verified with R 4.5.3, GOLDBioAge 0.2.0, flexsurv 2.3.2, and survival 3.8-6. The no-selection refit agrees with the previous implementation to numerical precision.

## 4. Preserve alignment with the original data

```r
scored <- NHANES4
scored$GOLDbioage_refit <- NA_real_
scored$GOLDbioage_refit[fit$rows_used] <- fit$GOLDBioage$GOLDbioage

# Apply the fitted coefficients in the original row order.
scored$GOLDbioage_frozen <- predict_gold_bioage(fit, scored)
stopifnot(isTRUE(all.equal(
  scored$GOLDbioage_refit, scored$GOLDbioage_frozen, tolerance = 1e-8
)))
```

Here the two columns agree because the bundled survival fields are complete, so missing biomarker rows coincide with omitted fitting rows. In a custom dataset, a participant missing only `time` or `status` can still be scored if all selected predictors are available. Use `rows_used` to align fitted rows; prediction always preserves the new data's row order.

## 5. Apply the published NHANES formula without refitting

```r
published_nhanes <- c(
  age = 1, albumin = -4.73575224, alp = 0.02604091, creat = 5.26909180,
  glucose_mmol = 0.57969768, lymph = -0.20318306, mcv = 0.33888497,
  rdw = 2.64452031, wbc = 0.44588830, ggt = 0.06075152, beta0 = -53.62873
)
scored$GOLDbioage_published <- predict_gold_bioage(
  published_nhanes, scored[, variables, drop = FALSE]
)
scored$age_difference <- scored$GOLDbioage_published - scored$age
summary(scored$age_difference)
```

The supplied published intercept is rounded to `-53.62873`; the refitted intercept above has more digits. This small precision difference is expected. [Coefficients.md](../Coefficients.md) displays biomarker coefficients rounded to four decimal places. The NHANES GOLD lymphocyte coefficient is **negative** and the GGT coefficient is **positive**.

For a single participant:

```r
person <- data.frame(age = 50, albumin = 4, alp = 80, creat = 1,
                     glucose_mmol = 5, lymph = 30, mcv = 90, rdw = 13, wbc = 7, ggt = 25)
person_score <- predict_gold_bioage(published_nhanes, person)
person_score
stopifnot(abs(person_score - 51.10203967) < 1e-8)
```

The result is **51.10203967 years**, displayed as **51.1 years** by the web calculator. All units are those listed above. No follow-up or event status is needed to apply this fixed score. Rename external columns and convert their units before prediction; the function does neither automatically.

## 6. Interpret and export

```r
# Optional local exports:
# write.csv(scored, "NHANES_GOLDBioAge.csv", row.names = FALSE)
# saveRDS(fit$coefficients, "NHANES_refitted_coefficients.rds")
sessionInfo()
```

`GOLDbioage - age` is a score difference in years, not a regression-residual age-acceleration measure. The fitted training mean difference is zero by construction. Keep frozen coefficients, including their intercept, unchanged when transferring the model to new data. Re-fitting the bundled cohort or applying its published coefficients back to it is a worked example, not independent validation or proof of survey-representative clinical performance.
