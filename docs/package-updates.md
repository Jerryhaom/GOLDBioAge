# Package updates and migration

[Getting started](../README.md) · [中文说明](../README.zh-CN.md) · [Full change log](../NEWS.md)

## Version 0.2.2: make the model definition usable

This release packages the five original GOLD/Light coefficient sets as `inst/extdata/published_coefficients.csv`: NHANES and UKB GOLD, plus NHANES, UKB, and CHARLS Light. Each row records the term, coefficient, raw input unit, transformation, and supplied decimal precision. The fitting and calibration algorithms are unchanged from 0.2.1. GOLD is evaluated by instantaneous hazard matching at baseline: both follow-up shape factors are then 1, even when their fitted shapes differ. A shared-shape constraint is not required. The [equations and code guide](mathematics-and-code.md) includes optional direct parameter-fixing code as a separate modeling example.

The documentation follows the order a user needs: choose an existing model or refit; read its equation; prepare the required units; run a complete calculation; then transfer frozen coefficients to new data. Start with:

| Task | Guide | Executable example |
|---|---|---|
| Apply original cohort coefficients | [Coefficient reference](../Coefficients.md) | [03_existing_coefficients.R](../examples/03_existing_coefficients.R) |
| Understand equations and calculation | [Mathematics and R code](mathematics-and-code.md) | Code alongside each equation |
| Develop a model with simulated survival | [Gompertz simulation](simulated-gompertz.md) | [01_simulated_gompertz.R](../examples/01_simulated_gompertz.R) |
| Refit bundled NHANES | [NHANES tutorial](nhanes.md) | [02_nhanes.R](../examples/02_nhanes.R) |
| Fit and transfer a Cox score | [Cox tutorial](cox-bioage.md) | Included in both fitting scripts |
| Check implementation agreement | [Consistency audit](consistency-audit.md) | R and calculator regression checks |

The fixed NHANES GOLD model uses negative LYM% and positive GGT coefficients. Its supplied biomarker weights retain eight decimals and its intercept five decimals. The UKB and Light sets retain their supplied four-decimal values. Tables display fewer digits where stated; refitting can produce small differences from a rounded fixed reference.

## Version 0.2.1: Cox intercept calibration

The default Cox score now calculates the intercept directly as

$$b_0=\operatorname{mean}_{\mathrm{training}}(A-s),\qquad
s=A+\sum_jw_jX_j.$$

Earlier code combined centered `predict.coxph(type = "lp")` predictions with uncentered age and biomarkers. The fitted Cox weights did not change, but the default intercept and score level can change. On the bundled GOLD-complete NHANES sample, the current intercept is about -53.56354467 and mean score about 49.49797 years. The earlier intercept is about -49.49796686 and mean score about 53.56354 years.

```r
library(GOLDBioAge)
data("NHANES4", package = "GOLDBioAge")
variables <- c("age", "albumin", "alp", "creat", "glucose_mmol",
               "lymph", "mcv", "rdw", "wbc", "ggt")

current <- cox_bioage(NHANES4, variables)
legacy <- cox_bioage(NHANES4, variables, calibration = "legacy")
stopifnot(isTRUE(all.equal(current$coefficients, legacy$coefficients)))
c(current = current$intercept, legacy = legacy$intercept)

# Keep the intercept when saving a portable model.
frozen <- current$score_coefficients
scores <- predict_cox_bioage(frozen, NHANES4)
```

`predict_cox_bioage()` also accepts a saved older fit and uses its stored intercept. This retains the original scores without silently recalibrating them. Existing matrix `fit$coefficients` and scalar `fit$intercept` remain available; `fit$score_coefficients` combines them.

## Version 0.2.0: fitting, selection, and prediction

The implementation validates predictors and survival fields, handles the age constraint by name, returns fitted-row positions, and provides `predict_gold_bioage()` for row-aligned frozen scoring. GOLD's intercept is calculated on the linear scale, avoiding unnecessary exponentiation. The no-selection NHANES GOLD refit retains the prior results to numerical precision.

Optional LASSO/elastic-net screening now retains every named nonzero Cox coefficient. The previous code could drop the first biomarker because Cox `glmnet` output has no intercept row. Screening uses a biomarker-only Cox model at `lambda.min`, followed by a final Gompertz fit with age retained. `selection_method = "none"` disables screening. Cross-validation has a reproducible seed and preserves the caller's random-number state. `glmnet` and `ggplot2` are needed only for this optional step.

When reporting analyses, record the package version, coefficient source, intercept convention, units, transformations, missing-data handling, and whether a model was refitted or transferred. An update that changes selection or calibration can change scores even if function names remain the same.
