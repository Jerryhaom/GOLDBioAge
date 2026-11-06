# GOLDBioAge 0.2.2

- Reorganize the English and Chinese READMEs around two tasks: apply existing cohort coefficients, or fit a new model. Pair the equations with runnable R code and worked numerical results.
- Add an installed reference file, `extdata/published_coefficients.csv`, containing the five existing GOLD/Light coefficient sets (37 terms) with input units, transformations, and supplied precision. Keep the original cohort values; retain the supplied higher precision for NHANES GOLD.
- Expand `Coefficients.md` to include the original NHANES/UKB GOLD and NHANES/UKB/CHARLS Light tables, fixed-score examples, CRP unit handling, and a separately labeled NHANES Cox refit example.
- Add a formula-to-code guide, a reader-oriented update/migration guide, and an executable example for all five existing coefficient sets. Fitting and scoring algorithms remain as in 0.2.1. Clarify that GOLD matches instantaneous hazards at baseline, where both follow-up shape factors equal 1 without a common-shape constraint.

# GOLDBioAge 0.2.1

- Add Cox fitting, frozen scoring, and Gompertz/Cox comparisons to the simulated and NHANES examples, both READMEs, and a dedicated Cox tutorial. Update the online methodological derivation with the package's explicit calibration convention and executable example.
- Add `predict_cox_bioage()`, `score_coefficients`, and the fitted Cox models. Preserve the existing coefficient-matrix and result-data fields.
- **Cox score calibration changes by default:** compute the intercept as `mean(age - linear_score)` on the training sample, avoiding inconsistent centering of Cox linear predictions and age offsets. This changes score levels, not the fitted age/biomarker coefficients. Use `cox_bioage(..., calibration = "legacy")` to reproduce the intercept and scores used through 0.2.0. Saved legacy objects retain their original intercept during prediction.
- Clarify that the Cox-derived score is based on relative risk, keeps its training intercept during transfer, and does not estimate an absolute mortality probability. Distinguish final Cox fitting from Cox-based screening inside the Gompertz workflow.

# GOLDBioAge 0.2.0

- Add executable simulated Gompertz and bundled NHANES examples, detailed Markdown tutorials, and a Chinese getting-started guide.
- Add `predict_gold_bioage()` to apply frozen fitted or published coefficients to new biomarker data without survival outcomes; preserve input row order and return `NA` for missing predictors.
- Accept `age` anywhere in the predictor list and match the constrained age parameter by name. The constraint is on the chronological-age coefficient, not the Gompertz shape.
- Correct Cox LASSO/elastic-net selection: Cox `glmnet` has no intercept row, so the first biomarker coefficient must not be discarded. **Selected variables and fitted scores may change for users of feature selection.** Screening remains biomarker-only, uses `lambda.min`, and always retains age in the final fit.
- Honor `selection_method = "none"`, declare optional dependencies, require them only for selection, and add a cross-validation seed that restores the caller's RNG state.
- Retain Breslow Cox tie handling when using newer `glmnet` versions whose defaults may change.
- Validate survival coding and model columns, report omitted row positions, and keep legacy result field names. GOLD scores are returned as a plain numeric column.
- Compute the existing GOLD cohort-centering intercept on the linear scale. Standard no-selection fits agree with the previous implementation; scoring new data retains the training intercept.
- Improve Cox variable-order handling and input checks while retaining its coefficient-matrix return format.
- Refresh installed help pages and add automated regression tests and package checks.
