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
