# Tutorial: GOLD BioAge with simulated Gompertz data

[Back to README](../README.md) · [中文说明](../README.zh-CN.md) · [Complete executable script](../examples/01_simulated_gompertz.R)

This example develops a model from artificial baseline biomarkers and right-censored survival data, then applies the frozen model to held-out participants. It illustrates the computational workflow, without claiming that simulated biomarker effects represent clinical biology.

## 1. Understand the survival model

Under the [flexsurv Gompertz parameterization](https://chjackson.github.io/flexsurv/reference/Gompertz.html), the follow-up hazard is

```text
h(t | age, X) = rate(age, X) × exp(shape × t)
log(rate)    = -9 + 0.08 × age + 0.4 × marker_risk - 0.3 × marker_protective
shape        = 0.08 per follow-up year
```

`age` is chronological age at baseline. `t` is time elapsed since baseline; it is not attained age. A positive shape gives increasing hazard with follow-up. `marker_risk` and `marker_protective` are arbitrary standardized numerical markers. `noise` has no direct effect in the data-generating model.

The [flexsurv regression documentation](https://chjackson.github.io/flexsurv/reference/flexsurvreg.html) specifies that Gompertz covariates act on log rate. GOLD fits an age-only model to estimate a reference age effect, fixes that **age coefficient** in the final age-plus-biomarker model, divides biomarker coefficients by the reference age effect, and estimates a cohort-centering intercept. The Gompertz follow-up shape is a different parameter and remains estimated.

## 2. Generate participants and censoring

```r
library(GOLDBioAge)
set.seed(20261106)
n <- 2000L
age <- runif(n, 40, 80)
marker_risk <- rnorm(n)
marker_protective <- rnorm(n)
noise <- rnorm(n)
individual_rate <- exp(-9 + 0.08 * age +
                       0.4 * marker_risk - 0.3 * marker_protective)
event_time <- flexsurv::rgompertz(n, shape = 0.08, rate = individual_rate)
censor_time <- runif(n, 5, 20)
simulated <- data.frame(
  time = pmin(event_time, censor_time),
  status = as.integer(event_time <= censor_time),
  age, marker_risk, marker_protective, noise
)

table(simulated$status)
summary(simulated$time)
```

Censoring is generated independently of age and biomarkers in this example. `status = 1` means the death/event time was observed; `0` means the participant was right-censored.

## 3. Split before model development

```r
train_id <- sample.int(n, size = floor(0.8 * n))
train <- simulated[train_id, , drop = FALSE]
test <- simulated[-train_id, , drop = FALSE]
variables <- c("age", "marker_risk", "marker_protective", "noise")

fit <- gold_bioage(train, variables)
fit$coefficients
head(fit$GOLDBioage[, c("age", "GOLDbioage")])
```

This creates 1,600 development participants and 400 held-out participants. There is no preprocessing fitted on the combined sample. If using transformations, imputation, or feature selection in a real study, learn them using the development sample only.

The fitted weights need not exactly recover `0.4/0.08` or `-0.3/0.08`: the reference age coefficient is estimated in a marginal age-only model, then constrained in the full model. Censoring and finite sample variation also affect the estimates. A noise coefficient may be nonzero in a finite sample.

## 4. Score held-out biomarkers without refitting

```r
test$GOLDbioage <- predict_gold_bioage(fit, test[, variables, drop = FALSE])
test$age_difference <- test$GOLDbioage - test$age
summary(test$age_difference)

# Numerical consistency on the training rows.
stopifnot(isTRUE(all.equal(
  predict_gold_bioage(fit, train),
  fit$GOLDBioage$GOLDbioage,
  tolerance = 1e-8
)))
mean(fit$GOLDBioage$GOLDbioage - train$age)
```

The last value is approximately zero because the training intercept centers the scores on training age. There is no requirement for the held-out mean age difference to be zero. Do not estimate a new intercept from the test sample if the goal is to apply the frozen training model.

Prediction uses only predictor columns and frozen coefficients. `time` and `status` are unnecessary for scoring; evaluating mortality prediction would still require independent outcomes and an explicit evaluation design. The example does not estimate clinical performance.

## 5. Optional LASSO and elastic net

Install `glmnet` and `ggplot2` as described in the README before running this block.

```r
lasso_fit <- gold_bioage(
  train, variables, feature_selection = TRUE,
  selection_method = "lasso", nfolds = 5, seed = 123
)
lasso_fit$selected_vars
print(lasso_fit$cv_plot)
test$GOLDbioage_lasso <- predict_gold_bioage(lasso_fit, test)

elasticnet_fit <- gold_bioage(
  train, variables, feature_selection = TRUE,
  selection_method = "elasticnet", alpha = 0.5,
  nfolds = 5, seed = 123
)
elasticnet_fit$selected_vars
```

The screening model is a biomarker-only Cox regression at `lambda.min`. Age is retained in the final Gompertz fit. Screening does not guarantee that the artificial noise variable is dropped. Rerunning with the same seed and software versions gives reproducible folds and selection. The function restores the caller's random-number state.

## 6. Fit the Cox alternative on the same training participants

```r
cox_fit <- cox_bioage(train, variables)
cox_fit$coefficients
cox_fit$intercept
test$Coxbioage <- predict_cox_bioage(cox_fit, test[, variables, drop = FALSE])
test$Cox_age_difference <- test$Coxbioage - test$age
c(mean_train_Cox_age_difference = mean(cox_fit$CoxBioage$Coxbioage - train$age),
  mean_test_Cox_age_difference = mean(test$Cox_age_difference),
  heldout_score_correlation = cor(test$GOLDbioage, test$Coxbioage))
```

Both fits use the same 1,600 development participants and fixed predictor set. The Cox baseline hazard is unspecified, while the Gompertz fit assumes a parametric follow-up hazard. The Cox age effect is fixed through an offset; its biomarker weights are converted to years, and its intercept is centered on training age by default. The held-out mean score difference is free to differ from zero. Correlation between the two held-out scores is descriptive, not a performance-equivalence test or clinical validation.

The final Cox model is distinct from biomarker-only Cox screening used for optional GOLD feature selection. For formulas and the legacy calibration option, see the [Cox tutorial](cox-bioage.md).

## 7. Freeze and reuse the model

```r
# Save the coefficients locally when you want a persistent deployment artifact:
# saveRDS(fit$coefficients, "gold_coefficients.rds")
# frozen <- readRDS("gold_coefficients.rds")
# future_scores <- predict_gold_bioage(frozen, future_biomarker_data)
# saveRDS(cox_fit$score_coefficients, "cox_coefficients.rds")
# future_cox_scores <- predict_cox_bioage(cox_fit$score_coefficients, future_biomarker_data)

sessionInfo()
```

Record the predictor names, units, transformations, event definition, follow-up unit, training population, and package versions alongside the coefficient vector. A reusable score definition requires these details as well as the numeric weights.
