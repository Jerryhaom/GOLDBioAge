# Cox biological age: fitting, calibration, and frozen scoring

[README](../README.md) · [中文说明](../README.zh-CN.md) · [Online Gompertz/Cox derivation](https://jerryhaom.github.io/GOLDBioAge/gompertz_cox_formulas.html)

This tutorial uses GOLDBioAge **0.2.1 or later**. It covers the final Cox-derived biological age model. Cox-LASSO screening inside `gold_bioage(feature_selection = TRUE)` is a separate operation whose final fitted model is Gompertz.

## 1. The two Cox models

First estimate the reference age effect from an age-only model. Then fix that effect through an offset in the biomarker model:

$$h_1(t\mid A)=h_{01}(t)\exp(\beta_A A),$$

$$h_2(t\mid A,X)=h_{02}(t)\exp\left(\beta_A A+\sum_j\theta_j X_j\right).$$

The Cox baseline hazards are unspecified and need not be equal. The package uses Efron handling of tied event times. Unlike the Gompertz approach, it does not estimate a parametric follow-up shape.

Convert the fitted effects to years and define a score on a common training reference:

$$s_i=A_i+\sum_j\frac{\theta_j}{\beta_A}X_{ij},\qquad
b_0=\frac{1}{n_{\mathrm{train}}}\sum_{i\in\mathrm{train}}(A_i-s_i),$$

$$\mathrm{CoxBioAge}_i=s_i+b_0.$$

Thus the training mean score equals training mean age. This identity comes from calibration; it is not a validation result. The reference age effect must be nonzero. The intercept remains fixed when scoring new participants.

Equivalently, calculate `b0 = mean(lp1_raw - lp2_raw) / beta_age` using **uncentered** linear predictors. [The survival documentation](https://stat.ethz.ch/R-manual/R-devel/library/survival/html/predict.coxph.html) explains that default `predict.coxph(type = "lp")` values are centered. Combining those values with an uncentered age offset and raw covariate weights can change the score origin. The default implementation avoids this ambiguity by directly computing `mean(age - linear_score)`.

This is an age-scaled relative-risk score, not an absolute death probability. Equating absolute hazards at a chosen time would additionally involve the ratio of the two baseline hazards. The package does not assume that ratio is one or that it is constant across follow-up.

## 2. Simulated Gompertz data: fit Cox and compare on held-out participants

```r
library(GOLDBioAge)
set.seed(20261106)
n <- 2000L
age <- runif(n, 40, 80)
marker_risk <- rnorm(n)
marker_protective <- rnorm(n)
noise <- rnorm(n)
rate <- exp(-9 + 0.08 * age + 0.4 * marker_risk - 0.3 * marker_protective)
event_time <- flexsurv::rgompertz(n, shape = 0.08, rate = rate)
censor_time <- runif(n, 5, 20)
simulated <- data.frame(time = pmin(event_time, censor_time),
                        status = as.integer(event_time <= censor_time),
                        age, marker_risk, marker_protective, noise)
train_id <- sample.int(n, floor(0.8 * n))
train <- simulated[train_id, ]
test <- simulated[-train_id, ]
variables <- c("age", "marker_risk", "marker_protective", "noise")

cox_fit <- cox_bioage(train, variables)
gold_fit <- gold_bioage(train, variables)
test$Coxbioage <- predict_cox_bioage(cox_fit, test[, variables])
test$GOLDbioage <- predict_gold_bioage(gold_fit, test[, variables])
cox_fit$score_coefficients
c(n_train = nrow(train), n_test = nrow(test),
  mean_train_difference = mean(cox_fit$CoxBioage$Coxbioage - train$age),
  mean_test_difference = mean(test$Coxbioage - test$age),
  heldout_score_correlation = cor(test$Coxbioage, test$GOLDbioage))
```

The training mean difference is approximately zero. The test mean need not be zero. Both methods use exactly the same training participants and candidate predictors; their score correlation describes agreement, not equivalent predictive accuracy. The artificial data do not establish clinical validity. See the [full simulation tutorial](simulated-gompertz.md) and [executable script](../examples/01_simulated_gompertz.R).

## 3. NHANES: fit on the same complete cases

```r
data("NHANES4", package = "GOLDBioAge")
variables <- c("age", "albumin", "alp", "creat", "glucose_mmol",
               "lymph", "mcv", "rdw", "wbc", "ggt")
cox_fit <- cox_bioage(NHANES4, variables)
gold_fit <- gold_bioage(NHANES4, variables)
stopifnot(identical(cox_fit$rows_used, gold_fit$rows_used))
w <- cox_fit$score_coefficients
round(cbind(Gompertz = gold_fit$coefficients[names(w)], Cox = w), 8)
c(n = nrow(cox_fit$CoxBioage), events = sum(cox_fit$CoxBioage$status),
  mean_age = mean(cox_fit$CoxBioage$age),
  mean_cox_age = mean(cox_fit$CoxBioage$Coxbioage))
```

Using the bundled data, the common complete-case sample contains **39,348 participants and 4,716 deaths**. The verified Cox weights, displayed to eight decimal places, are:

| Term | Cox weight |
|---|---:|
| `age` | 1.00000000 |
| `albumin` | -4.86554313 |
| `alp` | 0.02708691 |
| `creat` | 5.27188788 |
| `glucose_mmol` | 0.56861920 |
| `lymph` | -0.20286375 |
| `mcv` | 0.34176961 |
| `rdw` | 2.65618489 |
| `wbc` | 0.45024815 |
| `ggt` | 0.06130221 |
| `beta0` | -53.56354467 |

Mean chronological age and mean calibrated Cox BioAge are both approximately **49.49797 years**. Reference environment: R 4.5.3, GOLDBioAge 0.2.1, survival 3.8-6; small final-digit differences may occur across versions. Biomarker units and missingness are described in the [NHANES tutorial](nhanes.md). This is an unweighted complete-case analysis of the distributed file, without NHANES survey-design inference.

## 4. Apply the fitted Cox model without survival outcomes

```r
biomarkers <- NHANES4[, variables, drop = FALSE]
scores <- predict_cox_bioage(cox_fit, biomarkers)
stopifnot(length(scores) == nrow(NHANES4))
stopifnot(isTRUE(all.equal(scores[cox_fit$rows_used],
                           cox_fit$CoxBioage$Coxbioage, tolerance = 1e-8)))

# A portable named coefficient vector includes the training intercept beta0.
frozen <- cox_fit$score_coefficients
stopifnot(isTRUE(all.equal(predict_cox_bioage(frozen, biomarkers), scores)))
# saveRDS(frozen, "NHANES_Cox_coefficients.rds")
```

Prediction needs the numeric age and biomarker columns, with matching units and transformations. It retains original row order, returns `NA` for missing predictors, and keeps the training intercept. It does not require `time` or `status`.

`cox_fit$coefficients` retains the historical one-column matrix and does **not** include an intercept. Use `cox_fit$score_coefficients` to save a complete frozen score, or pass the whole fitted object to `predict_cox_bioage()`.

## 5. Reproduce results from version 0.2.0 or earlier

```r
legacy_fit <- cox_bioage(NHANES4, variables, calibration = "legacy")
stopifnot(isTRUE(all.equal(legacy_fit$coefficients, cox_fit$coefficients)))
c(centered_intercept = cox_fit$intercept,
  legacy_intercept = legacy_fit$intercept,
  centered_mean_score = mean(cox_fit$CoxBioage$Coxbioage),
  legacy_mean_score = mean(legacy_fit$CoxBioage$Coxbioage))
```

The legacy intercept is approximately **-49.49796686**, with mean score approximately **53.56354 years** on this NHANES sample. The corrected default intercept is approximately **-53.56354467**, with mean score **49.49797 years**. The fitted weights are identical; the difference is a constant score shift. Rank correlations and ordering are unchanged by that constant, but absolute age differences change.

Use the legacy option when reproducing old analyses and record the calibration convention with exported coefficients. `predict_cox_bioage()` accepts saved 0.2.0 objects that contain the historical matrix and intercept, and retains the saved intercept rather than changing their scores. The final Cox fit does not take the GOLD function's feature-selection arguments.
