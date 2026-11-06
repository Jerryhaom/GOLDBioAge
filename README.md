# GOLDBioAge

An R package for fitting **Gompertz law-based biological age (GOLD BioAge)** and a Cox-based biological age alternative, and for applying frozen GOLD BioAge coefficients to new biomarker data.

[中文使用说明](README.zh-CN.md) · [Simulated Gompertz tutorial](docs/simulated-gompertz.md) · [NHANES tutorial](docs/nhanes.md) · [Published coefficients](Coefficients.md)

## Choose your workflow

| Goal | Required data | Function |
|---|---|---|
| Develop or refit a GOLD BioAge model | Baseline age, biomarkers, survival follow-up, event status | `gold_bioage()` |
| Score new participants with a fitted model | Baseline age and the selected biomarkers, in the training units | `predict_gold_bioage(fit, newdata)` |
| Apply the published NHANES/UKB model | A named published coefficient vector, age, and the corresponding biomarkers | `predict_gold_bioage(coefficients, newdata)` |
| Fit the Cox-based alternative | Baseline age, biomarkers, survival follow-up, event status | `cox_bioage()` |

`gold_bioage()` **refits** the model. It does not automatically use the published NHANES or UKB coefficients. Scoring with frozen coefficients does not require `time` or `status` and does not recenter scores on the new sample.

## Installation

```r
install.packages("remotes")
remotes::install_github("Jerryhaom/GOLDBioAge")
library(GOLDBioAge)
```

Basic fitting uses `flexsurv` and `survival`, installed as package dependencies. Install the optional feature-selection packages only if needed:

```r
install.packages(c("glmnet", "ggplot2"))
```

Version 0.2.0 adds frozen-coefficient prediction, input checks, and runnable tutorials. See [NEWS.md](NEWS.md) for corrections that can affect older feature-selection results.

## Quick start: NHANES

```r
library(GOLDBioAge)
data("NHANES4", package = "GOLDBioAge")

variables <- c("age", "albumin", "alp", "creat", "glucose_mmol",
               "lymph", "mcv", "rdw", "wbc", "ggt")
fit <- gold_bioage(NHANES4, variables)

# Coefficients, including the training intercept beta0.
fit$coefficients

# Complete-case training scores; note the exact field names and capitalization.
head(fit$GOLDBioage[, c("age", "GOLDbioage")])

# Scores in the original input row order; missing predictors produce NA.
NHANES4$GOLDbioage <- predict_gold_bioage(fit, NHANES4)
NHANES4$age_difference <- NHANES4$GOLDbioage - NHANES4$age
```

The bundled data have **54,945 rows and 30 columns**. With the ten predictors above, the complete-case fit uses **39,348 observations and 4,716 deaths**, omitting 15,597 rows. It uses an unweighted analysis of the supplied processed data. See the [NHANES tutorial](docs/nhanes.md) for the data dictionary, missingness checks, fitted coefficients, and a worked example applying the published formula without refitting.

## Example: simulate Gompertz survival data

```r
set.seed(20261106)
n <- 2000L
age <- runif(n, 40, 80)
marker_risk <- rnorm(n)
marker_protective <- rnorm(n)
individual_rate <- exp(-9 + 0.08 * age +
                       0.4 * marker_risk - 0.3 * marker_protective)
event_time <- flexsurv::rgompertz(n, shape = 0.08, rate = individual_rate)
censor_time <- runif(n, 5, 20)
simulated <- data.frame(
  time = pmin(event_time, censor_time),
  status = as.integer(event_time <= censor_time),
  age, marker_risk, marker_protective
)

sim_fit <- gold_bioage(simulated, c("age", "marker_risk", "marker_protective"))
sim_fit$coefficients
```

Here `age` is baseline age in years, while `time` is follow-up in years. The Gompertz hazard increases with follow-up because the shape is positive. The [full simulated-data tutorial](docs/simulated-gompertz.md) explains this parameterization and adds participant-level training/test splitting, frozen prediction, LASSO, and elastic net. These simulated results demonstrate the algorithm and do not establish clinical validity.

## Input requirements and units

For fitting, supply a data frame with numeric columns and a predictor vector containing `age` exactly once. `age` may appear anywhere in the vector.

| Column | Meaning |
|---|---|
| `age` | Baseline chronological age in years |
| `time` | Non-negative follow-up time in one consistent unit |
| `status` | `0` = censored; `1` = observed event/death |
| Biomarkers | Numeric baseline measurements; units and preprocessing must be recorded |

Missing model values are omitted during fitting, with original positions recorded in `fit$rows_used` and `fit$rows_omitted`. Infinite values, factors/characters, duplicate predictors, invalid status coding, and constant predictors are rejected. No imputation, unit conversion, or automatic log transformation is performed. For frozen prediction, missing predictor values return `NA` in their original positions.

For published coefficients, match the units in [Coefficients.md](Coefficients.md). In particular, lymphocytes are entered as a percentage (`30`, not `0.30`), creatinine as mg/dL, glucose as mmol/L, and albumin as g/dL. The corrected NHANES GOLD LYM% coefficient is negative and the GGT coefficient is positive.

## Optional feature selection

```r
lasso_fit <- gold_bioage(
  simulated, c("age", "marker_risk", "marker_protective"),
  feature_selection = TRUE, selection_method = "lasso",
  nfolds = 5, seed = 123
)
lasso_fit$selected_vars
print(lasso_fit$cv_plot)

elasticnet_fit <- gold_bioage(
  simulated, c("age", "marker_risk", "marker_protective"),
  feature_selection = TRUE, selection_method = "elasticnet",
  alpha = 0.5, nfolds = 5, seed = 123
)
```

Selection uses a **biomarker-only Cox model at `lambda.min`** and retains age in the final Gompertz model. It needs at least two candidate biomarkers, strictly positive follow-up times, and at least `nfolds` events. `selection_method = "none"` disables selection. A fixed seed makes cross-validation reproducible while preserving the caller's random-number state. Feature selection and preprocessing belong inside the development sample, before frozen scoring of held-out participants.

## Outputs and interpretation

| Output | Meaning |
|---|---|
| `fit$coefficients` | Named weights in years, including `age` and `beta0` |
| `fit$GOLDBioage$GOLDbioage` | Numeric GOLD BioAge values for complete training rows |
| `fit$selected_vars` | Predictors retained in the fitted model |
| `fit$rows_used`, `fit$rows_omitted` | Original input row positions |
| `fit$models` | Age-only and constrained full Gompertz models for diagnostics |
| `fit$cv_fit`, `fit$cv_plot` | Selection diagnostics, or `NULL` without selection |

The chronological-age coefficient is fixed to 1 after conversion to years. The training intercept centers mean GOLD BioAge on mean chronological age in that complete-case training sample; this is a calibration identity, not evidence of predictive validity. Keep this intercept fixed when scoring new data.

`GOLDbioage - age` is an age difference in years, not an age-adjusted regression residual. A positive value represents a higher score relative to chronological age under the chosen model. Clinical utility and mortality prediction require separate evaluation. Re-fitting or scoring the bundled development data does not provide independent validation, and the package does not implement NHANES survey-design inference.

## Reproducible scripts

From a repository checkout, after installation:

```r
source("examples/01_simulated_gompertz.R")
source("examples/02_nhanes.R")
```

The scripts leave their fitted models and scored data in the R session. Local export commands are included but commented out. The `docs/` and `examples/` folders are available on GitHub; installed function help is available with `?gold_bioage` and `?predict_gold_bioage`.

## Development checks

```r
install.packages(c("roxygen2", "testthat"))
roxygen2::roxygenise()
testthat::test_local()
```

## Reference and contact

Hao M. et al. *Gompertz Law-Based Biological Age (GOLD BioAge): A Simple and Practical Measurement of Biological Ageing to Capture Morbidity and Mortality Risks.* Advanced Science (2025), e01765. [doi:10.1002/advs.202501765](https://doi.org/10.1002/advs.202501765).

[Mathematical derivation](https://jerryhaom.github.io/GOLDBioAge/gompertz_cox_formulas.html) · [Web calculator](https://jerryhaom.github.io/GOLDBioAgeWeb.io/) · Contact: haombio@gmail.com
