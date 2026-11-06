# GOLDBioAge

Fit **Gompertz law-based biological age (GOLD BioAge)** or a **Cox-derived age score**, and calculate scores for new participants with fixed coefficients.

[中文说明](README.zh-CN.md) · [Equations and code](docs/mathematics-and-code.md) · [Cohort coefficients](Coefficients.md) · [Updates](docs/package-updates.md) · [Online calculator](https://jerryhaom.github.io/GOLDBioAgeWeb.io/)

## Choose your workflow

| What you have | What to do | What you need |
|---|---|---|
| New biomarker measurements | Apply an existing NHANES/UKB GOLD or NHANES/UKB/CHARLS Light model | Age, matching biomarkers, units, and transformations |
| A survival cohort for model development | Fit a new Gompertz GOLD or Cox score | Predictors, follow-up `time`, and event `status` |
| A previously fitted model | Apply its frozen coefficients | Predictor measurements; retain its training intercept |

The calculation is a weighted sum in years:

$$\mathrm{BioAge}=\mathrm{age}+\sum_j w_jX_j+b_0,\qquad
\mathrm{BioAgeDiff}=\mathrm{BioAge}-\mathrm{age}.$$

Age has coefficient 1. $X_j$ is a biomarker in the fitted unit and transformation; $w_j$ is its weight in years per unit; $b_0$ is the intercept. A score does not directly give an absolute mortality probability or remaining lifespan.

## Install

```r
install.packages("remotes")
remotes::install_github("Jerryhaom/GOLDBioAge")
library(GOLDBioAge)
packageVersion("GOLDBioAge")
```

Install `glmnet` and `ggplot2` only if using optional LASSO/elastic-net screening:

```r
install.packages(c("glmnet", "ggplot2"))
```

## Calculate an existing GOLD score

This complete example applies the supplied NHANES coefficients without fitting a model:

```r
library(GOLDBioAge)
beta <- c(age = 1, creat = 5.26909180, glucose_mmol = 0.57969768,
          mcv = 0.33888497, rdw = 2.64452031, albumin = -4.73575224,
          alp = 0.02604091, lymph = -0.20318306, wbc = 0.44588830,
          ggt = 0.06075152, beta0 = -53.62873)
person <- data.frame(age = 50, creat = 1, glucose_mmol = 5,
                     mcv = 90, rdw = 13, albumin = 4,
                     alp = 80, lymph = 30, wbc = 7, ggt = 25)
score <- predict_gold_bioage(beta, person)
c(BioAge = score, BioAgeDiff = score - person$age)
# BioAge = 51.10203967; BioAgeDiff = 1.10203967 years.
```

Use creatinine in mg/dL, glucose in mmol/L, albumin in g/dL, MCV in fL, RDW/lymphocytes in percent, WBC in 1000 cells/uL, and ALP/GGT in U/L. Enter 30 for 30% lymphocytes. The NHANES lymphocyte weight is negative; the GGT weight is positive.

Pass a data frame with one participant per row to score a cohort. Survival outcomes are unnecessary. Prediction preserves row order, returns `NA` for missing predictors, and performs no automatic unit conversion or transformation.

Version 0.2.2 includes all five original coefficient sets:

```r
reference <- read.csv(system.file("extdata", "published_coefficients.csv",
                                  package = "GOLDBioAge"))
selected <- subset(reference, dataset == "ukb" & model == "gold")
beta_ukb <- setNames(selected$coefficient, selected$term)
predict_gold_bioage(beta_ukb, person)
# 50.7472 years.
```

See [Coefficients.md](Coefficients.md) for every coefficient, input unit, Light's `log1p(CRP)` transformation, and executable examples. NHANES Light uses CRP in mg/dL inside the logarithm; UKB/CHARLS use mg/L. Keep units, transformations, weights, and intercept together.

## Fit GOLD and Cox on NHANES

```r
data("NHANES4", package = "GOLDBioAge")
variables <- c("age", "albumin", "alp", "creat", "glucose_mmol",
               "lymph", "mcv", "rdw", "wbc", "ggt")
gold_fit <- gold_bioage(NHANES4, variables)
cox_fit <- cox_bioage(NHANES4, variables)

gold_fit$coefficients                 # Named vector, including beta0.
cox_fit$score_coefficients            # Named vector, including beta0.
NHANES4$GOLD_refit <- predict_gold_bioage(gold_fit, NHANES4)
NHANES4$Cox_refit <- predict_cox_bioage(cox_fit, NHANES4)
```

Both examples use 39,348 complete participants with 4,716 deaths from the processed file of 54,945 rows. GOLD fits Gompertz survival models; Cox uses an unspecified baseline hazard. Both fix the full-model age effect to the age-only estimate, divide biomarker effects by that estimate, and calculate one training intercept:

$$w_j=\theta_j/\beta_A,\qquad
b_0=\operatorname{mean}_{\mathrm{train}}\left(A-\left[A+\sum_jw_jX_j\right]\right).$$

Thus mean training BioAge equals mean training age. Keep this intercept fixed for new participants. Refit coefficients can differ slightly from a supplied rounded model. Refitting the bundled cohort is a worked example; this is an unweighted complete-case fit.

For your own cohort, use numeric predictors, baseline `age` in years, non-negative follow-up `time` in a documented unit, and `status` coded 0 = censored, 1 = event. Missing required fitting values are omitted; `rows_used` and `rows_omitted` identify their original positions. Restrict preprocessing and optional selection to training data.

## Read and run the tutorials

| Guide | Content | Runnable script |
|---|---|---|
| [Mathematics and R code](docs/mathematics-and-code.md) | Symbols, derivation, weighted sum, CRP units | Code beside equations |
| [Simulated Gompertz data](docs/simulated-gompertz.md) | Simulate events/censoring, split training/test, fit GOLD and Cox | [01_simulated_gompertz.R](examples/01_simulated_gompertz.R) |
| [NHANES example](docs/nhanes.md) | Data counts, reproduced weights, fitted and fixed scoring | [02_nhanes.R](examples/02_nhanes.R) |
| [Existing cohort coefficients](Coefficients.md) | NHANES/UKB GOLD and all three Light models | [03_existing_coefficients.R](examples/03_existing_coefficients.R) |
| [Cox fitting](docs/cox-bioage.md) | Offset, calibration, frozen prediction, older scores | Included in both fitting scripts |
| [Formula webpage](https://jerryhaom.github.io/GOLDBioAge/gompertz_cox_formulas.html) | Browser-readable derivation and R example | |

## What changed

Version **0.2.2** adds the installed coefficient reference and reorganized formula/code tutorials. GOLD remains defined by instantaneous hazard matching at baseline, where the follow-up shape factor is 1; no shared-shape constraint is needed. Version **0.2.1** adds frozen Cox prediction and corrects the default Cox intercept convention; `calibration = "legacy"` reproduces older scores. Version **0.2.0** improves validation, frozen GOLD prediction, variable ordering, and optional screening, including retaining the first Cox-selected biomarker. Screening results and Cox score levels can therefore change when updating older analyses.

Read the [update and migration guide](docs/package-updates.md), [NEWS](NEWS.md), and [consistency audit](docs/consistency-audit.md) for verified changes and limits. Standard no-selection GOLD results are retained to numerical precision.

## Citation

Hao M. et al. *Gompertz Law-Based Biological Age (GOLD BioAge): A Simple and Practical Measurement of Biological Ageing to Capture Morbidity and Mortality Risks.* Advanced Science (2025), e01765. [doi:10.1002/advs.202501765](https://doi.org/10.1002/advs.202501765).

Contact: haombio@gmail.com.
