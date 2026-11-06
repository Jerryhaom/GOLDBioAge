# Existing cohort coefficients: GOLD, Light, and the NHANES Cox example

[Getting started](README.md) · [中文说明](README.zh-CN.md) · [Equations and R code](docs/mathematics-and-code.md) · [Package updates](docs/package-updates.md)

Use this page when you want to calculate a score with an existing model. The coefficients in Tables 1 and 2 are the project's original cohort-specific GOLD and Light coefficients, including the corrected NHANES GOLD signs. Applying them requires baseline age and the model's biomarkers in the listed units. Survival outcomes are needed to train a new model, not to apply these fixed coefficients.

## 1. Available coefficient sets

| Model | NHANES | UK Biobank | CHARLS | Reference |
|---|---|---|---|---|
| GOLD, fitted with Gompertz regression | Original coefficient set | Original coefficient set | — | Table 1 |
| Light, fitted with Gompertz regression | Original coefficient set | Original coefficient set | Original coefficient set | Table 2 |
| Cox-derived age | Reproducible refit of bundled NHANES4 | — | — | Table 3 |

The original CHARLS Light model uses the 2011 and 2015 waves. Table 3 labels a new reproducible example from the package; Tables 1 and 2 preserve the existing cohort models. A dash means this repository does not supply that combination.

The five original GOLD/Light sets are also installed with version **0.2.2 or later** as [published_coefficients.csv](inst/extdata/published_coefficients.csv). Its 37 rows record dataset, model, R term name, coefficient, input unit, transformation, and the decimal precision supplied by the project. The NHANES GOLD vector retains the supplied eight-decimal biomarker weights and five-decimal intercept. Other original coefficients retain the precision available in Tables 1 and 2. No extra digits are inferred.

## 2. GOLD: equation, units, and original coefficients

$$
\mathrm{GOLD\ BioAge}=A+b_0+\sum_{j=1}^{9}w_jX_j.
$$

$A$ is chronological age in years, $X_j$ is a biomarker value in the stated unit, $w_j$ is its weight in years per biomarker unit, and $b_0$ is the intercept in years. The coefficient of age is 1. Add the intercept once, and apply all nine biomarker terms.

**Table 1. Original GOLD cohort coefficients, displayed to four decimal places.**

| Variable | R column | Category | Input unit | NHANES | UKB |
|---|---|---|---|---:|---:|
| Age | `age` | | years | 1.0000 | 1.0000 |
| CREA | `creat` | Kidney | mg/dL | 5.2691 | 5.3832 |
| GLU | `glucose_mmol` | Metabolic | mmol/L | 0.5797 | 1.4168 |
| MCV | `mcv` | Erythrocyte | fL | 0.3389 | 0.4206 |
| RDW | `rdw` | Erythrocyte | % | 2.6445 | 3.3162 |
| ALB | `albumin` | Liver | g/dL | -4.7358 | -5.0793 |
| ALP | `alp` | Liver | U/L | 0.0260 | 0.0385 |
| LYM% | `lymph` | Immune | % | -0.2032 | -0.1899 |
| WBC | `wbc` | Immune | 1000 cells/uL | 0.4459 | 0.9120 |
| GGT | `ggt` | Liver | U/L | 0.0608 | 0.1007 |
| Intercept | `beta0` | | years | -53.6287 | -78.6519 |

For NHANES, the supplied higher-precision equation is:

$$
\begin{aligned}
\mathrm{GOLD}_{N}={}&A+5.26909180\,\mathrm{CREA}+0.57969768\,\mathrm{GLU}\\
&+0.33888497\,\mathrm{MCV}+2.64452031\,\mathrm{RDW}\\
&-4.73575224\,\mathrm{ALB}+0.02604091\,\mathrm{ALP}\\
&-0.20318306\,\mathrm{LYM\%}+0.44588830\,\mathrm{WBC}\\
&+0.06075152\,\mathrm{GGT}-53.62873.
\end{aligned}
$$

LYM% is negative and GGT is positive in the NHANES model. Enter 30 for 30% lymphocytes. The complete NHANES refit estimates an intercept near -53.62873289; the supplied fixed-model intercept is -53.62873, so a small precision difference is expected.

## 3. Light: equation, CRP transformation, and original coefficients

$$
\mathrm{Light\ BioAge}=A+b_0+w_C\,\mathrm{CREA}
+w_G\,\mathrm{GLU}+w_R\ln(1+\mathrm{CRP}).
$$

**Table 2. Original Light cohort coefficients. CRP units are specified before transformation.**

| Variable / R term | NHANES input unit | NHANES | UKB input unit | UKB | CHARLS input unit | CHARLS |
|---|---|---:|---|---:|---|---:|
| Age / `age` | years | 1.0000 | years | 1.0000 | years | 1.0000 |
| CREA / `creat` | mg/dL | 8.3313 | mg/dL | 9.0873 | mg/dL | 5.0875 |
| GLU / `glucose_mmol` | mmol/L | 0.8270 | mmol/L | 1.7311 | mmol/L | 1.5766 |
| ln(1 + CRP) / `log_crp` | CRP in mg/dL | 5.7305 | CRP in mg/L | 4.1689 | CRP in mg/L | 4.2983 |
| Intercept / `beta0` | years | -13.5298 | years | -20.1395 | years | -16.7562 |

The logarithm is natural. In R, use `log1p(CRP)`. Convert the raw CRP concentration to the required unit **before** taking the logarithm:

$$
0.3\ \mathrm{mg/dL}=3\ \mathrm{mg/L},\qquad
\ln(1+0.3)\ne\ln(1+3).
$$

`predict_gold_bioage()` can evaluate these named fixed Light coefficient vectors after you supply the correctly transformed `log_crp` column. It does not transform CRP or convert units automatically. A refit using bundled `NHANES4$lncrp` gives closely agreeing NHANES Light weights, with small differences from the original fixed reference documented in the NHANES tutorial. Its inverse transformation is consistent with CRP in mg/dL: `10 * expm1(lncrp)` reproduces the reported CRP summary in mg/L. This numerical check does not reconstruct the original raw-data preprocessing. For new measurements, construct `log_crp` from a raw concentration with a known unit.

## 4. Read the coefficients and calculate a GOLD score

```r
library(GOLDBioAge)
coefficient_file <- system.file("extdata", "published_coefficients.csv", package = "GOLDBioAge")
reference <- read.csv(coefficient_file, stringsAsFactors = FALSE)

# A local helper to select one of the five original coefficient sets.
get_beta <- function(dataset, model) {
  selected <- reference[reference$dataset == dataset & reference$model == model, ]
  if (!nrow(selected)) stop("No supplied coefficients for this dataset/model.")
  setNames(selected$coefficient, selected$term)
}

person <- data.frame(age = 50, albumin = 4, alp = 80, creat = 1,
                     glucose_mmol = 5, lymph = 30, mcv = 90, rdw = 13, wbc = 7, ggt = 25)
beta_nhanes <- get_beta("nhanes", "gold")
predict_gold_bioage(beta_nhanes, person)
# Expected: 51.10203967 years (51.1 with one decimal place).

beta_ukb <- get_beta("ukb", "gold")
predict_gold_bioage(beta_ukb, person)
# Expected: 50.7472 years.
```

This performs the weighted sum above. To see the same calculation explicitly:

```r
markers <- setdiff(names(beta_nhanes), "beta0")
manual_score <- as.vector(as.matrix(person[, markers, drop = FALSE]) %*%
                           beta_nhanes[markers] + beta_nhanes[["beta0"]])
stopifnot(isTRUE(all.equal(manual_score, predict_gold_bioage(beta_nhanes, person))))
```

For a cohort, pass a data frame with one participant per row instead of `person`. Missing predictor values return `NA` in their original row positions.

## 5. Calculate all three Light versions

```r
# Same concentration: NHANES requires 0.3 mg/dL; UKB/CHARLS require 3 mg/L.
light_nhanes <- data.frame(age = 50, creat = 1, glucose_mmol = 5, log_crp = log1p(0.3))
light_mg_l <- data.frame(age = 50, creat = 1, glucose_mmol = 5, log_crp = log1p(3))

light_scores <- c(
  NHANES = predict_gold_bioage(get_beta("nhanes", "light"), light_nhanes),
  UKB = predict_gold_bioage(get_beta("ukb", "light"), light_mg_l),
  CHARLS = predict_gold_bioage(get_beta("charls", "light"), light_mg_l)
)
round(light_scores, 8)
# NHANES: 50.43997842; UKB: 53.38262256; CHARLS: 52.17300905 years.
```

These are different fitted models applied to the same example concentrations. Different results are expected. Choose the model definition you intend to apply and keep its coefficients, input units, transformation, and intercept together.

## 6. A separately labeled NHANES Cox refit example

The following coefficients come from `cox_bioage(NHANES4, variables)` with default training calibration in version 0.2.1 or later. This uses 39,348 complete participants and 4,716 deaths, and illustrates refitting the bundled dataset. It is separate from the original GOLD/Light models above.

**Table 3. Reproducible NHANES Cox example, displayed to eight decimal places.**

| R term | Cox weight |
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

```r
data("NHANES4", package = "GOLDBioAge")
variables <- c("age", "albumin", "alp", "creat", "glucose_mmol",
               "lymph", "mcv", "rdw", "wbc", "ggt")
cox_fit <- cox_bioage(NHANES4, variables)
round(cox_fit$score_coefficients, 8)
NHANES4$Coxbioage <- predict_cox_bioage(cox_fit, NHANES4)
```

The biomarker units are the same as Table 1. Keep the fitted intercept with the weights. See the [Cox tutorial](docs/cox-bioage.md) for the derivation and `calibration = "legacy"` when reproducing the older score convention.

[Run the five original-model examples](examples/03_existing_coefficients.R) · [Refit NHANES](docs/nhanes.md) · [Simulated-data tutorial](docs/simulated-gompertz.md)

Reference: Hao M. et al. *Advanced Science* (2025), e01765. [doi:10.1002/advs.202501765](https://doi.org/10.1002/advs.202501765). Contact: haombio@gmail.com.
