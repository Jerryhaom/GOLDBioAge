# Usage: reproduce NHANES GOLD and Light BioAge

[README](README.md) · [中文说明](README.zh-CN.md) · [Complete R script](examples/04_reproduce_nhanes.R) · [Original coefficient tables](Coefficients.md)

This guide fits **GOLD BioAge** and **Light BioAge** from the bundled `NHANES4` data and compares the results with the original NHANES coefficients. Run the R blocks below in order, or run `source("examples/04_reproduce_nhanes.R")` from a downloaded copy of this repository after installing the package.

Verified with GOLDBioAge **0.2.2**, R **4.5.3**, flexsurv **2.3.2**, and survival **3.8-6** on 12 November 2026:

| Model | Complete participants | Deaths | Agreement with original coefficients |
|---|---:|---:|---|
| GOLD BioAge | 39,348 | 4,716 | All nine biomarker weights agree to the supplied eight decimal places; the intercept agrees to five decimal places |
| Light BioAge | 38,001 | 6,373 | Closely reproduces the original model; the largest absolute coefficient difference is about 0.0000781 |

Both models use Gompertz regression through `gold_bioage()`. Light uses three biomarkers: creatinine, glucose, and transformed CRP. The reported BioAge coefficients are weights in years, available in `fit$coefficients`:

$$
\mathrm{BioAge}=\mathrm{age}+\sum_j w_jX_j+\beta_0.
$$

## 1. Install and load the data

```r
install.packages("remotes")
remotes::install_github("Jerryhaom/GOLDBioAge")
```

```r
library(GOLDBioAge)
data("NHANES4", package = "GOLDBioAge")
stopifnot(utils::packageVersion("GOLDBioAge") >= "0.2.2")
stopifnot(nrow(NHANES4) == 54945L)
```

Start each fit from the full `NHANES4` dataset. The function selects complete rows using that model's predictors, `time`, and `status`. Use the supplied data and units, with no train/test split, survey weighting, or biomarker selection for this reproduction. GOLD and Light have different complete samples; do not restrict Light to the rows used by GOLD.

## 2. Reproduce GOLD BioAge

```r
gold_variables <- c("age", "albumin", "alp", "creat", "glucose_mmol",
                    "lymph", "mcv", "rdw", "wbc", "ggt")
gold_fit <- gold_bioage(NHANES4, gold_variables, feature_selection = FALSE)
print(round(gold_fit$coefficients, 8))
gold_counts <- c(n = length(gold_fit$rows_used),
                 deaths = sum(NHANES4$status[gold_fit$rows_used]))
print(gold_counts)
stopifnot(gold_counts[["n"]] == 39348L, gold_counts[["deaths"]] == 4716L)
```

Expected fitted coefficients:

| Term | Input unit | Refitted coefficient |
|---|---|---:|
| `age` | years | 1.00000000 |
| `albumin` | g/dL | -4.73575224 |
| `alp` | U/L | 0.02604091 |
| `creat` | mg/dL | 5.26909180 |
| `glucose_mmol` | mmol/L | 0.57969768 |
| `lymph` | %, e.g. 30 for 30% | -0.20318306 |
| `mcv` | fL | 0.33888497 |
| `rdw` | % | 2.64452031 |
| `wbc` | 1000 cells/uL | 0.44588830 |
| `ggt` | U/L | 0.06075152 |
| `beta0` | intercept in years | -53.62873289 |

The original fixed intercept is **-53.62873**, the refitted value rounded to five decimal places. All nine biomarker weights agree to eight decimal places in the verified environment. Lymphocytes have a negative weight and GGT has a positive weight.

## 3. Reproduce Light BioAge

For NHANES Light, `log(CRP)` means the natural logarithm `ln(CRP + 1)`, with raw CRP in **mg/dL**. The bundled column `NHANES4$lncrp` is already transformed: use it directly, without taking another logarithm.

```r
light_variables <- c("age", "creat", "glucose_mmol", "lncrp")
light_fit <- gold_bioage(NHANES4, light_variables, feature_selection = FALSE)
print(round(light_fit$coefficients, 8))
light_counts <- c(n = length(light_fit$rows_used),
                  deaths = sum(NHANES4$status[light_fit$rows_used]))
print(light_counts)
stopifnot(light_counts[["n"]] == 38001L, light_counts[["deaths"]] == 6373L)
```

Expected fitted coefficients alongside the original fixed model:

| Term | Input unit / transformation | Refitted coefficient | Original fixed coefficient |
|---|---|---:|---:|
| `age` | years | 1.00000000 | 1.0000 |
| `creat` | mg/dL | 8.33137807 | 8.3313 |
| `glucose_mmol` | mmol/L | 0.82703148 | 0.8270 |
| `lncrp` | ln(CRP + 1), CRP in mg/dL | 5.73051692 | 5.7305 |
| `beta0` | intercept in years | -13.52984421 | -13.5298 |

Light closely reproduces the original coefficients, with a small last-digit difference, particularly for creatinine (**8.33137807** versus **8.3313**). The table preserves both values so users can distinguish a new fit from the original fixed model.

## 4. Check your fitted results against the original coefficients

The installed reference file is used only for comparison after fitting. Its CRP term is called `log_crp`; the corresponding fitting column in `NHANES4` is `lncrp`.

```r
reference <- read.csv(system.file("extdata", "published_coefficients.csv",
                                  package = "GOLDBioAge"))
compare_reference <- function(fit, model_name) {
  original <- reference[reference$dataset == "nhanes" &
                          reference$model == model_name, ]
  fitted_terms <- ifelse(original$term == "log_crp", "lncrp", original$term)
  refitted <- unname(fit$coefficients[fitted_terms])
  data.frame(term = original$term, original = original$coefficient,
             refitted = refitted, difference = refitted - original$coefficient)
}
gold_comparison <- compare_reference(gold_fit, "gold")
light_comparison <- compare_reference(light_fit, "light")
print(gold_comparison, digits = 10, row.names = FALSE)
print(light_comparison, digits = 10, row.names = FALSE)

stopifnot(all(abs(gold_comparison$difference) < 1e-5))
stopifnot(all(abs(light_comparison$difference) < 1e-4))
```

Both checks pass in the verified environment. They allow absolute coefficient differences below `1e-5` for GOLD and `1e-4` for Light, accommodating the supplied precision and the small Light refit difference. Small numerical changes can occur across optimizer or dependency versions; record your environment with `sessionInfo()`.

## 5. Calculate BioAge from the refitted models

Use `predict_gold_bioage()` for both fitted models. For a new raw NHANES CRP measurement, construct `lncrp` using `log1p(CRP_in_mg_dL)` to match the name and transformation used during fitting.

```r
person <- data.frame(age = 50, albumin = 4, alp = 80, creat = 1,
                     glucose_mmol = 5, lymph = 30, mcv = 90, rdw = 13, wbc = 7, ggt = 25,
                     lncrp = log1p(0.3))  # Raw CRP: 0.3 mg/dL.
scores <- c(GOLD = predict_gold_bioage(gold_fit, person),
            Light = predict_gold_bioage(light_fit, person))
print(round(scores, 4))
#    GOLD   Light
# 51.1020 50.4402

nhanes_scores <- data.frame(
  age = NHANES4$age,
  GOLD = predict_gold_bioage(gold_fit, NHANES4),
  Light = predict_gold_bioage(light_fit, NHANES4)
)
stopifnot(nrow(nhanes_scores) == nrow(NHANES4))
```

The individual uses the GOLD units listed above. Batch prediction keeps all 54,945 rows and returns `NA` where required predictors are missing. Keep the fitted coefficients and intercept fixed when scoring new participants. To use the original fixed coefficients instead, follow [the coefficient calculation examples](Coefficients.md).

```r
print(sessionInfo())
```

For the broader NHANES walkthrough, including Cox fitting, see [the NHANES tutorial](docs/nhanes.md).
