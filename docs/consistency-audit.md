# Consistency audit: R package, formulas, calculator, and coefficient tables

Checked on 2026-11-06 for GOLDBioAge 0.2.2 and the corresponding calculator update. This records implementation agreement and reproducibility checks, with the limits of the available data.

[中文入口](../README.zh-CN.md) · [Coefficient definitions](../Coefficients.md) · [Equations and code](mathematics-and-code.md) · [Updates](package-updates.md)

## Results and corrections

| Component | Finding | Result or correction |
|---|---|---|
| NHANES GOLD | Earlier tables/calculator had incorrect LYM% and GGT signs | LYM% = -0.20318306; GGT = +0.06075152; verified against the bundled-data refit |
| NHANES/UKB GOLD | Age, nine biomarkers, and intercept checked together | Packaged reference, both Markdown tables, calculator, and worked R scores agree at stated precision |
| NHANES/UKB Light | CRP scale differs by model | NHANES uses ln(1 + CRP in mg/dL); UKB uses ln(1 + CRP in mg/L). Calculator accepts mg/dL in both and multiplies UKB CRP by 10 before the logarithm |
| Calculator numerics | Stable evaluation near zero CRP | Uses `Math.log1p()`; zero is valid; logarithm is natural |
| Calculator content | CRP handling was only explained in repository prose; citation was a preprint; background image was absent | Adds model-specific unit text and accessible selectors, cites the 2025 publication, removes the missing asset reference |
| Gompertz equations | Baseline age and elapsed time were conflated in introductory text | Separates baseline age coefficient from follow-up shape; includes the zero-shape limit |
| Baseline hazard matching | At t = 0 both follow-up shape factors equal 1, even with different fitted shapes | GOLD needs no shared-shape constraint or follow-up correction; the cumulative-risk comparison is explained separately |
| Cox calibration | Earlier centered linear predictions produced a different intercept convention | Default training centering fixed in 0.2.1; `legacy` preserves the earlier scores; frozen scoring retains the stored intercept |
| Optional screening | Earlier code could drop the first biomarker in Cox `glmnet` output | Fixed in 0.2.0; selector is biomarker-only Cox at `lambda.min`, followed by a constrained Gompertz model |
| Descriptive CRP tables | Log-transformed rows were labeled as concentrations | Labels now specify unitless ln(1 + CRP in mg/L); explain why this differs from the NHANES Light predictor |

The original supplied values remain the fixed reference: NHANES GOLD retains the supplied higher precision; UKB GOLD and all three Light models retain the supplied four-decimal values. A fresh fit can differ slightly from these finite-precision models. For example, NHANES Light refits `creat = 8.3313781`, whereas its original fixed coefficient is 8.3313. The audit does not silently replace the original model with a new fit.

## Reproduced NHANES results

The processed `NHANES4` file contains 54,945 rows and 30 columns (MD5: `d4cbf3c774014e68b1041fa2a9aeccf0`).

| Fit | Complete participants | Deaths | Intercept | Calibration |
|---|---:|---:|---:|---|
| Gompertz GOLD, nine biomarkers | 39,348 | 4,716 | -53.62873289 | Mean score equals mean training age |
| Cox, same nine biomarkers | 39,348 | 4,716 | -53.56354467 | Same training centering convention |
| Gompertz Light, three biomarkers | 38,001 | 6,373 | -13.5298442 | Same training centering convention |

The supplied fixed NHANES GOLD intercept is -53.62873. The difference from the refit is due to supplied precision. The two default nine-biomarker fits have mean training age approximately 49.49797 years. Cox's legacy intercept is approximately -49.49796686, yielding mean score 53.56354 years; weights are unchanged.

For Light, `10 * expm1(lncrp)` gives mean CRP 4.0718323 mg/L and SD 5.6681251. Re-logging those mg/L concentrations gives mean 1.2570271. These reproduce the descriptive values 4.07, 5.67, and 1.26 and support mg/dL inside the stored NHANES fitting transform. See the [NHANES tutorial](nhanes.md) and [executable script](../examples/02_nhanes.R) for the calculation.

## Fixed-score reference cases

Inputs: age 50 years, creatinine 1 mg/dL, glucose 5 mmol/L, MCV 90 fL, RDW 13%, albumin 4 g/dL, ALP 80 U/L, lymphocytes 30%, WBC 7 (1000 cells/uL), GGT 25 U/L; CRP 0.3 mg/dL (= 3 mg/L).

| Dataset | GOLD | Light | Calculator display (GOLD / Light) |
|---|---:|---:|---|
| NHANES | 51.10203967 | 50.43997842 | 51.1 / 50.4 |
| UKB | 50.74720000 | 53.38262256 | 50.7 / 53.4 |
| CHARLS | Not supplied | 52.17300905 | Not offered by calculator |

Scores are in years. These four calculator cases are checked against independently calculated R values. The cross-repository script checks all 37 supplied terms against both Markdown coefficient tables, all 32 calculator terms, and 800 varied inputs across the four implemented calculator models. Browser checks cover all four selections, displayed results, and the model-specific unit explanation.

## Reproducible checks

From the R repository, after installing the package and test dependencies:

```sh
Rscript -e 'testthat::test_local()'
Rscript examples/01_simulated_gompertz.R
Rscript examples/02_nhanes.R
Rscript examples/03_existing_coefficients.R
R CMD build --no-build-vignettes .
R CMD check --no-manual GOLDBioAge_0.2.2.tar.gz
```

From the calculator repository:

```sh
node --test tests/calculator.test.cjs
node tests/cross-repository.cjs /path/to/GOLDBioAge
```

The R tests cover coefficient references, fixed prediction, missing-row alignment, age ordering, constraints, optional selection, Cox calibration, and input validation. All analysis code blocks in the Markdown tutorials and the formula page's R example were executed with the installed package. Installation commands and commented optional file exports are excluded from those documentation runs.

Verification result: **117 R assertions passed** (0 failures, warnings, or skips); **R CMD check --no-manual: Status OK**. All three executable examples, **36 Markdown analysis code blocks**, the HTML R example, **10 calculator tests**, and the **800-input cross-repository comparison** passed. The formula page rendered 43 mathematical expressions without a MathJax error in the browser.

## Items that cannot be settled from the available records

- **UKB and CHARLS training:** individual-level training records are not included. Their fixed coefficients, units, tables, and calculations can be checked for agreement; an independent refit of those cohorts has not been performed.
- **Supplementary Table 2 sex counts:** calculator README lists CHARLS male/female counts 8,038 + 9,116 = 17,154 against a total of 17,163. RuLAS lists 839 + 948 = 1,787 against 1,785, and percentages do not reproduce those counts. The existing values are flagged rather than guessed or overwritten. Correcting them requires the underlying cohort summary or a verified correction from its source; they are not calculation inputs.
- **Raw NHANES processing:** the distributed file and score fits reproduce the documented aggregate results, but a complete raw-to-analysis preprocessing pipeline is not supplied. Numerical consistency is narrower than raw-data provenance or independent validation.

Method sources: [Hao et al., Advanced Science (2025)](https://doi.org/10.1002/advs.202501765), [flexsurv Gompertz parameterization](https://chjackson.github.io/flexsurv/reference/Gompertz.html), [survival Cox fitting](https://stat.ethz.ch/R-manual/R-devel/library/survival/html/coxph.html), and [Levine et al., Aging (2018)](https://doi.org/10.18632/aging.101414).
