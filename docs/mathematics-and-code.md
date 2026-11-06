# GOLD and Cox BioAge: equations beside the R calculation

[Getting started](../README.md) · [中文说明](../README.zh-CN.md) · [Existing cohort coefficients](../Coefficients.md)

There are two workflows. To **calculate an existing score**, combine age and biomarkers with fixed coefficients. To **develop a new score**, first fit survival models in a training cohort. The resulting calculation has the same form, but the coefficients belong to that training cohort.

## 1. The calculation that every user needs

For participant $i$,

$$
\boxed{\mathrm{BioAge}_i=A_i+\sum_{j=1}^p w_jX_{ij}+b_0},\qquad
\mathrm{BioAgeDiff}_i=\mathrm{BioAge}_i-A_i.
$$

| Symbol | Meaning | In R |
|---|---|---|
| $A_i$ | Baseline chronological age, in years | `age` |
| $X_{ij}$ | Biomarker in the model's original unit and transformation | Named predictor column |
| $w_j$ | Effect expressed in years per biomarker unit | Named coefficient |
| $b_0$ | Fixed intercept in years | `beta0` |

Age has coefficient 1. BioAgeDiff is subtraction in years; it is not a regression residual. Positive values mean the model assigns an age above chronological age, without establishing an individual's lifespan or a treatment effect.

Here is a complete NHANES GOLD calculation using the supplied coefficients:

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

# The same weighted sum, written explicitly:
terms <- setdiff(names(beta), "beta0")
manual <- drop(as.matrix(person[, terms, drop = FALSE]) %*% beta[terms]) + beta[["beta0"]]
stopifnot(isTRUE(all.equal(score, unname(manual))))
```

Use creatinine in mg/dL, glucose in mmol/L, albumin in g/dL, lymphocytes and RDW in percent, WBC in 1000 cells/uL, MCV in fL, and ALP/GGT in U/L. For a whole dataset, provide one row per participant. Applying a fixed score requires neither follow-up time nor death status. Prediction retains row order and returns `NA` for missing predictors.

## 2. How Gompertz fitting produces those weights

Let $t$ be elapsed follow-up time and $A$ age at baseline. The [flexsurv Gompertz model](https://chjackson.github.io/flexsurv/reference/Gompertz.html) is

$$
h(t)=b e^{at},\quad
H(t)=\frac{b}{a}(e^{at}-1)\quad(a\ne0),\quad S(t)=e^{-H(t)}.
$$

At $a=0$, $H(t)=bt$. Positive $a$ means increasing hazard with follow-up. The follow-up shape $a$ and baseline-age coefficient are different parameters.

The package first estimates an age-only reference, then constrains the age coefficient in the full model:

$$
h_1(t\mid A)=b_1\exp(a_1t+\beta_A A),
$$
$$
h_2(t\mid A,X)=b_2\exp\left(a_2t+\beta_A A+\sum_j\theta_jX_j\right).
$$

The second model fixes **$\beta_A$**, while its follow-up shape remains estimated. GOLD matches instantaneous hazards at **$t=0$**, where $e^{a_1t}=e^{a_2t}=1$. Its baseline age equation therefore needs neither a common shape nor a follow-up time correction. Shape still participates in likelihood estimation, so changing a shape constraint can change the fitted coefficients. Effects are converted into years as $w_j=\theta_j/\beta_A$, requiring a finite, nonzero $\beta_A$.

At baseline, the uncorrected hazard-equivalent age would include $\ln(b_2/b_1)/\beta_A$. GOLD applies a training-sample correction

$$
\ln\gamma=\frac1n\sum_{i\in\mathrm{train}}
\ln\frac{h_1(0\mid A_i)}{h_2(0\mid A_i,X_i)}.
$$

Therefore the package's final intercept is

$$
b_0=\frac{\ln(\gamma b_2/b_1)}{\beta_A}
=\frac1n\sum_{i\in\mathrm{train}}\left(A_i-s_i\right),\qquad
s_i=A_i+\sum_j w_jX_{ij}.
$$

This explains why mean training BioAge equals mean training age. It also explains why the intercept is estimated once, then kept fixed for future data. The correction is a calibration convention, not proof of absolute-risk calibration in another population.

```r
data("NHANES4", package = "GOLDBioAge")
variables <- c("age", "albumin", "alp", "creat", "glucose_mmol",
               "lymph", "mcv", "rdw", "wbc", "ggt")
gold_fit <- gold_bioage(NHANES4, variables)
beta_A <- unname(gold_fit[["model1 coef"]]["age"])
theta <- gold_fit[["model2 coef"]][setdiff(variables, "age")]
weights <- c(age = 1, theta / beta_A)
training <- gold_fit$GOLDBioage
s <- drop(as.matrix(training[, names(weights), drop = FALSE]) %*% weights)
b0 <- mean(training$age - s)
stopifnot(abs(b0 - gold_fit$coefficients[["beta0"]]) < 1e-8)
stopifnot(abs(mean(training$GOLDbioage) - mean(training$age)) < 1e-8)
NHANES4$GOLD_refit <- predict_gold_bioage(gold_fit, NHANES4)
```

### Baseline evaluation and an optional constraint are separate choices

GOLD is defined at baseline. After fitting, its age equation contains no $a_1$ or $a_2$, even if the two shapes differ. The package therefore retains its original free-shape fitting and the supplied cohort coefficients.

If a separate analysis needs a common shape across follow-up, it can be fixed **before** the second regression. The direct `flexsurv` code is:

```r
initial <- gold_fit$models$full$res[, "est"]
initial[c("shape", "age")] <- gold_fit[["model1 coef"]][c("shape", "age")]
fixed <- match(c("shape", "age"), names(initial))
manual_shared <- flexsurv::flexsurvreg(
  survival::Surv(time, status) ~ age + albumin + alp + creat + glucose_mmol +
    lymph + mcv + rdw + wbc + ggt,
  data = training, dist = "gompertz", inits = initial, fixedpars = fixed
)
stopifnot(manual_shared$res["shape", "est"] == initial[["shape"]])
stopifnot(manual_shared$res["age", "est"] == initial[["age"]])
```

This demonstrates an optional constrained fit. It is not necessary for baseline GOLD calculation and can change the rate and biomarker coefficients. The follow-up shape affects estimation, even though it is absent from the final baseline scoring equation.

The complete-case fit has 39,348 participants and 4,716 deaths. It estimates an intercept near -53.62873289; the supplied fixed-model intercept is -53.62873. See the [NHANES tutorial](nhanes.md) for all fitted coefficients and the [simulation tutorial](simulated-gompertz.md) for a training/test example.

## 3. The Cox alternative

Cox fitting uses unspecified baseline hazards:

$$
h_1(t\mid A)=h_{01}(t)e^{\beta_A A},\qquad
h_2(t\mid A,X)=h_{02}(t)\exp\left(\beta_A A+\sum_j\theta_jX_j\right).
$$

The age-only estimate enters the full model through an offset. The package uses the same conversion to years and default training centering:

$$
w_j=\frac{\theta_j}{\beta_A},\quad
s_i=A_i+\sum_jw_jX_{ij},\quad
b_0=\operatorname{mean}_{\mathrm{train}}(A-s),\quad
\mathrm{CoxBioAge}_i=s_i+b_0.
$$

In standard Cox regression, the baseline hazard is unspecified: there is no Gompertz follow-up shape parameter to fix. The baseline-age coefficient can be fixed by an offset. The two fitted baseline hazards need not be equal. Equating absolute hazards at time $t$ would additionally require $\ln\{h_{02}(t)/h_{01}(t)\}/\beta_A$. The default score is a calibrated relative-risk-derived age, rather than an absolute mortality probability.

```r
cox_fit <- cox_bioage(NHANES4, variables)
cox_fit$score_coefficients             # Includes beta0.
NHANES4$Cox_refit <- predict_cox_bioage(cox_fit, NHANES4)
stopifnot(abs(mean(cox_fit$CoxBioage$Coxbioage) -
              mean(cox_fit$CoxBioage$age)) < 1e-8)

# A portable coefficient vector produces identical predictions.
stopifnot(isTRUE(all.equal(
  predict_cox_bioage(cox_fit$score_coefficients, NHANES4), NHANES4$Cox_refit
)))
```

The equivalent constrained Cox regression is explicit here:

```r
cox_training <- training
age_reference <- survival::coxph(survival::Surv(time, status) ~ age,
                                 data = cox_training, ties = "efron")
beta_age_cox <- unname(coef(age_reference)["age"])
cox_training$age_offset <- beta_age_cox * cox_training$age
manual_cox <- survival::coxph(
  survival::Surv(time, status) ~ offset(age_offset) + albumin + alp + creat +
    glucose_mmol + lymph + mcv + rdw + wbc + ggt,
  data = cox_training, ties = "efron"
)
stopifnot(isTRUE(all.equal(unname(coef(manual_cox)),
                           unname(cox_fit[["model2 coef"]][-1]), tolerance = 1e-8)))
```

The NHANES default Cox intercept is approximately -53.56354467. Its weights are listed separately in [Coefficients.md](../Coefficients.md). Use the [Cox tutorial](cox-bioage.md) for the full worked example and `calibration = "legacy"` for earlier score conventions.

## 4. Light BioAge and the CRP unit

Light uses the same weighted-sum calculation with three biomarkers:

$$
\mathrm{Light}=A+w_C\,\mathrm{CREA}+w_G\,\mathrm{GLU}
+w_R\ln(1+\mathrm{CRP})+b_0.
$$

The raw concentration inside the logarithm is **mg/dL for NHANES** and **mg/L for UKB/CHARLS**. Convert before transforming. Natural logarithms use `log1p()` in R.

```r
# The same concentration expressed in the two required units:
crp_mg_dl <- 0.3
crp_mg_l <- 10 * crp_mg_dl
nhanes_light <- 50 + 8.3313 * 1 + 0.8270 * 5 +
  5.7305 * log1p(crp_mg_dl) - 13.5298
ukb_light <- 50 + 9.0873 * 1 + 1.7311 * 5 +
  4.1689 * log1p(crp_mg_l) - 20.1395
c(NHANES = nhanes_light, UKB = ukb_light)
# 50.43997842 and 53.38262256 years.
```

The web calculator asks for mg/dL in both datasets and converts UKB CRP internally. The R prediction functions expect already prepared columns. See [all five supplied coefficient sets and runnable code](../Coefficients.md).

Reference: Hao et al., *Advanced Science* (2025), [doi:10.1002/advs.202501765](https://doi.org/10.1002/advs.202501765). Mathematical derivation and software documentation describe the implemented score; refitting a bundled cohort is a reproducibility example.
