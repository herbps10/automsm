
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TargetedMSM <img src="man/figures/logo.png" align="right" height="140" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/herbps10/TargetedMSM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/herbps10/TargetedMSM/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

The `TargetedMSM` package provides targeted estimation for a general
class of Non-parametric Marginal Structural Models (NP-MSMs). An NP-MSM
summarizes a high-dimensional target functional (such as a conditional
average treatment effect) by projecting it onto a lower-dimensional
working model with respect to a loss function, without assuming the
working model is correctly specified.

The package supports arbitrary user-supplied working models and loss
functions. This flexibility is made possible by **automatic
differentiation** via the [`torch`](https://torch.mlverse.org/) package,
which computes the derivatives of the working model and loss function
needed to form the efficient influence function of the target
parameters.

## Supported estimands and estimators

The package currently supports three estimands, each summarizing a
different target functional. For each estimand, one or more estimators
are available: a one-step estimator, a Targeted Minimum Loss-Based
Estimator (TMLE), and a Generalized Bayesian TMLE. All three estimators
are asymptotically efficient; see the vignettes for further discussion
of their differences.

| Estimand | Function | One-step | TMLE | Bayesian TMLE |
|----|----|:--:|:--:|:--:|
| Conditional Average Treatment Effect (CATE) | `treatment_effect_modification` | ✓ | ✓ | ✓ |
| Categorical dose-response function | `categorical_dose_response` | ✓ | ✓ | ✓ |
| Longitudinal treatment effects | `longitudinal_treatment` | ✓ | ✓ |  |

## Installation

You can install the development version of TargetedMSM from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("herbps10/TargetedMSM")
```

## Example

The following example estimates a linear working model with a
squared-error loss function to summarize how a Conditional Average
Treatment Effect (CATE) varies with a treatment effect modifier:

``` r
library(TargetedMSM)

set.seed(10016)
data <- simulate_treatment_effect_modification(N = 500)

fit <- treatment_effect_modification(
  data,
  c("X1", "X2"),
  "A",
  "Y",
  formula = ~1 + X2,
  outcome_type = "continuous",
  learners_trt = c("SL.glm.interaction"),
  learners_outcome = c("SL.glm.interaction"),
  loss = loss_squared_error,
  working_model = working_model_linear
)
#> Loading required package: nnls

summary(fit)
#> 
#> ── Targeted MSM: treatment effect modification ─────────────────────────────────
#> n = 500
#> 
#> ── One-step estimator ──
#> 
#>         Term    Est    SE   2.5%  97.5%
#>  (Intercept)  0.944 0.062  0.823  1.066
#>           X2 -1.929 0.115 -2.155 -1.704
#> 
#> ── TMLE estimator ──
#> 
#>         Term    Est    SE   2.5%  97.5%
#>  (Intercept)  0.948 0.062  0.827  1.069
#>           X2 -1.939 0.115 -2.164 -1.714
```
