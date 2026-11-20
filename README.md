# mfeargen

[![R-CMD-check](https://github.com/kennychechuanyu/mfeargen/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kennychechuanyu/mfeargen/actions/workflows/R-CMD-check.yaml)

`mfeargen` fits Bayesian mixture models to trial-level fear generalization data.
It combines Rescorla–Wagner learning with exponential similarity gradients and
estimates learning, generalization, response parameters, and posterior
probabilities for four latent groups: non-learners, overgeneralizers, physical
generalizers, and perceptual generalizers. Inspect the posterior probabilities
when interpreting group assignments.

The model is based on [Yu et al. (2023), *Communications Psychology*](https://doi.org/10.1038/s44271-023-00005-0).
The package includes simulation, parameter recovery, posterior predictive checks,
and models for repeated sessions.

## Installation

Use R 4.0 or later. Install [JAGS 4.3 or later](https://mcmc-jags.sourceforge.io/)
before installing the package. On Debian or Ubuntu, install the `jags` system
package. On macOS or Windows, use the JAGS installer for your R architecture.

```r
install.packages("remotes")
remotes::install_github("kennychechuanyu/mfeargen")

library(mfeargen)
rjags::jags.version()
```

If `rjags` cannot load JAGS, check that R and JAGS use the same architecture and
restart R after installing JAGS.

## Fit a model

```r
library(mfeargen)

simulated <- simulate_mfeargen(
  n_participants = 12,
  n_conditioning_trials = 24,
  n_generalization_trials = 2,
  model_type = "1cs",
  rating_scale = c(0, 100),
  seed = 42
)

fit <- fit_mfeargen(
  simulated,
  model_type = "1cs",
  rating_scale = c(0, 100),
  n.iter = 2000,
  n.burnin = 1000,
  n.chains = 3,
  n.thin = 1,
  parallel = FALSE,
  seed = 42
)

summary(fit)
fit$convergence$diagnostics
plot_mcmc_diagnostics(fit)
plot_group_probabilities(fit)
plot_generalization_gradient(fit)
plot_posterior_predictive(fit)
```

`n_generalization_trials` counts trials **per stimulus**. With the default ten
stimuli, this example creates 20 generalization trials per participant;
`n_conditioning_trials` counts all conditioning trials.

The short chains above demonstrate the workflow. For inference, increase retained
samples and inspect chain mixing. Diagnostics use rank-normalized split Rhat and
bulk and tail effective sample sizes. A single chain cannot establish convergence.
The defaults require Rhat below 1.01 and bulk and tail ESS of at least 100 per
chain. Constant stochastic parameters remain unassessed. Passing diagnostics does
not establish that the model describes the data well;
use the posterior predictive checks to examine that separately.

For differential conditioning, simulate with `model_type = "2cs"`, `cs_plus = 7`,
and `cs_minus = 3`, then fit with `model_type = "2cs"`. This design contains one
CS+ and one CS−.

## Predictions and uncertainty

```r
coef(fit)
fit$group_membership

# Uncertainty in expected ratings, in original input row order.
expected <- predict(fit, type = "response", interval = "credible",
                    format = "data.frame")
head(expected)

# Intervals for repeated observations include observation noise.
predicted <- predict(fit, type = "response", interval = "prediction",
                     format = "data.frame")

replicated <- posterior_predict(fit, n_draws = 100)
residuals(fit)
calculate_waic(fit)
```

Response tables retain participant IDs, trial numbers, observed ratings, and
missing observations. Expected ratings integrate over posterior uncertainty;
predictive intervals also include response noise. Predictions refer to fitted
participants and trial histories. Prediction for new participants or arbitrary
`newdata` is not implemented.

## Use your own data

Provide a data frame with one row per presented stimulus trial. Keep the actual
trial sequence and reinforcement history, including trials with missing ratings.
Remove inter-trial intervals before fitting.

| Column | Meaning |
| --- | --- |
| `participant` | Participant ID; numeric or text |
| `trials` | Positive integer trial order within each participant |
| `US_expect` | Numeric expectancy rating; `NA` is allowed |
| `Phy_size` | Physical stimulus value |
| `Per_size` | Measured perceptual stimulus value |
| `stimulus_type` | `"CS+"`, `"CS-"`, or `"GS"` |
| `US` | Actual outcome, `0` or `1`, for `1cs` |
| `USp`, `USm` | For `2cs`: CS+ targets use `0`/`1`; CS− targets use `-1` for inhibition, `0` off-cue, or `1` for reinforcement |
| `phase` | Optional `"conditioning"` or `"generalization"` |

When `phase` is supplied, learning updates occur during conditioning and are held
fixed during generalization. Without `phase`, CS trials update learning throughout
the sequence. Participants may have different trial counts. Duplicate trial
numbers within a participant are errors.

The size-column shortcut assumes physical and perceptual values share an origin
and units. If perceptual coordinates use another reference system, provide
precomputed physical and perceptual distances instead. For `1cs`, these are
`d_phy` and `d_per_old`, with the CS learning indicator `k`. For `2cs`, provide
`d_p_phy`, `d_m_phy`, `d_p_per`, and `d_m_per`, with `k_plus` and `k_minus`.
See `?prepare_mfeargen_data` for details. Missing covariates are not invented.

```r
trials <- read.csv("trials.csv", colClasses = c(participant = "character"),
                    stringsAsFactors = FALSE)
validate_mfeargen_data(trials, model_type = "1cs")
fit <- fit_mfeargen(trials, model_type = "1cs", rating_scale = c(0, 100),
                   observation_model = "censored", seed = 42)
```

Specify the instrument's known rating bounds. Inferring bounds from observed
responses changes the scale when participants do not use its full range. Fitting
uses an internal 0–10 scale; predictions and residuals use original rating units.

Use `observation_model = "censored"` when Gaussian responses are clipped to the
instrument's endpoints. It treats endpoint ratings as tail probabilities and
interior ratings as Gaussian densities. Use `"normal"` for an uncensored Gaussian
observation model; its predictive intervals may extend beyond the scale bounds.
Simulated data retain their observation model and scale as metadata. CSV files do
not retain those attributes, so specify both when fitting imported data.

An installed example reads a CSV and saves estimates, diagnostics, predictions,
plots, and the fitted model:

```r
source(system.file("examples", "analyze_trials.R", package = "mfeargen"))
fit <- analyze_trials("trials.csv", rating_scale = c(0, 100),
                      output_dir = "results", observation_model = "censored")
```

Predictive draws are generated from the stored parameter draws when needed,
which saves memory during fitting. Set `save_predictive = TRUE` to retain them
in the model object. Use `seed` in `posterior_predict()` or prediction intervals
for repeatable generated draws without changing the caller's random-number state.

Cross-validation stores fold scores by default. Set `keep_fits = TRUE` only when
you need the individual fold fits, which can use considerably more memory.

The example saves diagnostic results even when chains need more sampling. Inspect
those results before interpreting the saved parameter or group estimates.

## Priors and design checks

```r
priors <- set_priors_flexible(
  alpha_mu = prior_beta(shape1 = 2, shape2 = 2),
  pi = prior_dirichlet(alpha = c(1, 1, 1, 1))
)
fit_custom <- fit_mfeargen(simulated, priors = priors, seed = 42)
```

Custom priors specify population hyperparameters, rather than fixing individual
participant parameters. Prior templates are starting points. Check their
implications in your stimulus units and assess sensitivity to alternatives.
`run_parameter_recovery()` evaluates recovery for a simulated design; short
sampler runs only check that the procedure executes.

The repeated-session model supports `1cs` data with an additional `session`
column. It uses separate conditioning histories for each session and shares a
participant's latent group. Both observation models are supported. The current
model uses built-in population priors; custom hierarchical priors and `2cs`
requests produce an explicit error.

```r
# repeated contains participant, session, trials and the usual trial columns.
fit_repeated <- fit_mfeargen_hierarchical(
  repeated, rating_scale = c(0, 100), observation_model = "censored", seed = 42
)
predict(fit_repeated, type = "response", interval = "credible", format = "data.frame")
calculate_variance_decomposition(fit_repeated)
```

Repeated-session variance ratios describe the underlying normal parameter
scales, before group constraints and truncation. They are not empirical
intraclass correlations of the constrained parameters.

## Study data and citation

Use `citation("mfeargen")` for the package and model citations. The bundled
`exp1_data` and `exp2_data` retain the original variable names and measurements.
`replicate_paper_analysis()` reanalyzes these data with the current package model;
it does not claim an exact reproduction of the archived publication analysis.
Original study materials are available at [OSF](https://osf.io/sxjak/).

## Development

```r
install.packages(c("devtools", "roxygen2", "testthat"))
devtools::test()
devtools::check(args = "--no-manual")
```

CI installs JAGS and runs the package check on the release and previous-release
versions of R. The suite includes actual model fits, long trial histories,
missing observations, repeated sessions, prediction exports, and likelihood
checks. The suite verifies software behavior; assess convergence and model
adequacy for each analysis.

Report bugs with a reproducible example in the
[issue tracker](https://github.com/kennychechuanyu/mfeargen/issues).

## License

MIT. See [LICENSE](LICENSE).
