# mfeargen

[![R-CMD-check](https://github.com/kennychechuanyu/mfeargen/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kennychechuanyu/mfeargen/actions/workflows/R-CMD-check.yaml)

`mfeargen` fits Bayesian mixture models to trial-level fear generalization data.
It combines Rescorla–Wagner learning with exponential similarity gradients and
estimates posterior probabilities for four latent groups: non-learners,
overgeneralizers, physical generalizers, and perceptual generalizers. Group
membership is uncertain; use the posterior probabilities alongside assignments.

The model is based on [Yu et al. (2023), *Communications Psychology*](https://doi.org/10.1038/s44271-023-00005-0).
The package also provides simulation, parameter recovery, posterior summaries,
and diagnostic plots.

## Installation

Use R 4.0 or later. Install [JAGS 4.3 or later](https://mcmc-jags.sourceforge.io/)
before installing the package. On Debian or Ubuntu, install the `jags` system
package; on macOS or Windows, use the JAGS installer appropriate for your R
installation.

```r
install.packages("remotes")
remotes::install_github("kennychechuanyu/mfeargen")

library(mfeargen)
rjags::jags.version()
```

If `rjags` cannot load JAGS, check that R and JAGS use the same architecture and
restart R after installing JAGS.

## Simulate and fit a model

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
fit$convergence
fit$group_membership
plot_group_probabilities(fit)
plot_generalization_gradient(fit)
coef(fit)
fitted(fit)
```

`n_generalization_trials` is the number of trials **per stimulus**. With the
default ten stimuli, the example creates 20 generalization trials per participant.
`n_conditioning_trials` is the total conditioning trial count.

These short chains demonstrate the workflow. For inference, increase the number
of retained samples and check convergence and mixing with
`plot_mcmc_diagnostics(fit)`. A completed sampler run does not establish
convergence or reliable group classification.

For differential conditioning, simulate with `model_type = "2cs"`,
`cs_plus = 7`, and `cs_minus = 3`, then fit with `model_type = "2cs"`. The
standard `2cs` model represents one CS+ and one CS−. See
`?replicate_paper_analysis` for reanalysis of the bundled study datasets. The
current implementation is not a verified exact reproduction of the archived
analysis.

## Use your own data

Supply a data frame with one row per presented stimulus trial. Keep the actual
trial sequence and reinforcement outcomes, including trials with missing ratings.
Remove inter-trial intervals before fitting.

| Column | Meaning |
| --- | --- |
| `participant` | Participant ID; numeric or text |
| `trials` | Positive integer trial order within each participant |
| `US_expect` | Numeric expectancy rating; missing ratings may be `NA` |
| `Phy_size` | Physical stimulus value |
| `Per_size` | Measured perceptual stimulus value |
| `stimulus_type` | `"CS+"`, `"CS-"`, or `"GS"` |
| `US` | Actual outcome on each trial, `0` or `1`, for `1cs` |
| `USp`, `USm` | `2cs` learning targets: CS+ uses `0`/`1`; CS− uses `-1` for inhibition, `0` off-cue, or `1` for reinforcement |
| `phase` | Optional `"conditioning"` or `"generalization"` |

When `phase` is supplied, learning updates occur during conditioning and
associative strengths are held fixed during generalization. Without `phase`,
CS trials update learning throughout the sequence. Participants may have
different trial counts. Duplicate trial numbers within a participant are errors.
Physical and perceptual values must be observed; the package does not invent
perceptual measurements.

```r
trials <- read.csv("trials.csv")
validate_mfeargen_data(trials, model_type = "1cs")
fit <- fit_mfeargen(trials, model_type = "1cs", rating_scale = c(0, 100))
```

Specify the instrument's rating bounds with `rating_scale`. Inferring bounds
from observed responses can distort the scale when participants do not use its
full range. Simulated data retain their scale as metadata. Fitting uses a common
internal rating scale, while fitted responses and residuals use the original
scale.

For precomputed distances and CS indicators, see `?prepare_mfeargen_data`.

## Priors and model checks

```r
priors <- set_priors_flexible(
  alpha_mu = prior_beta(shape1 = 2, shape2 = 2),
  pi = prior_dirichlet(alpha = c(1, 1, 1, 1))
)
fit_custom <- fit_mfeargen(
  simulated,
  priors = priors,
  rating_scale = c(0, 100),
  seed = 42
)
```

Prior templates are starting points. Assess their implications on your stimulus
units and rating scale, and check sensitivity to reasonable alternatives.
`run_parameter_recovery()` can test whether an experimental design distinguishes
the parameters and groups you want to interpret. Recovery requires repeated
model fits and can take much longer than one analysis.

The repeated-session hierarchical model currently supports `1cs` data with an
additional `session` column. It estimates stable individual components and
session variation. It does not support custom priors or `2cs` designs; these
requests produce an explicit error.

## Development

```r
install.packages(c("devtools", "roxygen2", "testthat"))
devtools::test()
devtools::check(args = "--no-manual")
```

The test suite includes short JAGS runs for both conditioning designs. These
check installation and model execution; they are not parameter recovery studies.
CI runs the package check with JAGS on the release and previous-release versions
of R.

## Citation and data

Use `citation("mfeargen")` for the package and model citations.
Original study materials and data are available at [OSF](https://osf.io/sxjak/).
Bundled datasets are documented in `?exp1_data` and `?exp2_data`.

Report bugs with a reproducible example in the
[issue tracker](https://github.com/kennychechuanyu/mfeargen/issues).

## License

MIT. See [LICENSE](LICENSE) for copyright holders.
