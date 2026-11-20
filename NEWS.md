# mfeargen 1.1.0

* Model ratings clipped to instrument endpoints with a censored Gaussian
  likelihood; simulation records its observation model for consistent fitting.
* Add response intervals, posterior predictive draws, and original-row prediction
  tables for export, including missing ratings and repeated sessions.
* Check chain mixing with rank-normalized split Rhat and bulk and tail effective
  sample sizes; distinguish unavailable diagnostics from passing checks.
* Plot fitted gradients and posterior predictive discrepancies alongside the
  observations.
* Validate complete prepared-data mappings and reject ambiguous input columns.
* Include a CSV analysis example that saves estimates, diagnostics, predictions,
  plots, and the fitted model.

# mfeargen 1.0.1

* Repair simulation, data preparation, and JAGS fitting for single-cue and
  differential conditioning designs.
* Preserve participant IDs and trial order, validate outcomes and covariates,
  and support missing ratings and unequal trial counts.
* Make rating bounds explicit and keep response summaries on the original scale.
* Calculate group assignments from posterior group probabilities.
* Replace placeholder summaries and plots with calculations from posterior draws.
* Avoid dynamic variance indexing so long trial sequences compile promptly.
* Add executable regression tests, JAGS integration checks, and package checks in CI.
* Correct the installation instructions, examples, and supported-model documentation.

# mfeargen 1.0.0

Initial package source release, including simulation, Bayesian mixture modeling,
prior specifications, and the experimental datasets from Yu et al. (2023).
