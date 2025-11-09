# mfeargen

[![CRAN status](https://www.r-pkg.org/badges/version/mfeargen)](https://CRAN.R-project.org/package=mfeargen)

## Overview

**mfeargen** provides a comprehensive Bayesian framework for analyzing fear generalization data using latent mixture models. The package implements a two-phase experimental paradigm (conditioning and generalization) and identifies four distinct learning patterns in fear generalization behavior.

Based on the computational model described in:
> Yu, K., Tuerlinckx, F., Vanpaemel, W., & Zaman, J. (2023). Humans display interindividual differences in the latent mechanisms underlying fear generalization behaviour. *Communications Psychology*, 1, 5. https://doi.org/10.1038/s44271-023-00005-0

## Features

- **Rescorla-Wagner Learning**: Implements reinforcement learning during conditioning phase
- **Shepard's Similarity Law**: Models stimulus generalization using exponential decay
- **Four latent groups**: Automatically identifies Non-Learners, Overgeneralizers, Physical, and Perceptual Generalizers
- **Hierarchical Bayesian Inference**: Full Bayesian treatment of individual and group parameters
- **Parameter recovery**: Validate model performance with simulated data
- **Perceptual noise modeling**: Control for individual differences in perception
- **Automatic scale detection**: Works with any rating scale
- **Comprehensive visualization**: Built-in plotting functions for all results
- **Flexible prior specification**: Support for custom priors and prior templates
- **1CS and 2CS paradigms**: Supports both single CS+ and differential CS+/CS- conditioning

## Installation

Install from GitHub:

```r
# install.packages("devtools")
devtools::install_github("kennychechuanyu/mfeargen")
```

### System Requirements

- R (>= 4.0.0)
- JAGS (>= 4.3.0) - [Download JAGS](https://sourceforge.net/projects/mcmc-jags/)

## Quick Start

```r
library(mfeargen)

# 1. Simulate fear generalization data (1CS paradigm)
data <- simulate_mfeargen(
  n_participants = 50,
  n_conditioning_trials = 24,      # Learning phase
  n_generalization_trials = 160,   # Test phase
  model_type = "1cs",              # Single CS+ (or "2cs" for CS+/CS-)
  rating_scale = c(0, 100)
)

# For 2CS paradigm (CS+ and CS-):
# data <- simulate_mfeargen(
#   model_type = "2cs",
#   cs_plus = 7,
#   cs_minus = 3,
#   n_participants = 50
# )

# 2. Fit the Bayesian mixture model
fit <- fit_mfeargen(
  data = data,
  model_type = "1cs",   # Must match simulation
  n.iter = 10000,       # MCMC iterations
  n.burnin = 5000,      # Burn-in period
  n.chains = 3          # Parallel chains
)

# 3. Examine results
summary(fit)          # Parameter estimates
plot(fit)             # Visualization dashboard

# 4. Extract group membership
groups <- fit$group_membership
table(groups$group_label)
```

## Advanced Usage

### Custom Priors

```r
# Create custom priors using flexible specification
custom_priors <- set_priors_flexible(
  lambda_mu = prior_lognormal(meanlog = 0, sdlog = 0.5),
  alpha_mu = prior_beta(shape1 = 1, shape2 = 1),
  pi = prior_dirichlet(alpha = c(1, 1, 1, 1))
)
fit <- fit_mfeargen(data, priors = custom_priors)
```

### Parameter Recovery Analysis

```r
# Validate model performance
recovery <- run_parameter_recovery(
  n_participants = 30,
  n_simulations = 10,
  model_type = "1cs"
)

# Visualize recovery quality
plot_parameter_recovery(recovery)

# Check recovery metrics
print(recovery)  # Shows correlations and accuracy
```

## Citation

If you use mfeargen in your research, please cite:

```bibtex
@article{yu2023humans,
  title={Humans display interindividual differences in the latent mechanisms underlying fear generalization behaviour},
  author={Yu, Kenny and Tuerlinckx, Francis and Vanpaemel, Wolf and Zaman, Jonas},
  journal={Communications Psychology},
  volume={1},
  pages={5},
  year={2023},
  doi={10.1038/s44271-023-00005-0}
}
```

## Acknowledgments

This work was supported by:
- FWO Research Project (G079520N)
- Research Foundation Flanders (FWO) Postdoctoral Fellowship (12P8623N)
- Special Research Funds (FWO, 1500620N)
- Research Fund of KU Leuven (C14/19/054)
- Flemish Supercomputer Center (VSC)

## Data and Code Availability

- **Original Paper Data & Code**: https://osf.io/sxjak/
- **Experiment 1 Pre-registration**: https://osf.io/b4ngs
- **Experiment 2 Pre-registration**: https://osf.io/t4bzs

## Getting Help

- **Documentation**: See package documentation with `help(package = "mfeargen")`
- **Contact**: kenny.yu@kuleuven.be

## License

MIT © 2023-2024 Kenny Yu, Francis Tuerlinckx, Wolf Vanpaemel, Jonas Zaman
