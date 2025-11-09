# mfeargen 1.0.0

## Initial Release (2025-11-09)

This is the first release of the mfeargen package, providing a comprehensive Bayesian framework for analyzing fear generalization data using latent mixture models.

### Main Features

* **Bayesian Mixture Modeling**: Identifies four distinct learning patterns in fear generalization:
  - Non-Learners (no CS-US learning)
  - Overgeneralizers (excessive fear generalization)
  - Physical Generalizers (generalization based on physical stimulus properties)
  - Perceptual Generalizers (generalization based on perceived stimulus properties)

* **Two-Phase Experimental Paradigm**:
  - Conditioning phase with Rescorla-Wagner reinforcement learning
  - Generalization test phase without reinforcement
  - Adaptive lambda threshold calculation based on experiment design

* **Comprehensive Analysis Tools**:
  - Model fitting with MCMC via JAGS
  - Parameter recovery analysis for model validation
  - Automatic rating scale detection and normalization
  - Support for both simple (1CS) and differential (2CS) conditioning

* **Prior Specification System**:
  - Multiple prior distribution functions
  - Pre-configured prior templates (weak, strong, skeptical, literature-based)
  - Flexible prior specification for advanced users

* **Visualization and Diagnostics**:
  - Generalization gradients by group
  - MCMC convergence diagnostics
  - Parameter recovery plots
  - Group membership visualization

* **S3 Methods**: Full suite of S3 methods for model objects:
  - print, summary, plot
  - coef, fitted, residuals, predict
  - logLik for model comparison

### Documentation

* Comprehensive vignettes for getting started and advanced usage
* All functions fully documented with examples
* Based on Yu et al. (2023) Communications Psychology paper

### System Requirements

* R >= 4.0.0
* JAGS >= 4.3.0 for Bayesian inference

### Authors

* Kenny Yu (Maintainer)
* Francis Tuerlinckx
* Wolf Vanpaemel
* Jonas Zaman

### Citation

Yu, K., Tuerlinckx, F., Vanpaemel, W., & Zaman, J. (2023). Humans display interindividual differences in the latent mechanisms underlying fear generalization behaviour. *Communications Psychology*, 1(1), 1-15.

### Repository

* OSF: https://osf.io/sxjak/