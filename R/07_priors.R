#' @title Prior Specification Functions for mfeargen
#' @description Functions for setting and validating custom prior distributions
#' @name priors
NULL

#' Set Prior Distributions for mfeargen Models
#'
#' @description
#' Creates a prior specification object for use with fit_mfeargen.
#' This function uses traditional parameter specification.
#' For using prior objects, see set_priors_flexible().
#' All unspecified priors will use package defaults.
#'
#' @param lambda_mu_mean Mean for lambda_mu lognormal prior (default: 0)
#' @param lambda_mu_sd SD for lambda_mu lognormal prior (default: 0.5)
#' @param lambda_sigma_min Lower bound for lambda_sigma uniform prior (default: 1e-9)
#' @param lambda_sigma_max Upper bound for lambda_sigma uniform prior (default: 1)
#' @param alpha_mu_a Alpha parameter for alpha_mu beta prior (default: 1)
#' @param alpha_mu_b Beta parameter for alpha_mu beta prior (default: 1)
#' @param alpha_kappa_min Lower bound for alpha_kappa uniform prior (default: 1)
#' @param alpha_kappa_max Upper bound for alpha_kappa uniform prior (default: 10)
#' @param sigma1_min Lower bound for sigma[1] uniform prior (default: 1e-9)
#' @param sigma1_max Upper bound for sigma[1] uniform prior (default: 1.5)
#' @param sigma2_min Lower bound for sigma[2] uniform prior (default: 1.5)
#' @param sigma2_max Upper bound for sigma[2] uniform prior (default: 3)
#' @param pi_alpha Dirichlet alpha parameters for group probabilities (default: c(1,1,1,1))
#' @param w0_mu_mean Mean for w0_mu normal prior (default: 0)
#' @param w0_mu_sd SD for w0_mu normal prior (default: 10)
#' @param w0_sigma_shape Shape for w0_sigma scaled gamma prior (default: 2)
#' @param w0_sigma_rate Rate for w0_sigma scaled gamma prior (default: 1)
#' @param w1_a_shape Shape for w1_a scaled gamma prior (default: 2)
#' @param w1_a_rate Rate for w1_a scaled gamma prior (default: 1)
#' @param w1_b_shape Shape for w1_b scaled gamma prior (default: 2)
#' @param w1_b_rate Rate for w1_b scaled gamma prior (default: 1)
#' @param ... Additional prior specifications (currently unused)
#'
#' @return A list of prior specifications for use with fit_mfeargen
#' @export
#' @seealso set_priors_flexible for using prior distribution objects
#' @examples
#' # Use all defaults
#' priors_default <- set_priors()
#' 
#' # Make alpha (learning rate) more informative
#' priors_informative <- set_priors(
#'   alpha_mu_a = 2,
#'   alpha_mu_b = 2  # Beta(2,2) centered at 0.5
#' )
#' 
#' # Stronger prior on lambda (generalization)
#' priors_lambda <- set_priors(
#'   lambda_mu_mean = -2,  # Smaller lambda values
#'   lambda_mu_sd = 0.25   # Less variance
#' )
#' 
#' # Equal group probabilities prior
#' priors_equal <- set_priors(
#'   pi_alpha = c(10, 10, 10, 10)  # Strong prior for equal groups
#' )
set_priors <- function(
  ...,
  # Lambda (generalization) priors
  lambda_mu_mean = 0,
  lambda_mu_sd = 0.5,
  lambda_sigma_min = 1e-9,
  lambda_sigma_max = 1,
  
  # Alpha (learning rate) priors  
  alpha_mu_a = 1,
  alpha_mu_b = 1,
  alpha_kappa_min = 1,
  alpha_kappa_max = 10,
  
  # Sigma (noise) priors
  sigma1_min = 1e-9,
  sigma1_max = 1.5,
  sigma2_min = 1.5,
  sigma2_max = 3,
  
  # Group probability priors
  pi_alpha = c(1, 1, 1, 1),
  
  # W0 (intercept) priors
  w0_mu_mean = 0,
  w0_mu_sd = 10,
  w0_sigma_shape = 2,
  w0_sigma_rate = 1,
  
  # W1 (slope) priors
  w1_a_shape = 2,
  w1_a_rate = 1,
  w1_b_shape = 2,
  w1_b_rate = 1
) {
  overrides <- list(...)
  override_names <- character(0)
  priors <- list(
    lambda_mu = prior_lognormal(meanlog = lambda_mu_mean, sdlog = lambda_mu_sd),
    lambda_sigma = prior_uniform(lambda_sigma_min, lambda_sigma_max),
    alpha_mu = prior_beta(shape1 = alpha_mu_a, shape2 = alpha_mu_b),
    alpha_kappa = prior_uniform(alpha_kappa_min, alpha_kappa_max),
    sigma1 = prior_uniform(sigma1_min, sigma1_max),
    sigma2 = prior_uniform(sigma2_min, sigma2_max),
    pi = prior_dirichlet(alpha = pi_alpha),
    w0_mu = prior_normal(mean = w0_mu_mean, sd = w0_mu_sd),
    w0_sigma = prior_gamma(shape = w0_sigma_shape, rate = w0_sigma_rate),
    w1_a = prior_gamma(shape = w1_a_shape, rate = w1_a_rate),
    w1_b = prior_gamma(shape = w1_b_shape, rate = w1_b_rate)
  )

  if (length(overrides) > 0) {
    for (nm in names(overrides)) {
      value <- overrides[[nm]]
      if (!inherits(value, "mfeargen_prior")) {
        stop("Override values must be mfeargen_prior objects")
      }
      target <- switch(nm,
                       alpha = "alpha_mu",
                       sigma1 = "sigma1",
                       sigma2 = "sigma2",
                       nm)
      priors[[target]] <- value
      priors[[nm]] <- value
      override_names <- c(override_names, nm)
    }
  }

  if (is.null(priors$alpha)) priors$alpha <- priors$alpha_mu
  if (is.null(priors$lambda)) priors$lambda <- priors$lambda_mu
  if (is.null(priors$sigma)) priors$sigma <- priors$sigma1
  attr(priors, "overrides") <- unique(override_names)

  class(priors) <- c("mfeargen_priors", "list")
  return(priors)
}


#' Validate Prior Specifications
#'
#' @description
#' Validates that prior specifications are properly formatted and contain valid values.
#'
#' @param priors A prior specification object from set_priors()
#' @return TRUE if valid, otherwise throws an error with details
#' @export
validate_priors <- function(priors) {
  
  if (!inherits(priors, "mfeargen_priors")) {
    stop("Priors must be created with set_priors() function")
  }
  
  # Check lambda_mu
  if (priors$lambda_mu$precision <= 0) {
    stop("lambda_mu precision must be positive")
  }
  
  # Check lambda_sigma
  if (priors$lambda_sigma$min >= priors$lambda_sigma$max) {
    stop("lambda_sigma min must be less than max")
  }
  if (priors$lambda_sigma$min < 0) {
    stop("lambda_sigma min must be non-negative")
  }
  
  # Check alpha_mu
  if (priors$alpha_mu$a <= 0 || priors$alpha_mu$b <= 0) {
    stop("alpha_mu beta parameters must be positive")
  }
  
  # Check alpha_kappa
  if (priors$alpha_kappa$min >= priors$alpha_kappa$max) {
    stop("alpha_kappa min must be less than max")
  }
  if (priors$alpha_kappa$min <= 0) {
    stop("alpha_kappa min must be positive")
  }
  
  # Check sigma priors
  if (!is.null(priors$sigma1$min) && !is.null(priors$sigma1$max)) {
    if (priors$sigma1$min >= priors$sigma1$max) {
      stop("sigma1 min must be less than max")
    }
  } else if (!is.null(priors$sigma1$shape) && !is.null(priors$sigma1$rate)) {
    if (priors$sigma1$shape <= 0 || priors$sigma1$rate <= 0) {
      stop("sigma1 shape and rate must be positive")
    }
  }
  if (!is.null(priors$sigma2$min) && !is.null(priors$sigma2$max)) {
    if (priors$sigma2$min >= priors$sigma2$max) {
      stop("sigma2 min must be less than max")
    }
  } else if (!is.null(priors$sigma2$shape) && !is.null(priors$sigma2$rate)) {
    if (priors$sigma2$shape <= 0 || priors$sigma2$rate <= 0) {
      stop("sigma2 shape and rate must be positive")
    }
  }
  if (!is.null(priors$sigma1$max) && !is.null(priors$sigma2$min)) {
    if (priors$sigma1$max > priors$sigma2$min) {
      warning("sigma1 max overlaps with sigma2 min - groups may not be well separated")
    }
  }
  
  # Check pi (Dirichlet)
  if (length(priors$pi$alpha) != 4) {
    stop("pi alpha must have exactly 4 values (one per group)")
  }
  if (any(priors$pi$alpha <= 0)) {
    stop("All pi alpha values must be positive")
  }
  
  # Check w0_mu
  if (priors$w0_mu$precision <= 0) {
    stop("w0_mu precision must be positive")
  }
  
  # Check gamma shape/rate parameters
  gamma_params <- c("w0_sigma", "w1_a", "w1_b")
  for (param in gamma_params) {
    if (priors[[param]]$shape <= 0 || priors[[param]]$rate <= 0) {
      stop(paste(param, "shape and rate must be positive"))
    }
  }
  
  return(TRUE)
}


#' Get Default Priors
#'
#' @description
#' Returns the default prior specifications used by the package.
#'
#' @return A mfeargen_priors object with default values
#' @export
#' @examples
#' defaults <- get_default_priors()
#' print(defaults)
get_default_priors <- function() {
  set_priors()
}


#' Print Prior Specifications
#'
#' @description
#' Pretty prints the prior specifications.
#'
#' @param x A mfeargen_priors object
#' @param ... Additional arguments (ignored)
#' @return Invisible NULL
#' @export
print.mfeargen_priors <- function(x, ...) {
  cat("mfeargen Prior Specifications\n")
  cat("=============================\n\n")
  
  cat("Lambda (Generalization) Hyperpriors:\n")
  cat(sprintf("  lambda_mu ~ LogNormal(%.2f, %.2f)\n", 
              x$lambda_mu$mean, 1/sqrt(x$lambda_mu$precision)))
  cat(sprintf("  lambda_sigma ~ Uniform(%.3f, %.2f)\n", 
              x$lambda_sigma$min, x$lambda_sigma$max))
  
  cat("\nAlpha (Learning Rate) Hyperpriors:\n")
  cat(sprintf("  alpha_mu ~ Beta(%.1f, %.1f)\n", 
              x$alpha_mu$a, x$alpha_mu$b))
  cat(sprintf("  alpha_kappa ~ Uniform(%.1f, %.1f)\n", 
              x$alpha_kappa$min, x$alpha_kappa$max))
  
  cat("\nSigma (Noise) Priors:\n")
  cat(sprintf("  sigma[1] ~ Uniform(%.3f, %.2f)\n", 
              x$sigma1$min, x$sigma1$max))
  cat(sprintf("  sigma[2] ~ Uniform(%.2f, %.2f)\n", 
              x$sigma2$min, x$sigma2$max))
  
  cat("\nGroup Probability Prior:\n")
  cat(sprintf("  pi ~ Dirichlet(%s)\n", 
              paste(x$pi$alpha, collapse = ", ")))
  
  cat("\nW0 (Intercept) Hyperpriors:\n")
  cat(sprintf("  w0_mu ~ Normal(%.1f, %.1f)\n", 
              x$w0_mu$mean, 1/sqrt(x$w0_mu$precision)))
  cat(sprintf("  w0_sigma ~ ScaledGamma(%.1f, %.1f)\n", 
              x$w0_sigma$shape, x$w0_sigma$rate))
  
  cat("\nW1 (Slope) Hyperpriors:\n")
  cat(sprintf("  w1_a ~ ScaledGamma(%.1f, %.1f)\n", 
              x$w1_a$shape, x$w1_a$rate))
  cat(sprintf("  w1_b ~ ScaledGamma(%.1f, %.1f)\n", 
              x$w1_b$shape, x$w1_b$rate))
  
  invisible(NULL)
}
