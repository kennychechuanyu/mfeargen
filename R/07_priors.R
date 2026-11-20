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
#' @param sigma1_min Lower bound for `sigma[1]` uniform prior (default: 1e-9)
#' @param sigma1_max Upper bound for `sigma[1]` uniform prior (default: 1.5)
#' @param sigma2_min Lower bound for `sigma[2]` uniform prior (default: 1.5)
#' @param sigma2_max Upper bound for `sigma[2]` uniform prior (default: 3)
#' @param pi_alpha Dirichlet alpha parameters for group probabilities (default: c(1,1,1,1))
#' @param w0_mu_mean Mean for w0_mu normal prior (default: 0)
#' @param w0_mu_sd SD for w0_mu normal prior (default: 10)
#' @param w0_sigma_shape Shape for w0_sigma scaled gamma prior (default: 2)
#' @param w0_sigma_rate Rate for w0_sigma scaled gamma prior (default: 1)
#' @param w1_a_shape Shape for w1_a scaled gamma prior (default: 2)
#' @param w1_a_rate Rate for w1_a scaled gamma prior (default: 1)
#' @param w1_b_shape Shape for w1_b scaled gamma prior (default: 2)
#' @param w1_b_rate Rate for w1_b scaled gamma prior (default: 1)
#' @param ... Named prior objects overriding the corresponding hyperpriors.
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
  if (length(overrides) && (is.null(names(overrides)) || any(!nzchar(names(overrides))))) {
    stop("Prior overrides must be named")
  }
  if (anyDuplicated(names(overrides))) stop("Prior overrides must have distinct names")
  supplied <- names(match.call(expand.dots = FALSE))
  mapping <- c(lambda_mu_mean = "lambda_mu", lambda_mu_sd = "lambda_mu",
               lambda_sigma_min = "lambda_sigma", lambda_sigma_max = "lambda_sigma",
               alpha_mu_a = "alpha_mu", alpha_mu_b = "alpha_mu",
               alpha_kappa_min = "alpha_kappa", alpha_kappa_max = "alpha_kappa",
               sigma1_min = "sigma1", sigma1_max = "sigma1",
               sigma2_min = "sigma2", sigma2_max = "sigma2", pi_alpha = "pi",
               w0_mu_mean = "w0_mu", w0_mu_sd = "w0_mu",
               w0_sigma_shape = "w0_sigma", w0_sigma_rate = "w0_sigma",
               w1_a_shape = "w1_a", w1_a_rate = "w1_a",
               w1_b_shape = "w1_b", w1_b_rate = "w1_b")
  override_names <- unname(mapping[intersect(names(mapping), supplied)])
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
                       lambda = "lambda_mu",
                       sigma = "sigma1",
                       sigma1 = "sigma1",
                       sigma2 = "sigma2",
                       nm)
      if (!target %in% names(priors)) stop("Unknown prior: ", nm)
      priors[[target]] <- value
      override_names <- c(override_names, nm)
    }
  }

  priors$alpha <- priors$alpha_mu
  priors$lambda <- priors$lambda_mu
  priors$sigma <- priors$sigma1
  attr(priors, "overrides") <- unique(override_names)

  class(priors) <- c("mfeargen_priors", "list")
  validate_priors(priors)
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
    stop("Priors must be created with set_priors() or set_priors_flexible()")
  }
  families <- list(lambda_mu = c("norm", "lnorm"), lambda_sigma = "unif",
                   alpha_mu = "beta", alpha_kappa = "unif", sigma1 = "unif",
                   sigma2 = "unif", pi = "dirichlet", w0_mu = "norm",
                   w0_sigma = c("gamma", "half_cauchy"),
                   w1_a = c("gamma", "half_cauchy"), w1_b = c("gamma", "half_cauchy"))
  for (nm in names(families)) {
    prior <- priors[[nm]]
    if (!inherits(prior, "mfeargen_prior") || !prior$dist %in% families[[nm]]) {
      stop(nm, " requires a ", paste(families[[nm]], collapse = " or "), " prior")
    }
    check_prior_truncation(prior$truncate)
    fields <- switch(prior$dist, norm = c("mean", "precision"),
                     lnorm = c("mean", "precision"), unif = c("min", "max"),
                     beta = c("a", "b"), gamma = c("shape", "rate"),
                     half_cauchy = "scale")
    for (field in fields) {
      check_prior_number(prior[[field]], paste0(nm, "$", field),
                         positive = field %in% c("precision", "a", "b", "shape", "rate", "scale"))
    }
    if (prior$dist == "unif" && (prior$min >= prior$max || prior$min < 0)) {
      stop(nm, " requires non-negative increasing bounds")
    }
    if (prior$dist == "dirichlet" &&
        (!is.numeric(prior$alpha) || length(prior$alpha) != 4L ||
         any(!is.finite(prior$alpha)) || any(prior$alpha <= 0))) {
      stop("pi alpha must have exactly four finite positive values")
    }
    if (!is.null(prior$truncate)) {
      bounds <- prior$truncate
      if (nm == "lambda_mu" && bounds[1] < 0) {
        stop("lambda_mu truncation lower bound must be non-negative")
      }
      if (nm %in% c("w0_sigma", "w1_a", "w1_b") && bounds[1] < 0) {
        stop(nm, " truncation lower bound must be non-negative")
      }
      if (nm %in% c("alpha_mu", "pi")) stop("Custom truncation is not supported for ", nm)
    }
  }
  TRUE
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
  cat("mfeargen prior specifications\n")
  nodes <- c(lambda_mu = "lambda_mu", lambda_sigma = "lambda_sigma",
             alpha_mu = "alpha_mu", alpha_kappa = "alpha_kappa", sigma1 = "sigma[1]",
             sigma2 = "sigma[2]", pi = "pi[1:4]", w0_mu = "w0_mu",
             w0_sigma = "w0_sigma", w1_a = "w1_a", w1_b = "w1_b")
  for (nm in names(nodes)) cat("  ", prior_to_jags(x[[nm]], nodes[[nm]]), "\n", sep = "")
  invisible(x)
}
