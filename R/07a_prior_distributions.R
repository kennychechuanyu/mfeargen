#' @title Prior Distribution Specification Functions
#' @description Functions for specifying prior distributions in a flexible way
#' @name prior_distributions
#' @importFrom stats rgamma
NULL

#' Create a Normal Prior Distribution
#'
#' @description
#' Specifies a normal (Gaussian) prior distribution.
#'
#' @param mean Numeric. Mean of the distribution
#' @param sd Numeric. Standard deviation (must be positive)
#' @param truncate Numeric vector. Optional truncation bounds c(lower, upper)
#' @return A prior specification object
#' @export
#' @examples
#' # Standard normal prior
#' prior_normal(mean = 0, sd = 1)
#' 
#' # Truncated normal (positive values only)
#' prior_normal(mean = 5, sd = 2, truncate = c(0, Inf))
prior_normal <- function(mean = 0, sd = 1, truncate = NULL) {
  if (sd <= 0) stop("sd must be positive")
  
  prior <- list(
    dist = "norm",
    mean = mean,
    precision = 1 / sd^2,  # JAGS uses precision
    sd = sd,  # Keep for reference
    truncate = truncate
  )
  prior$distribution <- "normal"
  prior$parameters <- list(mean = mean, sd = sd, precision = 1 / sd^2, truncate = truncate)
  
  class(prior) <- c("mfeargen_prior", "prior_normal", "list")
  return(prior)
}


#' Create a Beta Prior Distribution
#'
#' @description
#' Specifies a beta prior distribution, typically used for probabilities.
#'
#' @param shape1 Numeric. First shape parameter (alpha, must be positive)
#' @param shape2 Numeric. Second shape parameter (beta, must be positive)
#' @param mode Numeric. Alternative parameterization by mode (optional)
#' @param concentration Numeric. Concentration around mode (optional)
#' @return A prior specification object
#' @export
#' @examples
#' # Uniform on [0,1]
#' prior_beta(shape1 = 1, shape2 = 1)
#' 
#' # Informative prior centered at 0.7
#' prior_beta(mode = 0.7, concentration = 10)
prior_beta <- function(shape1 = NULL, shape2 = NULL, mode = NULL, concentration = NULL) {
  
  # Allow parameterization by mode and concentration
  if (!is.null(mode) && !is.null(concentration)) {
    if (mode <= 0 || mode >= 1) stop("Mode must be between 0 and 1")
    if (concentration <= 2) stop("Concentration must be greater than 2")
    
    shape1 <- mode * (concentration - 2) + 1
    shape2 <- (1 - mode) * (concentration - 2) + 1
  }
  
  if (is.null(shape1) || is.null(shape2)) {
    stop("Must provide either shape1 and shape2, or mode and concentration")
  }
  
  if (shape1 <= 0 || shape2 <= 0) {
    stop("Shape parameters must be positive")
  }
  
  prior <- list(
    dist = "beta",
    a = shape1,
    b = shape2,
    mean = shape1 / (shape1 + shape2),
    mode = if (shape1 > 1 && shape2 > 1) (shape1 - 1) / (shape1 + shape2 - 2) else NA
  )
  prior$distribution <- "beta"
  prior$parameters <- list(shape1 = shape1, shape2 = shape2)
  
  class(prior) <- c("mfeargen_prior", "prior_beta", "list")
  return(prior)
}


#' Create a Gamma Prior Distribution
#'
#' @description
#' Specifies a gamma prior distribution, typically used for positive parameters.
#'
#' @param shape Numeric. Shape parameter (must be positive)
#' @param rate Numeric. Rate parameter (must be positive)
#' @param scale Numeric. Alternative scale parameterization (1/rate)
#' @param truncate Numeric vector. Optional truncation bounds c(lower, upper)
#' @return A prior specification object
#' @export
#' @examples
#' # Gamma with mean 2 and variance 1
#' prior_gamma(shape = 4, rate = 2)
#' 
#' # Using scale parameterization
#' prior_gamma(shape = 2, scale = 0.5)
prior_gamma <- function(shape, rate = NULL, scale = NULL, truncate = NULL) {
  
  if (!is.null(scale) && !is.null(rate)) {
    stop("Cannot specify both rate and scale")
  }
  
  if (!is.null(scale)) {
    rate <- 1 / scale
  }
  
  if (is.null(rate)) {
    stop("Must specify either rate or scale parameter")
  }
  
  if (shape <= 0 || rate <= 0) {
    stop("Shape and rate must be positive")
  }
  
  prior <- list(
    dist = "gamma",
    shape = shape,
    rate = rate,
    mean = shape / rate,
    variance = shape / rate^2,
    truncate = truncate
  )
  prior$distribution <- "gamma"
  prior$parameters <- list(shape = shape, rate = rate, truncate = truncate)
  
  class(prior) <- c("mfeargen_prior", "prior_gamma", "list")
  return(prior)
}


#' Create a Log-Normal Prior Distribution
#'
#' @description
#' Specifies a log-normal prior distribution.
#'
#' @param meanlog Numeric. Mean of the distribution on log scale
#' @param sdlog Numeric. SD of the distribution on log scale
#' @param mean Numeric. Alternative: mean on natural scale
#' @param sd Numeric. Alternative: SD on natural scale
#' @param truncate Numeric vector. Optional truncation bounds
#' @return A prior specification object
#' @export
#' @examples
#' # Log-normal with log-scale parameters
#' prior_lognormal(meanlog = 0, sdlog = 0.5)
#' 
#' # Log-normal with natural-scale parameters
#' prior_lognormal(mean = 2, sd = 1)
prior_lognormal <- function(meanlog = NULL, sdlog = NULL, mean = NULL, sd = NULL, truncate = NULL) {
  
  # Allow parameterization on natural scale
  if (!is.null(mean) && !is.null(sd)) {
    if (mean <= 0 || sd <= 0) stop("Mean and SD must be positive")
    
    # Convert to log scale
    cv <- sd / mean
    sdlog <- sqrt(log(1 + cv^2))
    meanlog <- log(mean) - sdlog^2 / 2
  }
  
  if (is.null(meanlog) || is.null(sdlog)) {
    stop("Must provide either meanlog and sdlog, or mean and sd")
  }
  
  if (sdlog <= 0) stop("SD on log scale must be positive")
  
  prior <- list(
    dist = "lnorm",
    mean = meanlog,
    precision = 1 / sdlog^2,  # JAGS parameterization
    sdlog = sdlog,
    truncate = truncate
  )
  prior$distribution <- "lognormal"
  prior$parameters <- list(meanlog = meanlog, sdlog = sdlog, truncate = truncate)
  
  class(prior) <- c("mfeargen_prior", "prior_lognormal", "list")
  return(prior)
}


#' Create a Uniform Prior Distribution
#'
#' @description
#' Specifies a uniform prior distribution.
#'
#' @param min Numeric. Lower bound
#' @param max Numeric. Upper bound
#' @return A prior specification object
#' @export
#' @examples
#' # Uniform on [0, 10]
#' prior_uniform(min = 0, max = 10)
prior_uniform <- function(min, max) {
  
  if (min >= max) stop("min must be less than max")
  
  prior <- list(
    dist = "unif",
    min = min,
    max = max,
    mean = (min + max) / 2,
    variance = (max - min)^2 / 12
  )
  prior$distribution <- "uniform"
  prior$parameters <- list(min = min, max = max)
  
  class(prior) <- c("mfeargen_prior", "prior_uniform", "list")
  return(prior)
}


#' Create a Dirichlet Prior Distribution
#'
#' @description
#' Specifies a Dirichlet prior distribution for probability vectors.
#'
#' @param alpha Numeric vector. Concentration parameters (all must be positive)
#' @param concentration Numeric vector. Alternative name for alpha parameter (optional)
#' @return A prior specification object
#' @export
#' @examples
#' # Uniform over 4-dimensional simplex
#' prior_dirichlet(alpha = c(1, 1, 1, 1))
#' 
#' # Informative prior favoring first category
#' prior_dirichlet(alpha = c(10, 2, 2, 2))
prior_dirichlet <- function(alpha = NULL, concentration = NULL) {
  
  if (is.null(alpha)) {
    alpha <- concentration
  }
  
  if (is.null(alpha)) {
    stop("Alpha/concentration vector must be provided")
  }
  
  if (any(alpha <= 0)) stop("All alpha values must be positive")
  
  prior <- list(
    dist = "dirichlet",
    alpha = alpha,
    mean = alpha / sum(alpha),
    k = length(alpha)
  )
  prior$distribution <- "dirichlet"
  prior$parameters <- list(alpha = alpha, concentration = alpha)
  
  class(prior) <- c("mfeargen_prior", "prior_dirichlet", "list")
  return(prior)
}


#' Create a Half-Cauchy Prior Distribution
#'
#' @description
#' Specifies a half-Cauchy prior distribution (positive values only).
#' Useful for scale parameters.
#'
#' @param scale Numeric. Scale parameter (must be positive)
#' @return A prior specification object
#' @export
#' @examples
#' # Weakly informative prior for standard deviation
#' prior_half_cauchy(scale = 2.5)
prior_half_cauchy <- function(scale) {
  
  if (scale <= 0) stop("Scale must be positive")
  
  prior <- list(
    dist = "half_cauchy",
    scale = scale
  )
  prior$distribution <- "half_cauchy"
  prior$parameters <- list(scale = scale)
  
  class(prior) <- c("mfeargen_prior", "prior_half_cauchy", "list")
  return(prior)
}


#' Create a Student-t Prior Distribution
#'
#' @description
#' Specifies a Student-t prior distribution.
#' More robust to outliers than normal distribution.
#'
#' @param df Numeric. Degrees of freedom (must be positive)
#' @param mean Numeric. Location parameter
#' @param scale Numeric. Scale parameter (must be positive)
#' @param truncate Numeric vector. Optional truncation bounds
#' @param location Numeric. Alternative name for mean parameter (optional)
#' @return A prior specification object
#' @export
#' @examples
#' # Heavy-tailed prior
#' prior_student_t(df = 3, mean = 0, scale = 1)
prior_student_t <- function(df, mean = 0, scale = 1, truncate = NULL, location = NULL) {
  if (!is.null(location)) {
    mean <- location
  }
  
  if (df <= 0) stop("Degrees of freedom must be positive")
  if (scale <= 0) stop("Scale must be positive")
  
  prior <- list(
    dist = "t",
    df = df,
    mean = mean,
    scale = scale,
    truncate = truncate
  )
  prior$distribution <- "student_t"
  prior$parameters <- list(df = df, mean = mean, scale = scale, truncate = truncate)
  
  class(prior) <- c("mfeargen_prior", "prior_student_t", "list")
  return(prior)
}


#' Create an Exponential Prior Distribution
#'
#' @description
#' Specifies an exponential prior distribution.
#'
#' @param rate Numeric. Rate parameter (must be positive)
#' @param truncate Numeric vector. Optional truncation bounds
#' @return A prior specification object
#' @export
#' @examples
#' # Exponential with mean 2
#' prior_exponential(rate = 0.5)
prior_exponential <- function(rate, truncate = NULL) {
  
  if (rate <= 0) stop("Rate must be positive")
  
  prior <- list(
    dist = "exp",
    rate = rate,
    mean = 1 / rate,
    variance = 1 / rate^2,
    truncate = truncate
  )
  prior$distribution <- "exponential"
  prior$parameters <- list(rate = rate, truncate = truncate)
  
  class(prior) <- c("mfeargen_prior", "prior_exponential", "list")
  return(prior)
}


#' Convert Prior to JAGS String
#'
#' @description
#' Converts a prior specification object to JAGS model syntax.
#'
#' @param prior A prior specification object
#' @param var_name Character. Variable name in JAGS model
#' @return Character string of JAGS syntax
#' @export
#' @keywords internal
prior_to_jags <- function(prior, var_name = NULL, param_name = NULL) {
  
  if (!inherits(prior, "mfeargen_prior")) {
    stop("Input must be a prior specification object")
  }
  if (is.null(var_name)) {
    var_name <- param_name
  }
  if (is.null(var_name)) {
    stop("var_name or param_name must be provided")
  }
  
  # Base distribution string
  jags_str <- switch(prior$dist,
    norm = sprintf("%s ~ dnorm(%.4f, %.4f)", var_name, prior$mean, prior$precision),
    beta = sprintf("%s ~ dbeta(%.4f, %.4f)", var_name, prior$a, prior$b),
    gamma = sprintf("%s ~ dgamma(%.4f, %.4f)", var_name, prior$shape, prior$rate),
    lnorm = sprintf("%s ~ dlnorm(%.4f, %.4f)", var_name, prior$mean, prior$precision),
    unif = sprintf("%s ~ dunif(%.4f, %.4f)", var_name, prior$min, prior$max),
    dirichlet = sprintf("%s ~ ddirch(c(%s))", var_name, paste(prior$alpha, collapse = ", ")),
    exp = sprintf("%s ~ dexp(%.4f)", var_name, prior$rate),
    stop("Unsupported distribution type: ", prior$dist)
  )
  
  # Add truncation if specified
  if (!is.null(prior$truncate)) {
    lower <- ifelse(is.infinite(prior$truncate[1]), "", prior$truncate[1])
    upper <- ifelse(is.infinite(prior$truncate[2]), "", prior$truncate[2])
    jags_str <- paste0(jags_str, sprintf(" T(%s,%s)", lower, upper))
  }
  
  return(jags_str)
}


#' Print Prior Distribution
#'
#' @description
#' Pretty prints a prior distribution specification.
#'
#' @param x A prior specification object
#' @param ... Additional arguments (ignored)
#' @return Invisible NULL
#' @export
print.mfeargen_prior <- function(x, ...) {
  
  cat("Prior Distribution:\n")
  
  # Distribution-specific formatting
  if (inherits(x, "prior_normal")) {
    cat(sprintf("  Normal(mean = %.3f, sd = %.3f)\n", x$mean, x$sd))
  } else if (inherits(x, "prior_beta")) {
    cat(sprintf("  Beta(a = %.3f, b = %.3f)\n", x$a, x$b))
    cat(sprintf("  Mean: %.3f\n", x$mean))
    if (!is.na(x$mode)) cat(sprintf("  Mode: %.3f\n", x$mode))
  } else if (inherits(x, "prior_gamma")) {
    cat(sprintf("  Gamma(shape = %.3f, rate = %.3f)\n", x$shape, x$rate))
    cat(sprintf("  Mean: %.3f, Variance: %.3f\n", x$mean, x$variance))
  } else if (inherits(x, "prior_lognormal")) {
    cat(sprintf("  LogNormal(meanlog = %.3f, sdlog = %.3f)\n", x$mean, x$sdlog))
  } else if (inherits(x, "prior_uniform")) {
    cat(sprintf("  Uniform(min = %.3f, max = %.3f)\n", x$min, x$max))
  } else if (inherits(x, "prior_dirichlet")) {
    cat(sprintf("  Dirichlet(alpha = [%s])\n", paste(round(x$alpha, 2), collapse = ", ")))
  } else if (inherits(x, "prior_exponential")) {
    cat(sprintf("  Exponential(rate = %.3f)\n", x$rate))
    cat(sprintf("  Mean: %.3f\n", x$mean))
  }
  
  if (!is.null(x$truncate)) {
    cat(sprintf("  Truncated to [%.3f, %.3f]\n", x$truncate[1], x$truncate[2]))
  }
  
  invisible(NULL)
}
