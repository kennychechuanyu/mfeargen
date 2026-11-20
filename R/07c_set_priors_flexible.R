#' @title Flexible Prior Setting Function
#' @description Alternative function for setting priors using prior objects
#' @name set_priors_flexible
NULL

#' Set Priors Using Prior Objects
#'
#' @description
#' Alternative to set_priors() that accepts prior distribution objects.
#' Unspecified hyperpriors use the package defaults. Each parameter accepts
#' the distribution families listed below.
#'
#' @param lambda_mu Prior for lambda_mu (from prior_lognormal or prior_normal)
#' @param lambda_sigma Prior for lambda_sigma (from prior_uniform)
#' @param alpha_mu Prior for alpha_mu (from prior_beta)
#' @param alpha_kappa Prior for alpha_kappa (from prior_uniform)
#' @param sigma1 Prior for `sigma[1]` (from prior_uniform)
#' @param sigma2 Prior for `sigma[2]` (from prior_uniform)
#' @param pi Prior for group probabilities (from prior_dirichlet)
#' @param w0_mu Prior for w0_mu (from prior_normal)
#' @param w0_sigma Prior for w0_sigma (from prior_gamma)
#' @param w1_a Prior for w1_a (from prior_gamma)
#' @param w1_b Prior for w1_b (from prior_gamma)
#' @return A mfeargen_priors object
#' @export
#' @examples
#' # Using prior objects
#' priors <- set_priors_flexible(
#'   lambda_mu = prior_lognormal(mean = 0.1, sd = 0.05),
#'   alpha_mu = prior_beta(mode = 0.5, concentration = 10),
#'   pi = prior_dirichlet(alpha = c(1, 2, 2, 2))
#' )
#' 
#' # All parameters are optional - unspecified ones use defaults
#' priors <- set_priors_flexible(
#'   alpha_mu = prior_beta(shape1 = 3, shape2 = 3)
#' )
set_priors_flexible <- function(
  lambda_mu = NULL,
  lambda_sigma = NULL,
  alpha_mu = NULL,
  alpha_kappa = NULL,
  sigma1 = NULL,
  sigma2 = NULL,
  pi = NULL,
  w0_mu = NULL,
  w0_sigma = NULL,
  w1_a = NULL,
  w1_b = NULL
) {
  if (inherits(lambda_mu, "mfeargen_priors") &&
      is.null(lambda_sigma) && is.null(alpha_mu) && is.null(alpha_kappa) &&
      is.null(sigma1) && is.null(sigma2) && is.null(pi) &&
      is.null(w0_mu) && is.null(w0_sigma) && is.null(w1_a) && is.null(w1_b)) {
    return(lambda_mu)
  }
  
  # Get defaults
  defaults <- get_default_priors()
  override_names <- character(0)
  
  # Process lambda_mu
  if (!is.null(lambda_mu)) {
    if (!inherits(lambda_mu, "mfeargen_prior")) {
      stop("lambda_mu must be a prior object (use prior_lognormal or prior_normal)")
    }
    processed_lambda_mu <- lambda_mu
    override_names <- c(override_names, "lambda")
  } else {
    processed_lambda_mu <- defaults$lambda_mu
  }
  
  # Process lambda_sigma
  if (!is.null(lambda_sigma)) {
    if (!inherits(lambda_sigma, "mfeargen_prior")) {
      stop("lambda_sigma must be a prior object (use prior_uniform)")
    }
    processed_lambda_sigma <- lambda_sigma
    override_names <- c(override_names, "lambda_sigma")
  } else {
    processed_lambda_sigma <- defaults$lambda_sigma
  }
  
  # Process alpha_mu
  if (!is.null(alpha_mu)) {
    if (!inherits(alpha_mu, "mfeargen_prior")) {
      stop("alpha_mu must be a prior object (use prior_beta)")
    }
    processed_alpha_mu <- alpha_mu
    override_names <- c(override_names, "alpha")
  } else {
    processed_alpha_mu <- defaults$alpha_mu
  }
  
  # Process alpha_kappa
  if (!is.null(alpha_kappa)) {
    if (!inherits(alpha_kappa, "mfeargen_prior")) {
      stop("alpha_kappa must be a prior object (use prior_uniform)")
    }
    processed_alpha_kappa <- alpha_kappa
    override_names <- c(override_names, "alpha_kappa")
  } else {
    processed_alpha_kappa <- defaults$alpha_kappa
  }
  
  # Process sigma priors
  if (!is.null(sigma1)) {
    if (!inherits(sigma1, "mfeargen_prior")) {
      stop("sigma1 must be a prior object (use prior_uniform)")
    }
    processed_sigma1 <- sigma1
    override_names <- c(override_names, "sigma1")
  } else {
    processed_sigma1 <- defaults$sigma1
  }

  if (!is.null(sigma2)) {
    if (!inherits(sigma2, "mfeargen_prior")) {
      stop("sigma2 must be a prior object (use prior_uniform)")
    }
    processed_sigma2 <- sigma2
    override_names <- c(override_names, "sigma2")
  } else {
    processed_sigma2 <- defaults$sigma2
  }
  
  # Process pi (Dirichlet)
  if (!is.null(pi)) {
    if (!inherits(pi, "mfeargen_prior")) {
      stop("pi must be a prior object (use prior_dirichlet)")
    }
    processed_pi <- pi
    override_names <- c(override_names, "pi")
  } else {
    processed_pi <- defaults$pi
  }
  
  # Process w0 priors
  if (!is.null(w0_mu)) {
    if (!inherits(w0_mu, "mfeargen_prior")) {
      stop("w0_mu must be a prior object (use prior_normal)")
    }
    processed_w0_mu <- w0_mu
    override_names <- c(override_names, "w0_mu")
  } else {
    processed_w0_mu <- defaults$w0_mu
  }
  
  if (!is.null(w0_sigma)) {
    if (!inherits(w0_sigma, "mfeargen_prior")) {
      stop("w0_sigma must be a prior object (use prior_gamma)")
    }
    processed_w0_sigma <- w0_sigma
    override_names <- c(override_names, "w0_sigma")
  } else {
    processed_w0_sigma <- defaults$w0_sigma
  }
  
  # Process w1 priors
  if (!is.null(w1_a)) {
    if (!inherits(w1_a, "mfeargen_prior")) {
      stop("w1_a must be a prior object (use prior_gamma)")
    }
    processed_w1_a <- w1_a
    override_names <- c(override_names, "w1_a")
  } else {
    processed_w1_a <- defaults$w1_a
  }
  
  if (!is.null(w1_b)) {
    if (!inherits(w1_b, "mfeargen_prior")) {
      stop("w1_b must be a prior object (use prior_gamma)")
    }
    processed_w1_b <- w1_b
    override_names <- c(override_names, "w1_b")
  } else {
    processed_w1_b <- defaults$w1_b
  }
  
  # Create prior specification list
  priors <- list(
    lambda_mu = processed_lambda_mu,
    lambda_sigma = processed_lambda_sigma,
    alpha_mu = processed_alpha_mu,
    alpha_kappa = processed_alpha_kappa,
    sigma1 = processed_sigma1,
    sigma2 = processed_sigma2,
    pi = processed_pi,
    w0_mu = processed_w0_mu,
    w0_sigma = processed_w0_sigma,
    w1_a = processed_w1_a,
    w1_b = processed_w1_b
  )
  priors$alpha <- priors$alpha_mu
  priors$lambda <- priors$lambda_mu
  priors$sigma <- priors$sigma1
  attr(priors, "overrides") <- unique(override_names)

  class(priors) <- c("mfeargen_priors", "list")
  validate_priors(priors)
  return(priors)
}


#' Combine Prior Specifications
#'
#' @description
#' Combines multiple prior specifications, with later ones overriding earlier ones.
#' Useful for starting with a template and modifying specific priors.
#'
#' @param ... Multiple prior specifications to combine
#' @return A combined mfeargen_priors object
#' @export
#' @examples
#' # Start with weak priors and modify alpha
#' priors <- combine_priors(
#'   get_weak_priors(),
#'   set_priors_flexible(
#'     alpha_mu = prior_beta(mode = 0.7, concentration = 10)
#'   )
#' )
combine_priors <- function(...) {
  prior_list <- list(...)
  
  if (length(prior_list) == 0) {
    stop("At least one prior specification must be provided")
  }
  
  # Start with first prior set
  combined <- prior_list[[1]]
  
  if (!inherits(combined, "mfeargen_priors")) {
    stop("All arguments must be mfeargen_priors objects")
  }
  combined_overrides <- attr(combined, "overrides")
  if (is.null(combined_overrides)) combined_overrides <- character(0)
  
  # Override with subsequent prior sets
  if (length(prior_list) > 1) {
    for (i in 2:length(prior_list)) {
      current <- prior_list[[i]]
      
      if (!inherits(current, "mfeargen_priors")) {
        stop("All arguments must be mfeargen_priors objects")
      }
      override_tags <- attr(current, "overrides")
      if (is.null(override_tags) || length(override_tags) == 0) {
        next
      }
      for (nm in override_tags) {
        value <- current[[nm]]
        if (is.null(value)) next
        target <- switch(nm,
                         alpha = "alpha_mu",
                         alpha_mu = "alpha_mu",
                         lambda = "lambda_mu",
                         lambda_mu = "lambda_mu",
                         sigma1 = "sigma1",
                         sigma2 = "sigma2",
                         nm)
        combined[[target]] <- value
        if (nm %in% c("alpha", "alpha_mu")) {
          combined$alpha <- value
        } else if (nm %in% c("lambda", "lambda_mu")) {
          combined$lambda <- value
        } else if (nm %in% c("sigma", "sigma1", "sigma2")) {
          combined$sigma <- value
          if (nm %in% c("sigma1")) {
            combined$sigma1 <- value
          }
          if (nm %in% c("sigma2")) {
            combined$sigma2 <- value
          }
        } else {
          combined[[nm]] <- value
        }
      }
      combined_overrides <- unique(c(combined_overrides, override_tags))
    }
  }
  
  combined$alpha <- combined$alpha_mu
  combined$lambda <- combined$lambda_mu
  combined$sigma <- combined$sigma1
  attr(combined, "overrides") <- combined_overrides
  class(combined) <- c("mfeargen_priors", "list")
  validate_priors(combined)
  return(combined)
}
