#' @title Prior Template Functions
#' @description Pre-configured prior sets for common modeling scenarios
#' @name prior_templates
NULL

#' Get Weakly Informative Priors
#'
#' @description
#' Returns weakly informative priors that provide minimal constraint
#' while ensuring computational stability. Suitable for exploratory analyses.
#'
#' @param model_type Character. "1cs" or "2cs". Both use the same template.
#' @return A mfeargen_priors object
#' @export
#' @examples
#' # Get default weakly informative priors
#' priors <- get_weak_priors()
#' 
#' \dontrun{
#' # Use in model fitting with real data
#' fit <- fit_mfeargen(data, priors = get_weak_priors())
#' }
get_weak_priors <- function(model_type = "1cs") {
  model_type <- match.arg(model_type, c("1cs", "2cs"))

  set_priors(
    # Lambda: weakly informative on log scale
    lambda_mu_mean = 0,
    lambda_mu_sd = 1,
    lambda_sigma_min = 1e-9,
    lambda_sigma_max = 2,
    
    # Alpha: uniform on probability scale
    alpha_mu_a = 1,
    alpha_mu_b = 1,
    alpha_kappa_min = 1,
    alpha_kappa_max = 20,
    
    # Sigma: wide range allowed
    sigma1_min = 1e-9,
    sigma1_max = 2,
    sigma2_min = 1.5,
    sigma2_max = 5,
    
    # Group probabilities: uniform
    pi_alpha = c(1, 1, 1, 1),
    
    # W0: centered at zero with wide variance
    w0_mu_mean = 0,
    w0_mu_sd = 10,
    w0_sigma_shape = 2,
    w0_sigma_rate = 0.5,
    
    # W1: weakly informative
    w1_a_shape = 2,
    w1_a_rate = 1,
    w1_b_shape = 2,
    w1_b_rate = 1
  )
}


#' Get Strongly Informative Priors
#'
#' @description
#' Constructs priors around user-supplied expectations for learning and
#' generalization rates. The confidence setting controls their concentration.
#'
#' @param expected_alpha Numeric. Expected learning rate (default: 0.5)
#' @param expected_lambda Numeric. Expected generalization parameter (default: 0.1)
#' @param confidence Character. "low", "medium", or "high" confidence in priors
#' @return A mfeargen_priors object
#' @export
#' @examples
#' # Strong priors with default expectations
#' priors <- get_strong_priors()
#' 
#' # Custom expectations with high confidence
#' priors <- get_strong_priors(
#'   expected_alpha = 0.7,
#'   expected_lambda = 0.15,
#'   confidence = "high"
#' )
get_strong_priors <- function(expected_alpha = 0.5, 
                             expected_lambda = 0.1,
                             confidence = c("medium", "low", "high")) {
  
  confidence <- match.arg(confidence)
  
  validate_scalar(expected_alpha, "expected_alpha", 0, 1, lower_open = TRUE, upper_open = TRUE)
  validate_scalar(expected_lambda, "expected_lambda", 0, lower_open = TRUE)

  # Calculate beta parameters for alpha_mu hyperprior
  # These control the distribution of the mean learning rate
  # More conservative values to avoid numerical issues
  alpha_prior_concentration <- switch(confidence,
    low = 2,
    medium = 4,
    high = 6
  )
  alpha_a <- expected_alpha * alpha_prior_concentration
  alpha_b <- (1 - expected_alpha) * alpha_prior_concentration
  
  # Calculate lambda parameters (on log scale)
  lambda_meanlog <- log(expected_lambda)
  lambda_sdlog <- switch(confidence,
    low = 0.5,
    medium = 0.25,
    high = 0.1
  )
  
  set_priors(
    # Lambda: informative based on expectation
    lambda_mu_mean = lambda_meanlog,
    lambda_mu_sd = lambda_sdlog,
    lambda_sigma_min = 0.01,
    lambda_sigma_max = switch(confidence, low = 0.5, medium = 0.3, high = 0.2),
    
    # Alpha: informative beta
    alpha_mu_a = alpha_a,
    alpha_mu_b = alpha_b,
    # Keep alpha_kappa in reasonable range
    alpha_kappa_min = switch(confidence, low = 1, medium = 2, high = 3),
    alpha_kappa_max = switch(confidence, low = 10, medium = 15, high = 20),
    
    # Sigma: tighter bounds
    sigma1_min = 0.1,
    sigma1_max = 1.2,
    sigma2_min = 1.5,
    sigma2_max = 3,
    
    # Group probabilities: slight preference for learners
    pi_alpha = switch(confidence,
      low = c(0.5, 1, 1.5, 1.5),
      medium = c(0.25, 1, 2, 2),
      high = c(0.1, 1, 3, 3)
    ),
    
    # W0: centered with moderate variance
    w0_mu_mean = 0,
    w0_mu_sd = switch(confidence, low = 5, medium = 3, high = 2),
    w0_sigma_shape = 3,
    w0_sigma_rate = 2,
    
    # W1: informative
    w1_a_shape = 3,
    w1_a_rate = 1.5,
    w1_b_shape = 3,
    w1_b_rate = 1.5
  )
}


#' Get Skeptical Priors
#'
#' @description
#' Returns skeptical priors that favor null effects and no group differences.
#' Useful for conservative analyses requiring strong evidence for effects.
#'
#' @return A mfeargen_priors object
#' @export
#' @examples
#' # Use skeptical priors for conservative analysis
#' priors <- get_skeptical_priors()
#' \dontrun{
#' # Use with real data
#' fit <- fit_mfeargen(data, priors = priors)
#' }
get_skeptical_priors <- function() {
  
  set_priors(
    # Lambda: favor broad generalization (small decay rates).
    lambda_mu_mean = log(0.01),
    lambda_mu_sd = 0.5,
    lambda_sigma_min = 1e-9,
    lambda_sigma_max = 0.1,
    
    # Alpha: favor values near 0 (slow learning)
    alpha_mu_a = 0.5,
    alpha_mu_b = 5,
    alpha_kappa_min = 5,
    alpha_kappa_max = 15,
    
    # Sigma: allow moderate noise
    sigma1_min = 0.5,
    sigma1_max = 1.5,
    sigma2_min = 1.5,
    sigma2_max = 2.5,
    
    # Group probabilities: favor non-learners
    pi_alpha = c(3, 1, 1, 1),
    
    # W0: tightly centered at zero
    w0_mu_mean = 0,
    w0_mu_sd = 2,
    w0_sigma_shape = 5,
    w0_sigma_rate = 5,
    
    # W1: favor small slopes
    w1_a_shape = 1,
    w1_a_rate = 2,
    w1_b_shape = 1,
    w1_b_rate = 2
  )
}




#' Get Priors for Parameter Recovery
#'
#' @description
#' Returns priors optimized for parameter recovery simulations.
#' These priors match typical simulation settings.
#'
#' @param true_values List. Optional true parameter values to center priors on
#' @return A mfeargen_priors object
#' @export
#' @examples
#' # Default recovery priors
#' priors <- get_recovery_priors()
#' 
#' # Recovery priors centered on true values
#' priors <- get_recovery_priors(
#'   true_values = list(alpha = 0.6, lambda = 0.12)
#' )
get_recovery_priors <- function(true_values = NULL) {
  
  if (!is.null(true_values)) {
    # Center priors on true values
    alpha_true <- true_values$alpha %||% 0.5
    lambda_true <- true_values$lambda %||% 0.1
    
    priors <- set_priors(
      lambda_mu_mean = log(lambda_true),
      lambda_mu_sd = 0.2,
      alpha_mu_a = alpha_true * 10,
      alpha_mu_b = (1 - alpha_true) * 10,
      alpha_kappa_min = 5,
      alpha_kappa_max = 15
    )
  } else {
    # Default recovery priors (mildly informative)
    priors <- set_priors(
      lambda_mu_mean = log(0.1),
      lambda_mu_sd = 0.3,
      lambda_sigma_min = 0.01,
      lambda_sigma_max = 0.5,
      alpha_mu_a = 2,
      alpha_mu_b = 2,
      alpha_kappa_min = 2,
      alpha_kappa_max = 20,
      sigma1_min = 0.1,
      sigma1_max = 1.5,
      sigma2_min = 1.5,
      sigma2_max = 3,
      pi_alpha = c(1, 1, 1, 1)
    )
  }
  
  return(priors)
}


#' Compare Prior Templates
#'
#' @description
#' Creates a comparison table of different prior templates.
#'
#' @param templates Character vector. Which templates to compare
#' @return A data frame comparing prior specifications
#' @export
#' @examples
#' # Compare all standard templates
#' compare_prior_templates()
#'
#' # Compare specific templates
#' compare_prior_templates(c("weak", "strong", "skeptical"))
compare_prior_templates <- function(templates = c("weak", "strong", "skeptical")) {

  # Get prior objects
  prior_list <- list()
  for (template in templates) {
    prior_list[[template]] <- switch(template,
      weak = get_weak_priors(),
      strong = get_strong_priors(),
      skeptical = get_skeptical_priors(),
      recovery = get_recovery_priors(),
      stop("Unknown template: ", template)
    )
  }
  
  comparison <- data.frame(
    Parameter = c("lambda_mu_mean", "lambda_mu_sd", 
                  "alpha_mu_a", "alpha_mu_b",
                  "pi_alpha[1]", "pi_alpha[2]", "pi_alpha[3]", "pi_alpha[4]",
                  "sigma1_min", "sigma1_max")
  )
  
  for (template in names(prior_list)) {
    prior <- prior_list[[template]]
    
    values <- c(
      prior$lambda_mu$mean,
      1/sqrt(prior$lambda_mu$precision),
      prior$alpha_mu$a,
      prior$alpha_mu$b,
      prior$pi$alpha[1],
      prior$pi$alpha[2],
      prior$pi$alpha[3],
      prior$pi$alpha[4],
      prior$sigma1$min,
      prior$sigma1$max
    )
    
    comparison[[template]] <- round(values, 3)
  }
  
  return(comparison)
}


#' Null coalescing operator
#' @noRd
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}