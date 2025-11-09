#' @title Hierarchical Model Fitting Functions
#' @description Functions for fitting the three-level hierarchical mfeargen model
#' @name hierarchical_fitting
#' @importFrom stats runif rnorm rgamma rlnorm
NULL

#' Fit Hierarchical mfeargen Model
#'
#' @description
#' Fits a three-level hierarchical Bayesian mixture model that separates
#' population-level distributions, individual stable components, and 
#' session-specific variation.
#'
#' @param data Data frame with hierarchical structure including session information
#' @param model_type Character. "1cs" or "2cs"
#' @param n.iter Integer. Number of MCMC iterations (default: 20000)
#' @param n.burnin Integer. Number of burn-in iterations (default: 10000)
#' @param n.chains Integer. Number of MCMC chains (default: 3)
#' @param n.thin Integer. Thinning interval (default: 5)
#' @param lambda_criterion Numeric. Criterion for overgeneralization boundary
#' @param priors mfeargen_priors object for population hyperpriors
#' @param parallel Logical. Run chains in parallel (default: TRUE)
#' @param verbose Logical. Print progress messages (default: TRUE)
#' @param ... Additional arguments passed to jagsUI
#' @return mfeargen_hierarchical object
#' @export
#' @examples
#' \dontrun{
#' # Simulate hierarchical data
#' data <- simulate_mfeargen_hierarchical(
#'   n_participants = 20,
#'   n_sessions = c(2, 3, 2)  # Variable sessions per participant
#' )
#' 
#' # Fit hierarchical model
#' fit <- fit_mfeargen_hierarchical(
#'   data = data,
#'   model_type = "1cs",
#'   n.iter = 10000
#' )
#' 
#' # Analyze variance components
#' var_decomp <- calculate_variance_decomposition(fit)
#' plot_hierarchical_structure(fit)
#' }
fit_mfeargen_hierarchical <- function(data, 
                                     model_type = c("1cs", "2cs"),
                                     n.iter = 20000,
                                     n.burnin = 10000, 
                                     n.chains = 3,
                                     n.thin = 5,
                                     lambda_criterion = 0.7,
                                     priors = NULL,
                                     parallel = TRUE,
                                     verbose = TRUE,
                                     ...) {
  
  # Validate arguments
  model_type <- match.arg(model_type)
  
  if (!is.data.frame(data)) {
    stop("Data must be a data frame")
  }
  
  if (n.iter <= n.burnin) {
    stop("n.iter must be greater than n.burnin")
  }
  
  # Check for required hierarchical structure
  required_cols <- c("participant", "session", "trials", "US_expect")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Data missing required columns for hierarchical model: ", 
         paste(missing_cols, collapse = ", "))
  }
  
  if (verbose) message("Preparing hierarchical data for JAGS...")
  
  # Prepare hierarchical data structure
  jags_data <- prepare_hierarchical_data(data, model_type, lambda_criterion)
  
  # Set model file
  model_file <- system.file("jags", 
                           paste0("model_", model_type, "_hierarchical.txt"), 
                           package = "mfeargen")
  
  if (!file.exists(model_file)) {
    stop("Hierarchical JAGS model file not found: ", model_file)
  }
  
  if (verbose) message("Initializing hierarchical MCMC parameters...")
  
  # Initialize parameters for hierarchical model
  inits <- initialize_hierarchical_parameters(jags_data, n.chains)
  
  # Parameters to monitor
  parameters <- c(
    # Population parameters
    "mu_alpha", "mu_lambda", "mu_w0", "mu_w1",
    
    # Variance components
    "tau_alpha_sq", "tau_lambda_sq", "tau_w0_sq", "tau_w1_sq",
    "sigma_alpha_sq", "sigma_lambda_sq", "sigma_w0_sq", "sigma_w1_sq",
    
    # Stability metrics
    "stability_alpha", "stability_lambda", "stability_w0", "stability_w1",
    
    # Individual stable components
    "alpha_stable", "lambda_stable", "w0_stable", "w1_stable",
    
    # Group indicators and probabilities
    "gp", "pi", "sigma_obs"
  )
  
  # Run MCMC
  if (verbose) message("Running hierarchical MCMC sampling...")
  
  start_time <- Sys.time()
  
  jags_fit <- jagsUI::jags(
    data = jags_data,
    inits = inits,
    parameters.to.save = parameters,
    model.file = model_file,
    n.iter = n.iter,
    n.burnin = n.burnin,
    n.chains = n.chains,
    n.thin = n.thin,
    parallel = parallel,
    verbose = verbose,
    ...
  )
  
  end_time <- Sys.time()
  
  if (verbose) {
    message(sprintf("MCMC completed in %.1f minutes", 
                   as.numeric(end_time - start_time, units = "mins")))
  }
  
  # Check convergence
  if (verbose) message("Checking convergence...")
  convergence <- check_hierarchical_convergence(jags_fit)
  
  if (verbose) {
    message(sprintf("Convergence: %.1f%% of parameters have Rhat < 1.1", 
                   convergence$percent_converged))
  }
  
  # Extract variance decomposition
  if (verbose) message("Calculating variance decomposition...")
  var_decomp <- extract_variance_components(jags_fit)
  
  # Create results object
  results <- list(
    jags_fit = jags_fit,
    model_type = model_type,
    data = data,
    convergence = convergence,
    variance_components = var_decomp,
    data_info = list(
      n_participants = jags_data$Nparticipants,
      total_sessions = sum(jags_data$Nsessions),
      total_trials = sum(jags_data$Ntrials, na.rm = TRUE)
    ),
    population_parameters = extract_population_parameters(jags_fit),
    runtime = as.numeric(end_time - start_time, units = "mins")
  )
  
  class(results) <- c("mfeargen_hierarchical", "mfeargen", "list")
  
  if (verbose) message("Hierarchical model fitting completed successfully!")
  
  return(results)
}

#' Prepare Hierarchical Data for JAGS
#' @keywords internal
prepare_hierarchical_data <- function(data, model_type, lambda_criterion) {
  
  # Get unique participants and sessions
  participants <- sort(unique(data$participant))
  n_participants <- length(participants)
  
  # Initialize data structures
  max_sessions <- max(table(data$participant))
  max_trials <- max(table(paste(data$participant, data$session)))
  
  # Create arrays
  y <- array(NA, dim = c(n_participants, max_sessions, max_trials))
  k <- array(0, dim = c(n_participants, max_sessions, max_trials))  # Trial indicators
  r <- array(0, dim = c(n_participants, max_sessions, max_trials))  # US presence
  d_phy <- array(0, dim = c(n_participants, max_sessions, max_trials))
  d_per <- array(0, dim = c(n_participants, max_sessions, max_trials))
  
  # Session counts per participant
  nsessions <- numeric(n_participants)
  ntrials <- matrix(0, n_participants, max_sessions)
  
  # Fill arrays
  for (i in seq_along(participants)) {
    pid <- participants[i]
    p_data <- data[data$participant == pid, ]
    
    sessions <- sort(unique(p_data$session))
    nsessions[i] <- length(sessions)
    
    for (s in seq_along(sessions)) {
      sid <- sessions[s]
      s_data <- p_data[p_data$session == sid, ]
      n_trials_session <- nrow(s_data)
      ntrials[i, s] <- n_trials_session
      
      if (n_trials_session > 0) {
        # Fill trial data
        y[i, s, 1:n_trials_session] <- s_data$US_expect
        k[i, s, 1:n_trials_session] <- 1
        
        # Handle US presence
        if ("US" %in% names(s_data)) {
          r[i, s, 1:n_trials_session] <- s_data$US
        } else {
          # Infer from phase and stimulus type
          r[i, s, 1:n_trials_session] <- ifelse(
            s_data$phase == "conditioning" & s_data$stimulus_type == "CS+", 1, 0
          )
        }
        
        # Calculate distances
        if ("Phy_size" %in% names(s_data) && "cs_plus" %in% names(s_data)) {
          d_phy[i, s, 1:n_trials_session] <- abs(s_data$Phy_size - s_data$cs_plus)
        }
        
        if ("Per_size" %in% names(s_data)) {
          d_per[i, s, 1:n_trials_session] <- abs(s_data$Per_size - s_data$cs_plus)
        }
      }
    }
  }
  
  # Calculate lambda threshold
  if ("stimulus_range" %in% names(attributes(data))) {
    stimulus_range <- attr(data, "stimulus_range")
    max_distance <- stimulus_range[2] - stimulus_range[1]
    lambda_threshold <- -log(lambda_criterion) / max_distance
  } else {
    lambda_threshold <- 0.08  # Default fallback
  }
  
  # Rating scale bounds
  if ("rating_scale" %in% names(attributes(data))) {
    rating_scale <- attr(data, "rating_scale")
  } else {
    rating_scale <- c(min(data$US_expect, na.rm = TRUE),
                     max(data$US_expect, na.rm = TRUE))
  }
  
  # Create JAGS data list
  jags_data <- list(
    Nparticipants = n_participants,
    Nsessions = nsessions,
    Ntrials = ntrials,
    y = y,
    k = k,
    r = r,
    d_phy = d_phy,
    d_per = d_per,
    lambda_threshold = lambda_threshold,
    A = rating_scale[1],
    K = rating_scale[2],
    pi_alpha = c(1, 1, 1, 1)  # Dirichlet prior for group probabilities
  )
  
  return(jags_data)
}

#' Initialize Hierarchical MCMC Parameters
#' @keywords internal
initialize_hierarchical_parameters <- function(jags_data, n.chains) {
  
  inits_list <- list()
  
  for (chain in 1:n.chains) {
    # Set different seeds for each chain
    set.seed(chain * 123)
    
    # Population parameters
    mu_alpha <- runif(1, 0.2, 0.8)
    mu_lambda <- rlnorm(1, -2, 0.5)
    mu_w0 <- rnorm(1, 0, 2)
    mu_w1 <- rgamma(1, 2, 1)
    
    # Variance components
    tau_alpha_sq <- rgamma(1, 2, 4)
    tau_lambda_sq <- rgamma(1, 2, 4) 
    tau_w0_sq <- rgamma(1, 2, 2)
    tau_w1_sq <- rgamma(1, 2, 4)
    
    sigma_alpha_sq <- rgamma(1, 2, 8)
    sigma_lambda_sq <- rgamma(1, 2, 8)
    sigma_w0_sq <- rgamma(1, 2, 4)
    sigma_w1_sq <- rgamma(1, 2, 8)
    
    # Individual stable components
    n_participants <- jags_data$Nparticipants
    
    alpha_stable_raw <- rnorm(n_participants, mu_alpha, sqrt(tau_alpha_sq))
    alpha_stable_raw <- pmax(0.01, pmin(0.99, alpha_stable_raw))
    
    lambda_stable_raw <- rlnorm(n_participants, log(mu_lambda), sqrt(tau_lambda_sq))
    lambda_stable_raw <- pmax(0.001, lambda_stable_raw)
    
    w0_stable <- rnorm(n_participants, mu_w0, sqrt(tau_w0_sq))
    w1_stable_raw <- rgamma(n_participants, mu_w1^2 / tau_w1_sq, mu_w1 / tau_w1_sq)
    w1_stable_raw <- pmax(0.1, w1_stable_raw)
    
    # Group indicators
    gp <- sample(1:4, n_participants, replace = TRUE)
    
    # Observation error
    sigma_obs <- c(runif(1, 0.5, 1.5), runif(1, 2, 3))
    
    # Group probabilities
    pi <- as.numeric(rdirichlet(1, c(1, 1, 1, 1)))
    
    inits_list[[chain]] <- list(
      mu_alpha = mu_alpha,
      mu_lambda = mu_lambda, 
      mu_w0 = mu_w0,
      mu_w1 = mu_w1,
      tau_alpha_sq = tau_alpha_sq,
      tau_lambda_sq = tau_lambda_sq,
      tau_w0_sq = tau_w0_sq,
      tau_w1_sq = tau_w1_sq,
      sigma_alpha_sq = sigma_alpha_sq,
      sigma_lambda_sq = sigma_lambda_sq,
      sigma_w0_sq = sigma_w0_sq,
      sigma_w1_sq = sigma_w1_sq,
      alpha_stable_raw = alpha_stable_raw,
      lambda_stable_raw = lambda_stable_raw,
      w0_stable = w0_stable,
      w1_stable_raw = w1_stable_raw,
      gp = gp,
      sigma_obs = sigma_obs,
      pi = pi
    )
  }
  
  return(inits_list)
}

#' Helper function for Dirichlet random sampling
#' @keywords internal
rdirichlet <- function(n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x / as.vector(sm))
}

#' Check Hierarchical Model Convergence
#' @keywords internal
check_hierarchical_convergence <- function(jags_fit) {
  
  # Get Rhat values
  rhat_vals <- jags_fit$Rhat
  
  # Remove NA values
  rhat_vals <- rhat_vals[!is.na(rhat_vals)]
  
  # Calculate convergence metrics
  n_params <- length(rhat_vals)
  n_converged <- sum(rhat_vals < 1.1)
  percent_converged <- (n_converged / n_params) * 100
  
  # Identify problematic parameters
  problematic <- names(rhat_vals)[rhat_vals >= 1.1]
  
  results <- list(
    n_parameters = n_params,
    n_converged = n_converged,
    percent_converged = percent_converged,
    max_rhat = max(rhat_vals, na.rm = TRUE),
    problematic_parameters = problematic,
    rhat_summary = summary(rhat_vals)
  )
  
  return(results)
}

#' Extract Variance Components from MCMC
#' @keywords internal
extract_variance_components <- function(jags_fit) {
  
  # Get posterior means for variance components
  tau_alpha_sq <- mean(jags_fit$sims.list$tau_alpha_sq)
  tau_lambda_sq <- mean(jags_fit$sims.list$tau_lambda_sq)
  tau_w0_sq <- mean(jags_fit$sims.list$tau_w0_sq)
  tau_w1_sq <- mean(jags_fit$sims.list$tau_w1_sq)
  
  sigma_alpha_sq <- mean(jags_fit$sims.list$sigma_alpha_sq)
  sigma_lambda_sq <- mean(jags_fit$sims.list$sigma_lambda_sq)
  sigma_w0_sq <- mean(jags_fit$sims.list$sigma_w0_sq)
  sigma_w1_sq <- mean(jags_fit$sims.list$sigma_w1_sq)
  
  # Calculate stability metrics
  stability_alpha <- tau_alpha_sq / (tau_alpha_sq + sigma_alpha_sq)
  stability_lambda <- tau_lambda_sq / (tau_lambda_sq + sigma_lambda_sq)
  stability_w0 <- tau_w0_sq / (tau_w0_sq + sigma_w0_sq)
  stability_w1 <- tau_w1_sq / (tau_w1_sq + sigma_w1_sq)
  
  results <- list(
    between_variance = list(
      alpha = tau_alpha_sq,
      lambda = tau_lambda_sq,
      w0 = tau_w0_sq,
      w1 = tau_w1_sq
    ),
    within_variance = list(
      alpha = sigma_alpha_sq,
      lambda = sigma_lambda_sq,
      w0 = sigma_w0_sq,
      w1 = sigma_w1_sq
    ),
    stability = list(
      alpha = stability_alpha,
      lambda = stability_lambda,
      w0 = stability_w0,
      w1 = stability_w1
    )
  )
  
  return(results)
}

#' Extract Population Parameters
#' @keywords internal
extract_population_parameters <- function(jags_fit) {
  
  params <- list(
    mu_alpha = mean(jags_fit$sims.list$mu_alpha),
    mu_lambda = mean(jags_fit$sims.list$mu_lambda),
    mu_w0 = mean(jags_fit$sims.list$mu_w0),
    mu_w1 = mean(jags_fit$sims.list$mu_w1)
  )
  
  return(params)
}

#' S3 Print Method for Hierarchical Model
#' @param x A mfeargen_hierarchical object
#' @param ... Additional arguments (not used)
#' @export
print.mfeargen_hierarchical <- function(x, ...) {
  cat("Hierarchical mfeargen Model Results\n")
  cat("===================================\n\n")
  
  cat(sprintf("Model Type: %s\n", toupper(x$model_type)))
  cat(sprintf("Participants: %d\n", x$data_info$n_participants))
  cat(sprintf("Total Sessions: %d\n", x$data_info$total_sessions))
  cat(sprintf("Total Trials: %d\n", x$data_info$total_trials))
  cat(sprintf("Runtime: %.1f minutes\n\n", x$runtime))
  
  cat("Convergence:\n")
  cat(sprintf("  %.1f%% of parameters converged (Rhat < 1.1)\n", 
              x$convergence$percent_converged))
  cat(sprintf("  Max Rhat: %.3f\n\n", x$convergence$max_rhat))
  
  cat("Parameter Stability (ICC):\n")
  for (param in names(x$variance_components$stability)) {
    stability <- x$variance_components$stability[[param]]
    cat(sprintf("  %s: %.3f (%s)\n", 
                param, stability,
                ifelse(stability > 0.7, "High",
                      ifelse(stability > 0.3, "Moderate", "Low"))))
  }
  
  cat("\nUse summary() for detailed results\n")
  cat("Use plot() for diagnostic plots\n")
  cat("Use calculate_variance_decomposition() for full variance analysis\n")
}