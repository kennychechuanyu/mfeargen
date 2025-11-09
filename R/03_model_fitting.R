#' @title Model Fitting Functions for mfeargen
#' @description Core functions for fitting Bayesian mixture models
#' @name model_fitting
#' @importFrom stats runif rnorm rgamma
NULL

#' Fit mfeargen Bayesian Model
#'
#' @description
#' Fits a Bayesian mixture model to fear generalization data using JAGS.
#' Supports both 1CS (simple conditioning) and 2CS (differential conditioning) designs.
#'
#' @param data Data frame with required columns: participant, trials, US_expect
#' @param model_type Character. "1cs" or "2cs" 
#' @param n.iter Integer. Number of MCMC iterations (default: 10000)
#' @param n.burnin Integer. Number of burn-in iterations (default: 5000)
#' @param n.chains Integer. Number of MCMC chains (default: 3)
#' @param n.thin Integer. Thinning interval (default: 2)
#' @param lambda_criterion Numeric. Criterion for overgeneralization boundary (default: 0.7)
#' @param priors mfeargen_priors object. Custom priors from set_priors() (optional)
#' @param prepared_data List. Pre-prepared JAGS data from prepare_mfeargen_data() (optional)
#' @param parallel Logical. Run chains in parallel (default: TRUE)
#' @param verbose Logical. Print progress messages (default: TRUE)
#' @param ... Additional arguments passed to jagsUI
#' 
#' @return mfeargen object containing:
#'   - jags_fit: MCMC output from JAGS
#'   - model_type: Model specification used
#'   - data: Original data
#'   - convergence: Convergence diagnostics
#'   - group_membership: Predicted group assignments
#'   
#' @export
#' @import jagsUI
#' @examples
#' \dontrun{
#' # Simulate data
#' data <- simulate_mfeargen(n_participants = 20)
#' 
#' # Fit model
#' fit <- fit_mfeargen(data, model_type = "1cs", n.iter = 5000)
#' 
#' # View results
#' summary(fit)
#' plot(fit)
#' }
fit_mfeargen <- function(data, 
                        model_type = c("1cs", "2cs"), 
                        n.iter = 10000,
                        n.burnin = 5000,
                        n.chains = 3,
                        n.thin = 2,
                        lambda_criterion = 0.7,
                        priors = NULL,
                        prepared_data = NULL,
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
  
  # Determine if using flexible priors
  use_flexible_priors <- !is.null(priors)
  
  # Prepare data for JAGS
  if (!is.null(prepared_data)) {
    if (!is.list(prepared_data)) {
      stop("prepared_data must be a list created by prepare_mfeargen_data()")
    }
    jags_data <- prepared_data
    if (verbose) message("Using supplied prepared data for JAGS.")
  } else {
    if (verbose) message("Preparing data for JAGS...")
    jags_data <- prepare_mfeargen_data(data, model_type, lambda_criterion, 
                                      rating_scale = NULL,
                                      priors = priors,
                                      use_flexible_priors = use_flexible_priors)
  }
  
  # Set up model file based on whether using flexible priors
  if (use_flexible_priors) {
    model_file <- system.file("jags", paste0("model_", model_type, "_flexible.txt"), 
                             package = "mfeargen")
    if (verbose) message("Using flexible prior model")
  } else {
    model_file <- system.file("jags", paste0("model_", model_type, ".txt"), 
                             package = "mfeargen")
  }
  
  if (!file.exists(model_file)) {
    stop("JAGS model file not found: ", model_file)
  }
  
  # Initialize parameters
  if (verbose) message("Initializing MCMC parameters...")
  inits <- initialize_mcmc_parameters(jags_data, n.chains)
  
  # Parameters to monitor
  parameters <- c("alpha", "lambda", "gp", "pi", "w0", "w1", "sigma")
  
  # Run MCMC
  if (verbose) message("Running MCMC sampling...")
  
  fit <- jags(
    data = jags_data,
    inits = inits,
    parameters.to.save = parameters,
    model.file = model_file,
    n.chains = n.chains,
    n.iter = n.iter,
    n.burnin = n.burnin,
    n.thin = n.thin,
    parallel = parallel,
    DIC = TRUE,
    ...
  )
  
  # Post-processing
  if (verbose) message("Post-processing results...")
  
  # Extract group membership
  group_membership <- extract_group_membership(fit)
  
  # Check convergence
  convergence <- check_convergence(fit, verbose = verbose)
  
  # Create mfeargen object
  mfeargen_fit <- list(
    jags_fit = fit,
    model_type = model_type,
    data = data,
    jags_data = jags_data,
    group_membership = group_membership,
    lambda_threshold = jags_data$lambda_threshold,
    convergence = convergence,
    n.iter = n.iter,
    n.burnin = n.burnin,
    n.chains = n.chains,
    call = match.call()
  )
  
  class(mfeargen_fit) <- "mfeargen"
  
  if (verbose) {
    message("Model fitting complete!")
    message("  Model type: ", model_type)
    message("  Participants: ", jags_data$Nparticipants)
    message("  Trials: ", jags_data$Ntrials)
    message("  Convergence: ", ifelse(convergence$converged, "Yes", "No"))
  }
  
  return(mfeargen_fit)
}


#' Initialize MCMC Parameters
#'
#' @description
#' Creates initial values for MCMC chains to improve convergence.
#'
#' @param jags_data List. Prepared JAGS data
#' @param n.chains Integer. Number of chains
#' @return List of initial values for each chain
#' @keywords internal
initialize_mcmc_parameters <- function(jags_data, n.chains) {
  
  n_participants <- jags_data$Nparticipants
  
  inits <- list()
  
  for (chain in 1:n.chains) {
    inits[[chain]] <- list(
      alpha_1 = runif(n_participants, 0.1, 0.9),
      lambda_1 = runif(n_participants, 
                      jags_data$lambda_threshold, 
                      jags_data$lambda_threshold * 3),
      lambda_2 = runif(n_participants, 0, jags_data$lambda_threshold),
      w0 = rnorm(n_participants, 0, 1),
      w1_1 = rgamma(n_participants, 2, 1) + 1,
      gp = sample(1:4, n_participants, replace = TRUE),
      pi = rep(0.25, 4),
      alpha_mu = runif(1, 0.3, 0.7),
      alpha_kappa = runif(1, 2, 8),
      lambda_mu = runif(1, 0.05, 0.15),
      lambda_sigma = runif(1, 0.01, 0.1),
      w0_mu = rnorm(1, 0, 0.5),
      w0_sigma = rgamma(1, 2, 1),
      w1_a = rgamma(1, 2, 1),
      w1_b = rgamma(1, 2, 1),
      sigma = c(runif(1, 0.5, 1.5), runif(1, 1.5, 3))
    )
  }
  
  return(inits)
}


#' Extract Group Membership
#'
#' @description
#' Extracts the most probable group membership for each participant.
#'
#' @param fit JAGS output object
#' @return Data frame with group assignments and probabilities
#' @keywords internal
extract_group_membership <- function(fit) {
  
  if ("gp" %in% names(fit$mean)) {
    # Get posterior mean group assignments
    group_means <- fit$mean$gp
    n_participants <- length(group_means)
    
    # Round to get most probable group
    group_assignment <- round(group_means)
    
    # Get posterior samples for uncertainty
    if ("gp" %in% names(fit$sims.list)) {
      gp_samples <- fit$sims.list$gp
      
      # Calculate probability of each group for each participant
      group_probs <- matrix(0, n_participants, 4)
      for (i in 1:n_participants) {
        for (g in 1:4) {
          group_probs[i, g] <- mean(gp_samples[, i] == g)
        }
      }
      
      colnames(group_probs) <- paste0("prob_group_", 1:4)
      
    } else {
      group_probs <- NULL
    }
    
    # Create output data frame
    membership <- data.frame(
      participant = 1:n_participants,
      group = group_assignment,
      group_label = c("Non-Learner", "Overgeneralizer", 
                     "Physical", "Perceptual")[group_assignment],
      assigned_group = c("Non-Learner", "Overgeneralizer", 
                        "Physical", "Perceptual")[group_assignment]
    )
    
    if (!is.null(group_probs)) {
      membership <- cbind(membership, group_probs)
    }
    
    return(membership)
    
  } else {
    warning("Group membership (gp) not found in model output")
    return(NULL)
  }
}


#' Check MCMC Convergence
#'
#' @description
#' Checks convergence of MCMC chains using Rhat statistics.
#'
#' @param fit JAGS output object
#' @param threshold Numeric. Rhat threshold for convergence (default: 1.1)
#' @param verbose Logical. Print diagnostic messages
#' @return List with convergence information
#' @export
check_convergence <- function(fit, threshold = 1.1, verbose = TRUE) {
  
  # If fit is a numeric vector (Rhat values directly), use the helper function
  if (is.numeric(fit)) {
    result <- check_mcmc_convergence(fit, threshold)
    return(result)
  }
  
  # Handle mfeargen objects
  if (inherits(fit, "mfeargen")) {
    fit <- fit$jags_fit
  }
  
  # Extract Rhat values from JAGS fit object
  rhat_values <- NULL
  
  # Try different locations for Rhat values
  if ("Rhat" %in% names(fit)) {
    # jagsUI output format - Rhat is a list
    rhat_values <- unlist(fit$Rhat)
  } else if (!is.null(fit$BUGSoutput) && !is.null(fit$BUGSoutput$summary)) {
    # Standard JAGS/BUGS output location
    if ("Rhat" %in% colnames(fit$BUGSoutput$summary)) {
      rhat_values <- fit$BUGSoutput$summary[, "Rhat"]
    }
  } else if (!is.null(fit$summary) && "Rhat" %in% colnames(fit$summary)) {
    # In summary
    rhat_values <- fit$summary[, "Rhat"]
  }
  
  # Process Rhat values
  if (!is.null(rhat_values)) {
    rhat_values <- rhat_values[!is.na(rhat_values)]
    
    if (length(rhat_values) == 0) {
      if (verbose) warning("No valid Rhat values found")
      return(list(converged = NA, max_rhat = NA))
    }
    
    max_rhat <- max(rhat_values, na.rm = TRUE)
    converged <- all(rhat_values < threshold, na.rm = TRUE)
    
    if (verbose && !converged) {
      n_poor <- sum(rhat_values >= threshold, na.rm = TRUE)
      message("Warning: ", n_poor, " parameters have Rhat >= ", threshold)
      message("         Maximum Rhat: ", round(max_rhat, 3))
      message("         Consider increasing n.iter or n.burnin")
    }
    
    return(list(
      converged = converged,
      max_rhat = max_rhat,
      n_parameters = length(rhat_values),
      n_poor_convergence = sum(rhat_values >= threshold, na.rm = TRUE)
    ))
    
  } else {
    if (verbose) warning("Rhat values not found in model output")
    return(list(converged = NA, max_rhat = NA))
  }
}


#' Extract Model Parameters
#'
#' @description
#' Extracts posterior distributions of model parameters.
#'
#' @param fit mfeargen object or JAGS output
#' @param parameter Character. Which parameter to extract
#' @param summary_only Logical. Return summary statistics only
#' @return Data frame or matrix of parameter values
#' @export
#' @examples
#' \dontrun{
#' fit <- fit_mfeargen(data)
#' alpha_params <- extract_parameters(fit, "alpha")
#' }
extract_parameters <- function(fit, parameter = "all", summary_only = TRUE) {
  
  # Handle mfeargen object
  if (inherits(fit, "mfeargen")) {
    jags_output <- fit$jags_fit
  } else {
    jags_output <- fit
  }
  
  # Available parameters - check what exists in the output
  if (!is.null(jags_output$mean)) {
    available_params <- names(jags_output$mean)
  } else if (!is.null(jags_output$summary)) {
    # Extract parameter names from summary
    param_names <- rownames(jags_output$summary)
    available_params <- unique(gsub("\\[.*\\]", "", param_names))
  } else {
    stop("Cannot extract parameters: model output format not recognized")
  }
  
  if (parameter == "all") {
    parameters <- available_params
  } else {
    if (!(parameter %in% available_params)) {
      stop("Parameter '", parameter, "' not found in model output")
    }
    parameters <- parameter
  }
  
  # Extract parameters
  results <- list()
  
  for (param in parameters) {
    if (summary_only) {
      # Return summary statistics
      results[[param]] <- data.frame(
        mean = as.vector(jags_output$mean[[param]]),
        sd = as.vector(jags_output$sd[[param]]),
        q2.5 = as.vector(jags_output$q2.5[[param]]),
        q50 = as.vector(jags_output$q50[[param]]),
        q97.5 = as.vector(jags_output$q97.5[[param]])
      )
    } else {
      # Return full posterior samples
      if (param %in% names(jags_output$sims.list)) {
        results[[param]] <- jags_output$sims.list[[param]]
      }
    }
  }
  
  # Add group membership if available and requested
  if (parameter == "all" && inherits(fit, "mfeargen") && !is.null(fit$group_membership)) {
    results$group_membership <- fit$group_membership
  }
  
  if (length(results) == 1) {
    return(results[[1]])
  } else {
    # For multiple parameters, add summary components for backwards compatibility
    results$mean <- jags_output$mean
    results$sd <- jags_output$sd
    results$quantiles <- list(
      q2.5 = jags_output$q2.5,
      q50 = jags_output$q50,
      q97.5 = jags_output$q97.5
    )
    return(results)
  }
}
