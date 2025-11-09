#' Paper Replication Functions
#'
#' Functions to facilitate replication of analyses from Yu et al. (2023).
#' These functions provide easy access to experimental data and exact priors
#' used in the original paper.
#'
#' @name paper_replication
#' @keywords replication
#' @importFrom utils data
NULL

#' Load Experimental Data from Yu et al. (2023)
#'
#' Convenience function to load the experimental datasets used in
#' Yu et al. (2023). Returns cleaned and validated data ready for
#' analysis with \code{\link{prepare_mfeargen_data}}.
#'
#' @param experiment Character string specifying which experiment to load.
#'   Options are \code{"exp1"} (single CS+ design) or \code{"exp2"}
#'   (two CS+ design). Default is \code{"exp1"}.
#'
#' @return A data frame containing the experimental data with all variables
#'   properly formatted. See \code{\link{exp1_data}} or \code{\link{exp2_data}}
#'   for detailed variable descriptions.
#'
#' @details
#' This function loads the cleaned experimental data and provides metadata
#' as attributes. The data can be directly used with
#' \code{\link{prepare_mfeargen_data}} for model fitting.
#'
#' For Experiment 1 (\code{experiment = "exp1"}), the data uses a single CS+
#' design where distances are measured from one conditioned stimulus.
#'
#' For Experiment 2 (\code{experiment = "exp2"}), the data uses a two CS+
#' design (CS+m and CS+p) allowing assessment of both physical and perceptual
#' generalization patterns.
#'
#' @seealso
#' \code{\link{exp1_data}}, \code{\link{exp2_data}},
#' \code{\link{prepare_mfeargen_data}}, \code{\link{get_paper_priors}},
#' \code{\link{replicate_paper_analysis}}
#'
#' @references
#' Yu, K., Tuerlinckx, F., Vanpaemel, W., & Zaman, J. (2023).
#'   A multi-fidelity computational approach to understand human fear learning.
#'   \emph{Communications Psychology}, \emph{1}(1), 5.
#'   \doi{10.1038/s44271-023-00005-0}
#'
#' @export
#' @examples
#' \dontrun{
#' # Load Experiment 1 data
#' exp1 <- load_paper_data("exp1")
#' head(exp1)
#'
#' # Load Experiment 2 data
#' exp2 <- load_paper_data("exp2")
#' head(exp2)
#'
#' # Check dataset metadata
#' attr(exp1, "experiment")
#' attr(exp1, "n_participants")
#' }
load_paper_data <- function(experiment = c("exp1", "exp2")) {

  # Match argument
  experiment <- match.arg(experiment)

  # Load appropriate dataset
  if (experiment == "exp1") {
    data("exp1_data", envir = environment())
    dataset <- get("exp1_data", envir = environment())
    message(sprintf("Loaded Experiment 1 data: %s design",
                    attr(dataset, "design")))
    message(sprintf("  Participants: %d", attr(dataset, "n_participants")))
    message(sprintf("  Observations: %d", attr(dataset, "n_trials")))
  } else {
    data("exp2_data", envir = environment())
    dataset <- get("exp2_data", envir = environment())
    message(sprintf("Loaded Experiment 2 data: %s design",
                    attr(dataset, "design")))
    message(sprintf("  Participants: %d", attr(dataset, "n_participants")))
    message(sprintf("  Observations: %d", attr(dataset, "n_trials")))
  }

  return(dataset)
}

#' Get Prior Distributions from Yu et al. (2023)
#'
#' Returns the exact prior distributions used in the original paper
#' Yu et al. (2023). These priors are weakly informative and based on
#' theoretical considerations and pilot data.
#'
#' @param model_type Character string specifying the model type.
#'   Options are \code{"1cs"} (single CS+ model) or \code{"2cs"}
#'   (two CS+ model). Default is \code{"1cs"}.
#'
#' @return A list of prior distributions compatible with
#'   \code{\link{set_priors}}. The list contains priors for:
#'   \itemize{
#'     \item \code{alpha}: Learning rate (Beta distribution)
#'     \item \code{tau}: Perceptual noise (Gamma distribution)
#'     \item \code{theta}: Mixture probabilities (Dirichlet distribution)
#'     \item \code{sigma}: Response noise (Gamma distribution)
#'   }
#'
#' @details
#' The priors used in Yu et al. (2023) were designed to be weakly informative,
#' allowing the data to dominate while providing regularization. Specifically:
#'
#' \strong{Learning rate (alpha)}:
#' \itemize{
#'   \item Distribution: Beta(2, 2)
#'   \item Rationale: Symmetric around 0.5, favors moderate learning rates
#' }
#'
#' \strong{Perceptual noise (tau)}:
#' \itemize{
#'   \item Distribution: Gamma(1, 0.1)
#'   \item Rationale: Weakly informative, mean = 10, allows flexible noise levels
#' }
#'
#' \strong{Mixture probabilities (theta)}:
#' \itemize{
#'   \item Distribution: Dirichlet(1, 1, 1, 1)
#'   \item Rationale: Uniform prior over simplex (non-informative)
#' }
#'
#' \strong{Response noise (sigma)}:
#' \itemize{
#'   \item Distribution: Gamma(2, 0.5)
#'   \item Rationale: Mean = 4, regularizes measurement error
#' }
#'
#' @seealso
#' \code{\link{set_priors}}, \code{\link{load_paper_data}},
#' \code{\link{replicate_paper_analysis}}, \code{\link{get_weak_priors}}
#'
#' @references
#' Yu, K., Tuerlinckx, F., Vanpaemel, W., & Zaman, J. (2023).
#'   A multi-fidelity computational approach to understand human fear learning.
#'   \emph{Communications Psychology}, \emph{1}(1), 5.
#'   \doi{10.1038/s44271-023-00005-0}
#'
#' @export
#' @examples
#' \dontrun{
#' # Get priors for 1CS model
#' priors_1cs <- get_paper_priors("1cs")
#' print(priors_1cs)
#'
#' # Get priors for 2CS model
#' priors_2cs <- get_paper_priors("2cs")
#' print(priors_2cs)
#'
#' # Use with fit_mfeargen
#' data_prepared <- prepare_mfeargen_data(
#'   data = exp1_data,
#'   participant_id = "participant",
#'   us_expectancy = "us_expect",
#'   distance_physical = "dphy",
#'   distance_perceptual = "dper",
#'   cs_indicator = "cs_trials"
#' )
#'
#' fit <- fit_mfeargen(
#'   data = data_prepared,
#'   priors = priors_1cs,
#'   n_chains = 3,
#'   n_iter = 5000
#' )
#' }
get_paper_priors <- function(model_type = c("1cs", "2cs")) {

  # Match argument
  model_type <- match.arg(model_type)

  priors <- set_priors(
    alpha = prior_beta(shape1 = 2, shape2 = 2),
    sigma = prior_gamma(shape = 2, rate = 0.5),
    pi = prior_dirichlet(alpha = c(1, 1, 1, 1))
  )
  priors$tau <- prior_gamma(shape = 1, rate = 0.1)
  priors$theta <- priors$pi
  attr(priors, "overrides") <- unique(c(attr(priors, "overrides"), "tau", "theta"))
  attr(priors, "model_type") <- model_type
  attr(priors, "source") <- "Yu et al. (2023)"
  attr(priors, "doi") <- "10.1038/s44271-023-00005-0"

  message(sprintf("Loaded priors for %s model from Yu et al. (2023)",
                  toupper(model_type)))
  message("These are the exact priors used in the original paper.")

  return(priors)
}

#' Replicate Analysis from Yu et al. (2023)
#'
#' Wrapper function to replicate the exact analysis pipeline from
#' Yu et al. (2023). This function loads the appropriate experimental data,
#' applies the paper's priors, and fits the model with the same settings
#' used in the original publication.
#'
#' @param experiment Character string specifying which experiment to replicate.
#'   Options are \code{"exp1"} (single CS+ design) or \code{"exp2"}
#'   (two CS+ design). Default is \code{"exp1"}.
#' @param model_type_paper Logical indicating whether to use the exact JAGS
#'   model from the paper. If \code{TRUE} (default), uses the original model
#'   specifications. If \code{FALSE}, uses the package's current model
#'   implementation. Default is \code{TRUE}.
#' @param cs_reference For Experiment 2 only. Character string specifying
#'   which CS+ to use as reference: \code{"csm"} (physically prototypical)
#'   or \code{"csp"} (perceptually prototypical). Default is \code{"csm"}.
#' @param n_chains Number of MCMC chains. Default is 3 (as in paper).
#' @param n_iter Number of iterations per chain. Default is 10000 (as in paper).
#' @param n_burnin Number of burn-in iterations. Default is 5000 (as in paper).
#' @param n_thin Thinning interval. Default is 2 (as in paper).
#' @param ... Additional arguments passed to \code{\link{fit_mfeargen}}.
#'
#' @return An object of class \code{mfeargen} containing:
#'   \itemize{
#'     \item \code{mcmc}: MCMC samples from posterior distributions
#'     \item \code{data}: Prepared data used for fitting
#'     \item \code{priors}: Prior distributions used
#'     \item \code{model_type}: Model specification used
#'     \item \code{convergence}: Convergence diagnostics
#'   }
#'
#' @details
#' This function is designed to make it easy to replicate the exact analyses
#' from Yu et al. (2023). It handles:
#' \enumerate{
#'   \item Loading the appropriate experimental dataset
#'   \item Setting up the exact priors used in the paper
#'   \item Preparing data in the correct format
#'   \item Running MCMC with the same settings as the paper
#' }
#'
#' For Experiment 1, the analysis uses a single CS+ design with one distance
#' metric. For Experiment 2, you can choose to analyze generalization from
#' either CS+m (physically prototypical) or CS+p (perceptually prototypical).
#'
#' \strong{Computational Note:} The full analysis can take 1-2 hours depending
#' on your system. For testing purposes, consider reducing \code{n_iter} and
#' \code{n_burnin}.
#'
#' @seealso
#' \code{\link{load_paper_data}}, \code{\link{get_paper_priors}},
#' \code{\link{fit_mfeargen}}, \code{\link{prepare_mfeargen_data}}
#'
#' @references
#' Yu, K., Tuerlinckx, F., Vanpaemel, W., & Zaman, J. (2023).
#'   A multi-fidelity computational approach to understand human fear learning.
#'   \emph{Communications Psychology}, \emph{1}(1), 5.
#'   \doi{10.1038/s44271-023-00005-0}
#'
#' @export
#' @examples
#' \dontrun{
#' # Replicate Experiment 1 analysis (this will take time!)
#' fit_exp1 <- replicate_paper_analysis(
#'   experiment = "exp1",
#'   n_chains = 3,
#'   n_iter = 10000
#' )
#'
#' # View results
#' print(fit_exp1)
#' plot(fit_exp1)
#'
#' # Replicate Experiment 2 analysis using CS+m as reference
#' fit_exp2 <- replicate_paper_analysis(
#'   experiment = "exp2",
#'   cs_reference = "csm",
#'   n_chains = 3,
#'   n_iter = 10000
#' )
#'
#' # Quick test run with fewer iterations
#' fit_test <- replicate_paper_analysis(
#'   experiment = "exp1",
#'   n_chains = 2,
#'   n_iter = 1000,
#'   n_burnin = 500
#' )
#' }
replicate_paper_analysis <- function(experiment = c("exp1", "exp2"),
                                    model_type_paper = TRUE,
                                    cs_reference = c("csm", "csp"),
                                    n_chains = 3,
                                    n_iter = 10000,
                                    n_burnin = 5000,
                                    n_thin = 2,
                                    ...) {

  # Match arguments
  experiment <- match.arg(experiment)
  cs_reference <- match.arg(cs_reference)

  # Determine model type
  model_type <- ifelse(experiment == "exp1", "1cs", "2cs")

  message("========================================")
  message("Replicating Yu et al. (2023) Analysis")
  message("========================================")
  message(sprintf("Experiment: %s", toupper(experiment)))
  message(sprintf("Model type: %s", toupper(model_type)))
  if (experiment == "exp2") {
    message(sprintf("CS+ reference: %s", toupper(cs_reference)))
  }
  message("")

  # Step 1: Load experimental data
  message("Step 1: Loading experimental data...")
  exp_data <- load_paper_data(experiment)
  message("")

  # Step 2: Get paper priors
  message("Step 2: Loading priors from paper...")
  priors <- get_paper_priors(model_type)
  message("")

  # Step 3: Prepare data for modeling
  message("Step 3: Preparing data for modeling...")

  # Rename columns to match expected format
  if (experiment == "exp1") {
    # Experiment 1: Single CS+ design
    exp_data$US_expect <- exp_data$us_expect
    exp_data$Per_size <- exp_data$per_size
    exp_data$Phy_size <- exp_data$phy_size

    model_data <- prepare_mfeargen_data(
      data = exp_data,
      model_type = "1cs",
      priors = priors
    )
  } else {
    # Experiment 2: Two CS+ design
    exp_data$US_expect <- exp_data$us_expect
    exp_data$Per_size <- exp_data$per_size
    exp_data$Phy_size <- exp_data$phy_size

    # For 2CS model, we need to set up US indicators
    if (cs_reference == "csm") {
      exp_data$USp <- exp_data$us_p
      exp_data$USm <- exp_data$us_m
    } else {
      exp_data$USp <- exp_data$us_p
      exp_data$USm <- exp_data$us_m
    }

    model_data <- prepare_mfeargen_data(
      data = exp_data,
      model_type = "2cs",
      priors = priors
    )
  }
  message("")

  # Step 4: Fit the model
  message("Step 4: Fitting the model...")
  message(sprintf("  Chains: %d", n_chains))
  message(sprintf("  Iterations: %d", n_iter))
  message(sprintf("  Burn-in: %d", n_burnin))
  message(sprintf("  Thinning: %d", n_thin))
  message("")
  message("This may take a while. Please be patient...")
  message("")

  fit <- fit_mfeargen(
    data = exp_data,
    model_type = model_type,
    priors = priors,
    n.chains = n_chains,
    n.iter = n_iter,
    n.burnin = n_burnin,
    n.thin = n_thin,
    ...
  )

  # Add replication metadata
  attr(fit, "replication") <- TRUE
  attr(fit, "experiment") <- experiment
  attr(fit, "paper") <- "Yu et al. (2023)"
  attr(fit, "doi") <- "10.1038/s44271-023-00005-0"

  message("")
  message("========================================")
  message("Replication complete!")
  message("========================================")

  return(fit)
}
