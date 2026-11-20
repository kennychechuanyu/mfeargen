#' Load the bundled experimental data
#' @param experiment Either "exp1" or "exp2".
#' @return The documented data frame, retaining its source metadata.
#' @seealso exp1_data, exp2_data
#' @export
load_paper_data <- function(experiment = c("exp1", "exp2")) {
  experiment <- match.arg(experiment)
  environment <- new.env(parent = emptyenv())
  utils::data(list = paste0(experiment, "_data"), package = "mfeargen", envir = environment)
  environment[[paste0(experiment, "_data")]]
}

#' Hyperpriors corresponding to Table 1 of Yu et al. (2023)
#' @param model_type Either "1cs" or "2cs".
#' @return A mfeargen_priors object. These match the published hyperprior
#'   families; they do not make the package likelihood an exact reproduction
#'   of the archived analysis. Stimulus units and rating bounds still matter.
#' @references Yu et al. (2023), Communications Psychology 1, 5.
#'   \doi{10.1038/s44271-023-00005-0}, Table 1.
#' @export
get_paper_priors <- function(model_type = c("1cs", "2cs")) {
  model_type <- match.arg(model_type)
  priors <- set_priors_flexible(
    lambda_mu = prior_normal(mean = 0.1, sd = 1, truncate = c(0, Inf)),
    lambda_sigma = prior_uniform(1e-9, 1),
    alpha_mu = prior_beta(1, 1),
    alpha_kappa = prior_uniform(1, 10),
    sigma1 = prior_uniform(1e-9, 1.5),
    sigma2 = prior_uniform(1.5, 3),
    pi = prior_dirichlet(rep(1, 4)),
    w0_mu = prior_normal(0, 10),
    w0_sigma = prior_half_cauchy(scale = 2),
    w1_a = prior_half_cauchy(scale = 2),
    w1_b = prior_half_cauchy(scale = 2)
  )
  attr(priors, "model_type") <- model_type
  attr(priors, "source") <- "Yu et al. (2023), Table 1"
  priors
}

prepare_paper_trials <- function(experiment) {
  data <- load_paper_data(experiment)
  # ITI rows have no presented stimulus and cannot define a learning trial.
  keep <- !is.na(data$phy_size)
  if (experiment == "exp1") {
    keep <- keep & !is.na(data$dphy) & !is.na(data$dper)
  } else {
    distances <- c("d_phy_p", "d_phy_m", "d_per_p", "d_per_m")
    keep <- keep & stats::complete.cases(data[distances])
  }
  data <- data[keep, , drop = FALSE]
  data$US_expect <- data$us_expect
  data$Phy_size <- data$phy_size
  data$Per_size <- data$per_size
  if (experiment == "exp1") {
    data$k <- data$cs_trials
    data$d_phy <- data$dphy
    data$d_per_old <- data$dper
  } else {
    data$k_plus <- data$csp_trials
    data$k_minus <- data$csm_trials
    data$USp <- ifelse(data$k_plus == 1, data$us_p, 0)
    data$USm <- ifelse(data$k_minus == 1, data$us_m, 0)
    data$d_p_phy <- data$d_phy_p
    data$d_m_phy <- data$d_phy_m
    data$d_p_per <- data$d_per_p
    data$d_m_per <- data$d_per_m
  }
  attr(data, "rating_scale") <- c(1, 10)
  data
}

#' Reanalyze the bundled study data
#' @description Fits the package's current model to the bundled trial data.
#'   Exact replication of the archived publication model is not available.
#' @param experiment Either "exp1" or "exp2".
#' @param model_type_paper Must be FALSE. TRUE requests an unsupported exact
#'   replication and raises an error before sampling.
#' @param cs_reference The historical argument; only "csm" is supported.
#'   The 2cs model jointly uses both cue histories, so switching a single
#'   reference is not a separate analysis.
#' @param n_chains Number of chains.
#' @param n_iter Iterations per chain.
#' @param n_burnin Burn-in iterations.
#' @param n_thin Thinning interval.
#' @param ... Additional arguments to fit_mfeargen; priors and rating_scale
#'   may override the published hyperpriors and 1-10 rating bounds.
#' @return A fitted mfeargen model with reanalysis metadata.
#' @export
replicate_paper_analysis <- function(experiment = c("exp1", "exp2"),
                                     model_type_paper = FALSE,
                                     cs_reference = c("csm", "csp"),
                                     n_chains = 3, n_iter = 10000,
                                     n_burnin = 5000, n_thin = 2, ...) {
  experiment <- match.arg(experiment)
  cs_reference <- match.arg(cs_reference)
  if (!identical(model_type_paper, FALSE)) {
    stop("Exact paper replication is unavailable. Use model_type_paper = FALSE for a reanalysis with the current model; archived code is at https://osf.io/sxjak/.")
  }
  if (experiment == "exp2" && cs_reference != "csm") {
    stop("The current 2cs model jointly uses both cue histories; cs_reference = 'csp' is unsupported.")
  }
  model_type <- if (experiment == "exp1") "1cs" else "2cs"
  args <- list(...)
  if (is.null(args$priors)) args$priors <- get_paper_priors(model_type)
  if (is.null(args$rating_scale)) args$rating_scale <- c(1, 10)
  fit <- do.call(fit_mfeargen, c(list(
    data = prepare_paper_trials(experiment), model_type = model_type,
    n.chains = n_chains, n.iter = n_iter, n.burnin = n_burnin, n.thin = n_thin), args))
  attr(fit, "replication") <- FALSE
  attr(fit, "reanalysis") <- TRUE
  attr(fit, "experiment") <- experiment
  attr(fit, "doi") <- "10.1038/s44271-023-00005-0"
  fit
}
