#' Compare models fitted to the same observations
#' @param ... Two or more fitted mfeargen models.
#' @param criterion Either "DIC" or "WAIC". WAIC uses individual observed
#'   ratings as predictive units, conditional on participant effects.
#' @return A comparison data frame, ordered by the selected criterion.
#' @export
compare_models <- function(..., criterion = c("DIC", "WAIC")) {
  criterion <- match.arg(criterion)
  models <- list(...)
  if (length(models) < 2L || !all(vapply(models, inherits, logical(1), "mfeargen"))) {
    stop("Supply at least two fitted mfeargen models.")
  }
  reference <- models[[1]]
  same_observations <- vapply(models, function(x) {
    identical(x$jags_data$y, reference$jags_data$y) &&
      identical(attr(x$jags_data, "participants"), attr(reference$jags_data, "participants")) &&
      identical(attr(x$jags_data, "row_index"), attr(reference$jags_data, "row_index")) &&
      identical(x$rating_scale, reference$rating_scale)
  }, logical(1))
  if (!all(same_observations)) stop("Models must be fitted to the same observations and rating scale.")
  model_names <- names(models)
  if (is.null(model_names)) model_names <- rep("", length(models))
  expressions <- vapply(as.list(substitute(list(...)))[-1], deparse1, character(1))
  model_names[!nzchar(model_names)] <- expressions[!nzchar(model_names)]
  values <- vapply(models, function(x) {
    if (criterion == "WAIC") calculate_waic(x)$waic else {
      if (is.null(x$jags_fit$DIC)) NA_real_ else x$jags_fit$DIC
    }
  }, numeric(1))
  if (any(!is.finite(values))) stop(criterion, " is unavailable or non-finite for one or more models.")
  result <- data.frame(
    Model = make.unique(model_names),
    Type = vapply(models, `[[`, character(1), "model_type"),
    N_Participants = vapply(models, function(x) x$jags_data$Nparticipants, numeric(1)),
    N_Trials = vapply(models, function(x) sum(!is.na(x$jags_data$y)), numeric(1)),
    Converged = vapply(models, function(x) isTRUE(x$convergence$converged), logical(1))
  )
  result[[criterion]] <- values
  delta <- values - min(values)
  result[[paste0("Delta_", criterion)]] <- delta
  result[[paste0(criterion, "_Weight")]] <- exp(-delta / 2) / sum(exp(-delta / 2))
  result$Rank <- rank(values, ties.method = "min")
  result <- result[order(result$Rank), , drop = FALSE]
  rownames(result) <- NULL
  structure(result, class = c("mfeargen_comparison", "data.frame"),
            criterion = criterion, best_model = result$Model[1], n_models = length(models))
}

#' Print model comparisons
#' @param x A mfeargen_comparison object.
#' @param ... Arguments passed to print.data.frame.
#' @return The input, invisibly.
#' @export
print.mfeargen_comparison <- function(x, ...) {
  cat("Model comparison using", attr(x, "criterion"), "\n")
  print.data.frame(x, ...)
  invisible(x)
}

log_mean_exp <- function(x) {
  m <- max(x)
  if (is.infinite(m)) return(m)
  m + log(mean(exp(x - m)))
}

#' Calculate conditional trial-level WAIC
#' @param fit_object A fitted mfeargen model with monitored loglik draws.
#' @return WAIC, its standard error, effective parameter count, and pointwise
#'   contributions. Missing and padded observations are excluded.
#' @details Each observed rating is a predictive unit. This estimates prediction
#'   for fitted participants, conditional on their latent effects; it does not
#'   estimate prediction for new participants. Likelihoods are on the normalized
#'   0-10 rating scale and must only be compared across identical observations.
#' @export
calculate_waic <- function(fit_object) {
  if (!inherits(fit_object, "mfeargen")) stop("fit_object must be a fitted mfeargen model.")
  samples <- fit_object$jags_fit$sims.list$loglik
  if (is.null(samples)) stop("Log likelihood is unavailable; refit with loglik monitoring.")
  if (length(dim(samples)) != 3L) stop("Expected loglik dimensions: draws by participants by trials.")
  observed <- !is.na(fit_object$jags_data$y)
  ll <- matrix(samples, nrow = dim(samples)[1])[, as.vector(observed), drop = FALSE]
  if (nrow(ll) < 2L || !ncol(ll) || any(!is.finite(ll))) {
    stop("WAIC requires at least two finite draws for every observed rating.")
  }
  lppd <- apply(ll, 2, log_mean_exp)
  penalty <- apply(ll, 2, stats::var)
  pointwise <- -2 * (lppd - penalty)
  structure(list(waic = sum(pointwise), lppd = sum(lppd), p_waic = sum(penalty),
                 se_waic = if (length(pointwise) > 1L) sqrt(length(pointwise) * stats::var(pointwise)) else NA_real_,
                 pointwise = pointwise, unit = "observed rating, conditional on participant effects"),
            class = "mfeargen_waic")
}

#' Cross-validate held-out ratings within participants
#' @param data A trial-level data frame accepted by fit_mfeargen.
#' @param model_type Either "1cs" or "2cs".
#' @param k_folds Number of rating folds.
#' @param seed Optional integer seed for fold assignment and MCMC.
#' @param ... Arguments passed to fit_mfeargen, including rating_scale.
#' @return Fold fits, held-out log predictive densities and mean squared errors,
#'   fold assignments, and aggregate summaries.
#' @details Ratings are randomly held out within each participant; all stimulus
#'   and reinforcement histories remain available. Each fold refits the model
#'   with its held-out ratings set to NA. This evaluates missing-rating prediction
#'   for existing participants, not forecasts of future trials or new participants.
#' @export
cross_validate_mfeargen <- function(data, model_type = c("1cs", "2cs"),
                                   k_folds = 10, seed = NULL, ...) {
  model_type <- match.arg(model_type)
  if (!is.data.frame(data) || !all(c("participant", "US_expect") %in% names(data))) {
    stop("data must contain participant and US_expect columns.")
  }
  if (!is.numeric(k_folds) || length(k_folds) != 1L || !is.finite(k_folds) || k_folds < 2L || k_folds != as.integer(k_folds)) {
    stop("k_folds must be an integer of at least two.")
  }
  if (!is.null(seed)) {
    validate_scalar(seed, "seed", lower = 0, upper = .Machine$integer.max)
    if (seed != floor(seed) || seed + k_folds > .Machine$integer.max) {
      stop("seed must be an integer allowing a distinct seed for each fold")
    }
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
  observed <- which(!is.na(data$US_expect))
  counts <- table(as.character(data$participant[observed]))
  if (!length(counts) || any(counts < k_folds)) {
    stop("Each participant must have at least k_folds observed ratings.")
  }
  assignments <- rep(NA_integer_, nrow(data))
  for (id in unique(data$participant)) {
    rows <- observed[data$participant[observed] == id]
    assignments[rows] <- sample(rep(seq_len(k_folds), length.out = length(rows)))
  }
  args <- list(...)
  if (!is.null(args$prepared_data)) {
    stop("prepared_data cannot be reused across folds; provide trial-level data only")
  }
  if (is.null(args$rating_scale)) {
    args$rating_scale <- attr(data, "rating_scale")
    if (is.null(args$rating_scale)) args$rating_scale <- detect_rating_scale(data$US_expect)
  }
  fits <- vector("list", k_folds)
  metrics <- vector("list", k_folds)
  for (k in seq_len(k_folds)) {
    rows <- which(assignments == k)
    training <- data
    training$US_expect[rows] <- NA_real_
    fold_seed <- if (is.null(seed)) NULL else seed + k
    fit <- do.call(fit_mfeargen, c(list(data = training, model_type = model_type, seed = fold_seed), args))
    map <- attr(fit$jags_data, "row_index")
    cells <- match(rows, as.vector(map))
    if (anyNA(cells)) stop("Held-out rows could not be matched to prepared trials.")
    theta_draws <- fit$jags_fit$sims.list$theta
    theta <- matrix(theta_draws, nrow = dim(theta_draws)[1])[, cells, drop = FALSE]
    gp <- matrix(fit$jags_fit$sims.list$gp, ncol = fit$jags_data$Nparticipants)
    sigma <- matrix(fit$jags_fit$sims.list$sigma, ncol = 2L)
    ids <- ((cells - 1L) %% fit$jags_data$Nparticipants) + 1L
    bounds <- args$rating_scale
    y <- (data$US_expect[rows] - bounds[1]) / diff(bounds) * 10
    log_density <- vapply(seq_along(rows), function(j) {
      sd <- sigma[cbind(seq_len(nrow(theta)), ifelse(gp[, ids[j]] == 1L, 2L, 1L))]
      log_mean_exp(stats::dnorm(y[j], theta[, j], sd, log = TRUE)) + log(10 / diff(bounds))
    }, numeric(1))
    predictions <- colMeans(theta) / 10 * diff(bounds) + bounds[1]
    metrics[[k]] <- data.frame(fold = k, n_test = length(rows),
                              loglik_test = sum(log_density),
                              mse_test = mean((data$US_expect[rows] - predictions)^2),
                              converged = isTRUE(fit$convergence$converged))
    fits[[k]] <- list(fold = k, fit = fit, test_rows = rows,
                      predictions = predictions, log_predictive_density = log_density)
  }
  fold_metrics <- do.call(rbind, metrics)
  structure(list(fold_results = fits, fold_metrics = fold_metrics,
                 fold_assignment = assignments,
                 summary = list(mean_loglik = sum(fold_metrics$loglik_test) / length(observed),
                                mean_mse = stats::weighted.mean(fold_metrics$mse_test, fold_metrics$n_test),
                                convergence_rate = mean(fold_metrics$converged),
                                successful_folds = k_folds),
                 k_folds = k_folds, model_type = model_type), class = "mfeargen_cv")
}

#' Describe the model comparison result
#' @param comparison_result A mfeargen_comparison object.
#' @return A brief description of the criterion difference.
#' @export
model_selection_advice <- function(comparison_result) {
  if (!inherits(comparison_result, "mfeargen_comparison")) stop("Expected a mfeargen_comparison object.")
  criterion <- attr(comparison_result, "criterion")
  delta <- comparison_result[[paste0("Delta_", criterion)]][2]
  paste0(comparison_result$Model[1], " has the lowest ", criterion,
         "; the next model differs by ", round(delta, 2),
         ". Review convergence and predictive uncertainty before selecting a model.")
}
