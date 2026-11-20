#' @title Plotting Functions for mfeargen
#' @description Visualization functions for model results
#' @name plotting
#' @importFrom stats aggregate cor setNames
NULL

# Global variables for ggplot2 (to avoid R CMD check notes)
utils::globalVariables(c(
  "..density..", "Count", "Recovered", "True", "US_expect",
  "alpha_recovered", "alpha_true", "count", "group", "group_label",
  "iteration", "label", "lambda_recovered", "lambda_true",
  "participant", "stimulus", "value", "chain", "density", "probability",
  "statistic", "estimate", "observed", "lower", "upper", "predicted"
))

#' Plot mfeargen Results
#'
#' @description
#' Main plotting function for mfeargen objects. Creates diagnostic and result plots.
#'
#' @param x mfeargen object
#' @param type Character. Plot type: "diagnostic", "parameters", "groups", "gradient", "predictive", or "all"
#' @param ... Additional arguments passed to specific plot functions
#' @return ggplot object or list of plots
#' @export
#' @import ggplot2
#' @examples
#' \dontrun{
#' fit <- fit_mfeargen(data)
#' plot(fit)
#' plot(fit, type = "groups")
#' }
plot.mfeargen <- function(x, type = "diagnostic", ...) {
  
  type <- match.arg(type, c("diagnostic", "parameters", "groups", 
                            "gradient", "predictive", "all"))
  
  if (type == "diagnostic") {
    return(plot_mcmc_diagnostics(x, ...))
  } else if (type == "parameters") {
    return(plot_parameters(x, ...))
  } else if (type == "groups") {
    return(plot_group_probabilities(x, ...))
  } else if (type == "gradient") {
    return(plot_generalization_gradient(x, ...))
  } else if (type == "predictive") {
    return(plot_posterior_predictive(x, ...))
  } else if (type == "all") {
    plots <- list(
      diagnostic = plot_mcmc_diagnostics(x, ...),
      parameters = plot_parameters(x, ...),
      groups = plot_group_probabilities(x, ...),
      gradient = plot_generalization_gradient(x, ...),
      predictive = plot_posterior_predictive(x, ...)
    )
    return(plots)
  }
}

#' Plot MCMC Diagnostics
#'
#' @description Creates diagnostic plots for MCMC convergence
#' @param fit mfeargen object
#' @param parameters Character vector of parameters to plot
#' @return ggplot object
#' @export
plot_mcmc_diagnostics <- function(fit, parameters = c("alpha", "lambda")) {
  
  if (!inherits(fit, "mfeargen")) {
    stop("Input must be a mfeargen object")
  }
  
  chains <- fit$jags_fit$samples
  if (is.null(chains)) stop("Chain samples are unavailable for this fit.")
  plot_data <- list()
  for (chain in seq_along(chains)) {
    samples <- as.matrix(chains[[chain]])
    columns <- colnames(samples)
    keep <- columns[vapply(columns, function(nm) {
      any(vapply(parameters, function(p) nm == p || startsWith(nm, paste0(p, "[")), logical(1)))
    }, logical(1))]
    for (nm in utils::head(keep, 10L)) {
      plot_data[[length(plot_data) + 1L]] <- data.frame(
        iteration = as.numeric(stats::time(chains[[chain]])), value = samples[, nm],
        parameter = nm, chain = factor(chain))
    }
  }
  if (!length(plot_data)) stop("None of the requested parameters were monitored.")
  ggplot(do.call(rbind, plot_data), aes(x = iteration, y = value, color = chain)) +
    geom_line(alpha = 0.7) + facet_wrap(~ parameter, scales = "free_y") +
    labs(title = "MCMC traces", x = "Iteration", y = "Value", color = "Chain") + theme_minimal()
}

#' Plot Group Probabilities
#'
#' @description Creates bar plot of group membership probabilities
#' @param fit mfeargen object
#' @return ggplot object
#' @export
plot_group_probabilities <- function(fit) {
  
  if (!inherits(fit, "mfeargen")) {
    stop("Input must be a mfeargen object")
  }
  
  membership <- extract_group_membership(fit$jags_fit,
                                         attr(fit$jags_data, "participants"))
  probabilities <- as.matrix(membership[paste0("prob_group_", 1:4)])
  labels <- get_group_labels(FALSE)
  plot_data <- data.frame(
    participant = factor(rep(as.character(membership$participant), times = 4),
                         levels = as.character(membership$participant)),
    group_label = factor(rep(labels, each = nrow(membership)), levels = labels),
    probability = as.vector(probabilities))
  ggplot(plot_data, aes(x = participant, y = probability, fill = group_label)) +
    geom_col() +
    scale_fill_manual(values = c("#CC79A7", "#F0E442", "#56B4E9", "#D55E00")) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "Posterior group probabilities", x = "Participant", y = "Probability", fill = "Group") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Plot Parameters
#'
#' @description Shows the distribution of posterior participant means.
#' @param fit mfeargen object
#' @param parameter Character. Which parameter to plot
#' @return ggplot object
#' @export
plot_parameters <- function(fit, parameter = "alpha") {
  
  if (!inherits(fit, "mfeargen")) {
    stop("Input must be a mfeargen object")
  }
  
  mcmc_output <- fit$jags_fit
  
  # Try different locations for the parameter
  param_values <- NULL
  
  if (!is.null(mcmc_output$mean) && parameter %in% names(mcmc_output$mean)) {
    param_values <- mcmc_output$mean[[parameter]]
  } else if (!is.null(mcmc_output$summary)) {
    # Extract from summary
    param_rows <- grepl(paste0("^", parameter, "\\["), rownames(mcmc_output$summary))
    if (any(param_rows)) {
      param_values <- mcmc_output$summary[param_rows, "mean"]
    }
  }
  
  if (is.null(param_values)) {
    stop("Parameter '", parameter, "' not found")
  }
  
  plot_data <- data.frame(value = as.vector(param_values))
  p <- ggplot(plot_data, aes(x = value))
  if (length(param_values) > 1L && stats::sd(as.vector(param_values)) > 0) {
    p <- p + geom_histogram(aes(y = after_stat(density)), bins = 30,
                            fill = "lightblue", alpha = 0.7) +
      geom_density(color = "darkblue", linewidth = 1)
  } else {
    p <- p + geom_point(aes(y = 1), size = 3)
  }
  p + labs(title = paste("Posterior means of", parameter), x = parameter,
           y = if (length(param_values) > 1L) "Density" else "Participant") + theme_minimal()
}

#' Plot Observed and Fitted Generalization Gradients
#'
#' @param fit A fitted mfeargen model.
#' @param group_separate Plot separate panels for posterior modal groups.
#' @param level Credible interval width for the mean gradient.
#' @return A ggplot object. Points are observed means; lines and ribbons summarize
#'   the posterior expected mean of the same trial histories. Group panels use
#'   modal assignments for display; consult group probabilities before interpreting
#'   uncertain assignments. With phase labels, only generalization trials are used.
#' @export
plot_generalization_gradient <- function(fit, group_separate = TRUE, level = 0.95) {
  if (!inherits(fit, "mfeargen")) stop("fit must be a fitted mfeargen model.")
  validate_scalar(level, "level", lower = 0, upper = 1, lower_open = TRUE, upper_open = TRUE)
  if (!is.logical(group_separate) || length(group_separate) != 1L || is.na(group_separate)) {
    stop("group_separate must be TRUE or FALSE.")
  }
  data <- fit$data
  data$.input_row <- seq_len(nrow(data))
  if ("phase" %in% names(data)) data <- data[data$phase == "generalization", , drop = FALSE]
  if (!nrow(data)) stop("No generalization trials are available to plot.")
  if ("Phy_size" %in% names(data)) data$stimulus <- data$Phy_size
  if (!"stimulus" %in% names(data) || anyNA(data$stimulus)) {
    stop("Provide stimulus or Phy_size values to plot a generalization gradient.")
  }
  membership <- fit$group_membership
  if (group_separate) {
    if (is.null(membership)) stop("Group membership is unavailable for this fit.")
    data$group_label <- membership$group_label[match(data$participant, membership$participant)]
  } else data$group_label <- "All participants"
  map <- attr(fit$jags_data, "row_index")
  if (is.null(map)) stop("Trial-row mappings are unavailable; refit the model.")
  cells <- match(data$.input_row, as.vector(map))
  if (anyNA(cells) || anyNA(data$group_label)) stop("Trial or participant mappings could not be matched.")
  draws <- posterior_response_draws(fit)
  predictions <- matrix(draws, nrow = dim(draws)[1])[, cells, drop = FALSE]
  grouping <- interaction(data$stimulus, data$group_label, drop = TRUE, lex.order = TRUE)
  tail <- (1 - level) / 2
  rows <- lapply(split(seq_len(nrow(data)), grouping), function(index) {
    predicted <- rowMeans(predictions[, index, drop = FALSE])
    interval <- stats::quantile(predicted, c(tail, 1 - tail))
    observed <- data$US_expect[index]
    data.frame(stimulus = data$stimulus[index[1]], group_label = data$group_label[index[1]],
               observed = if (any(!is.na(observed))) mean(observed, na.rm = TRUE) else NA_real_,
               predicted = mean(predicted), lower = interval[1], upper = interval[2], row.names = NULL)
  })
  plot_data <- do.call(rbind, rows)
  line_data <- plot_data[stats::ave(seq_len(nrow(plot_data)), plot_data$group_label, FUN = length) > 1L, , drop = FALSE]
  p <- ggplot(plot_data, aes(x = stimulus, group = group_label)) +
    geom_ribbon(data = line_data, aes(ymin = lower, ymax = upper), fill = "#56B4E9", alpha = 0.25) +
    geom_line(data = line_data, aes(y = predicted), color = "#0072B2", linewidth = 0.9) +
    geom_point(aes(y = predicted), color = "#0072B2", size = 1.5) +
    geom_point(aes(y = observed), na.rm = TRUE, size = 2) +
    labs(title = "Generalization gradients", x = "Stimulus", y = "US expectancy",
         subtitle = "Points: observed means; line and interval: expected response") + theme_minimal()
  if (group_separate) p <- p + facet_wrap(~ group_label)
  p
}

#' Plot Parameter Recovery
#'
#' @description Creates diagnostic plots for parameter recovery
#' @param recovery mfeargen_recovery object
#' @param parameters Character vector of parameters to plot
#' @return ggplot object
#' @export
plot_parameter_recovery <- function(recovery, 
                                   parameters = c("alpha", "lambda", "group")) {
  
  if (!inherits(recovery, "mfeargen_recovery")) {
    stop("Input must be a mfeargen_recovery object")
  }
  
  data <- recovery$metrics$combined_data
  plots <- list()
  
  # Alpha recovery
  if ("alpha" %in% parameters && all(c("alpha_true", "alpha_recovered") %in% names(data))) {
    alpha_data <- subset(data, group != 1)  # Exclude non-learners
    
    plots$alpha <- ggplot(alpha_data, aes(x = alpha_true, y = alpha_recovered)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      geom_point(size = 3, alpha = 0.7) +
      labs(title = sprintf("Alpha Recovery (r = %.3f)", recovery$metrics$alpha_correlation),
           x = "True Alpha", y = "Recovered Alpha") +
      coord_fixed() +
      theme_minimal()
  }
  
  # Lambda recovery
  if ("lambda" %in% parameters && all(c("lambda_true", "lambda_recovered") %in% names(data))) {
    lambda_data <- subset(data, group != 1)  # Exclude non-learners
    
    plots$lambda <- ggplot(lambda_data, aes(x = lambda_true, y = lambda_recovered)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      geom_point(size = 3, alpha = 0.7) +
      labs(title = sprintf("Lambda Recovery (r = %.3f)", recovery$metrics$lambda_correlation),
           x = "True Lambda", y = "Recovered Lambda") +
      coord_fixed() +
      theme_minimal()
  }
  
  # Group recovery
  if ("group" %in% parameters && !is.null(recovery$metrics$confusion_matrix)) {
    conf_matrix <- as.data.frame(as.table(recovery$metrics$confusion_matrix))
    names(conf_matrix) <- c("True", "Recovered", "Count")
    
    plots$group <- ggplot(conf_matrix, aes(x = Recovered, y = True, fill = Count)) +
      geom_tile(color = "white") +
      geom_text(aes(label = Count), size = 5) +
      scale_fill_gradient(low = "white", high = "steelblue") +
      labs(title = sprintf("Group Recovery (Accuracy = %.1f%%)", 
                          recovery$metrics$group_accuracy * 100),
           x = "Recovered Group", y = "True Group") +
      coord_fixed() +
      theme_minimal()
  }
  
  if (length(plots) == 1) {
    return(plots[[1]])
  } else if (length(plots) > 1) {
    return(plots)
  } else {
    warning("No valid parameters to plot")
    return(NULL)
  }
}

#' mfeargen Theme
#'
#' @description A clean theme for mfeargen plots
#' @return ggplot2 theme object
#' @export
theme_mfeargen <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 11),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}
