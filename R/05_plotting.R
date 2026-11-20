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
  "participant", "stimulus", "value", "chain", "density", "probability"
))

#' Plot mfeargen Results
#'
#' @description
#' Main plotting function for mfeargen objects. Creates diagnostic and result plots.
#'
#' @param x mfeargen object
#' @param type Character. Plot type: "diagnostic", "parameters", "groups", "gradient", or "all"
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
                            "gradient", "recovery", "all"))
  
  if (type == "diagnostic") {
    return(plot_mcmc_diagnostics(x, ...))
  } else if (type == "parameters") {
    return(plot_parameters(x, ...))
  } else if (type == "groups") {
    return(plot_group_probabilities(x, ...))
  } else if (type == "gradient") {
    return(plot_generalization_gradient(x, ...))
  } else if (type == "recovery" && inherits(x, "mfeargen_recovery")) {
    return(plot_parameter_recovery(x, ...))
  } else if (type == "all") {
    plots <- list(
      diagnostic = plot_mcmc_diagnostics(x, ...),
      parameters = plot_parameters(x, ...),
      groups = plot_group_probabilities(x, ...),
      gradient = plot_generalization_gradient(x, ...)
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

#' Plot Generalization Gradient
#'
#' @description Creates plot showing generalization gradient across stimuli
#' @param fit mfeargen object
#' @param group_separate Logical. Plot groups separately?
#' @return ggplot object
#' @export
plot_generalization_gradient <- function(fit, group_separate = TRUE) {
  
  if (!inherits(fit, "mfeargen")) {
    stop("Input must be a mfeargen object")
  }
  
  data <- fit$data
  if ("phase" %in% names(data)) data <- data[data$phase == "generalization", , drop = FALSE]
  if (!"stimulus" %in% names(data) && "Phy_size" %in% names(data)) data$stimulus <- data$Phy_size
  
  # Calculate mean response by stimulus
  if (!"stimulus" %in% names(data)) {
    warning("Stimulus column not found in data")
    return(NULL)
  }
  
  gradient_data <- aggregate(US_expect ~ stimulus, data = data, mean)
  
  membership <- fit$group_membership
  has_membership <- !is.null(membership) &&
    all(c("participant", "group_label") %in% names(membership))
  
  if (group_separate && has_membership) {
    participant_lookup <- setNames(membership$group_label, membership$participant)
    data$group_label <- participant_lookup[as.character(data$participant)]
    
    if (!any(!is.na(data$group_label))) {
      warning("Group labels could not be matched to participants; plotting overall gradient instead.")
    } else {
      gradient_data <- aggregate(US_expect ~ stimulus + group_label, data = data, mean)
      p <- ggplot(gradient_data, aes(x = stimulus, y = US_expect, color = group_label)) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 3) +
        scale_color_manual(values = c("#CC79A7", "#F0E442", "#56B4E9", "#D55E00")) +
        labs(title = "Generalization Gradients by Group",
             x = "Stimulus",
             y = "Mean US Expectancy",
             color = "Group") +
        theme_minimal()
      return(p)
    }
  }
  
  p <- ggplot(gradient_data, aes(x = stimulus, y = US_expect)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    labs(title = "Overall Generalization Gradient",
         x = "Stimulus",
         y = "Mean US Expectancy") +
    theme_minimal()
  return(p)
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
    lambda_data <- subset(data, group > 2)  # Exclude non-learners
    
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
