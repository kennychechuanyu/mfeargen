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
  "participant", "stimulus", "value"
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
  
  # Extract MCMC samples
  mcmc_output <- fit$jags_fit
  
  # Create trace plots
  plot_data <- data.frame()
  
  for (param in parameters) {
    if (param %in% names(mcmc_output$sims.list)) {
      samples <- mcmc_output$sims.list[[param]]
      
      # Sample a few participants for clarity
      n_show <- min(5, ncol(samples))
      
      for (i in 1:n_show) {
        plot_data <- rbind(plot_data, data.frame(
          iteration = 1:nrow(samples),
          value = samples[, i],
          parameter = param,
          participant = i
        ))
      }
    }
  }
  
  if (nrow(plot_data) == 0) {
    warning("No data to plot")
    return(NULL)
  }
  
  # Create trace plot
  p <- ggplot(plot_data, aes(x = iteration, y = value, color = factor(participant))) +
    geom_line(alpha = 0.7) +
    facet_wrap(~ parameter, scales = "free_y") +
    labs(title = "MCMC Trace Plots",
         x = "Iteration",
         y = "Parameter Value",
         color = "Participant") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
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
  
  if (is.null(fit$group_membership)) {
    warning("No group membership data found")
    return(NULL)
  }
  
  # Get group counts - ensure all 4 groups are represented
  if ("group" %in% names(fit$group_membership)) {
    group_counts <- table(factor(fit$group_membership$group, levels = 1:4))
  } else if ("assigned_group" %in% names(fit$group_membership)) {
    # Map group labels to numbers
    group_map <- c("Non-Learner" = 1, "Overgeneralizer" = 2, "Physical" = 3, "Perceptual" = 4)
    group_nums <- group_map[fit$group_membership$assigned_group]
    group_counts <- table(factor(group_nums, levels = 1:4))
  } else {
    # Default if no group info
    group_counts <- rep(0, 4)
  }
  
  group_data <- data.frame(
    group = factor(1:4),
    count = as.numeric(group_counts),
    label = c("Non-Learner", "Overgeneralizer", "Physical", "Perceptual")
  )
  
  # Create bar plot
  p <- ggplot(group_data, aes(x = label, y = count, fill = label)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = count), vjust = -0.5) +
    scale_fill_manual(values = c("#CC79A7", "#F0E442", "#56B4E9", "#D55E00")) +
    labs(title = "Group Membership Distribution",
         x = "Group",
         y = "Number of Participants") +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

#' Plot Parameters
#'
#' @description Creates density plots of parameter distributions
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
  
  # Create density plot
  plot_data <- data.frame(value = param_values)
  
  p <- ggplot(plot_data, aes(x = value)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", alpha = 0.7) +
    geom_density(color = "darkblue", size = 1) +
    labs(title = paste("Distribution of", parameter),
         x = parameter,
         y = "Density") +
    theme_minimal()
  
  return(p)
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
        geom_line(size = 1.2) +
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
    geom_line(size = 1.2) +
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
    lambda_data <- subset(data, group > 1)  # Exclude non-learners
    
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
