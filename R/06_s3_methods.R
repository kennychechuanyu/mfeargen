#' Print Method for Recovery Results
#'
#' @param x mfeargen_recovery object
#' @param ... Additional arguments
#' @export
print.mfeargen_recovery <- function(x, ...) {
  cat("mfeargen Parameter Recovery Analysis\n")
  cat("====================================\n")
  cat("Participants:", x$settings$n_participants, "\n")
  cat("Model type:", x$settings$model_type, "\n")
  cat("MCMC iterations:", x$settings$n.iter, "\n\n")
  
  if (!is.null(x$metrics$alpha_correlation)) {
    cat("Recovery Correlations:\n")
    cat(sprintf("  Alpha:  r = %.3f\n", x$metrics$alpha_correlation))
    cat(sprintf("  Lambda: r = %.3f\n", x$metrics$lambda_correlation))
  }
  
  if (!is.null(x$metrics$group_accuracy)) {
    cat(sprintf("\nGroup Accuracy: %.1f%%\n", x$metrics$group_accuracy * 100))
  }
  
  invisible(x)
}
