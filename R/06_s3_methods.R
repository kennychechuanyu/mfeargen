#' Print Parameter Recovery Results
#' @param x An mfeargen_recovery object.
#' @param ... Additional arguments (ignored).
#' @return The input, invisibly.
#' @export
print.mfeargen_recovery <- function(x, ...) {
  cat("Parameter recovery (", x$settings$model_type, ")\n", sep = "")
  cat(x$settings$n_participants, "participants;", x$settings$n.iter, "iterations per chain\n")
  parameters <- c("alpha", "lambda", "w0", "w1")
  table <- data.frame(parameter = parameters)
  for (metric in c("correlation", "rmse", "coverage_95")) {
    table[[metric]] <- vapply(parameters, function(parameter) {
      x$metrics[[paste0(parameter, "_", metric)]] %||% NA_real_
    }, numeric(1))
  }
  print(table, row.names = FALSE)
  if (!is.null(x$metrics$group_accuracy)) cat("Group accuracy:", round(100 * x$metrics$group_accuracy, 1), "%\n")
  if (!isTRUE(x$fit$convergence$converged)) {
    cat("Sampling diagnostics require review; recovery metrics are provisional.\n")
  }
  invisible(x)
}
