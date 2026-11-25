
# Constructor for the S3 result class.
new_rpmmd_result <- function(x) {

  x <- structure(x, class = "rpmmd_result")

  validate_rpmmd_result(x)
}

# Validate the internal structure of the result.
validate_rpmmd_result <- function(x) {

  stopifnot(
    inherits(x, "rpmmd_result"),
    x$statistic %in% c("cvm", "ks"),
    length(x$combined_statistic) == 1,
    is.finite(x$combined_statistic),
    length(x$mmd_by_projection) == x$nproj,
    all(is.finite(x$mmd_by_projection)),
    nrow(x$directions) == x$nproj,
    ncol(x$directions) == x$dimension,
    length(x$bandwidths) == x$nproj,
    all(is.finite(x$bandwidths)),
    all(x$bandwidths > 0),
    length(x$perm_statistics) == x$nperm,
    all(is.finite(x$perm_statistics)),
    length(x$p_value) == 1,
    is.finite(x$p_value),
    x$p_value >= 1 / (x$nperm + 1),
    x$p_value <= 1
  )

  x
}

#' Print a random-projection MMD test result
#'
#' @param x An object of class \code{rpmmd_result}.
#' @param ... Additional arguments.
#'
#' @return The original object, invisibly.
#' @export
print.rpmmd_result <- function(x, ...) {

  cat("\nRandom-projection MMD two-sample test\n\n")

  cat("Statistic type:", toupper(x$statistic), "\n")
  cat("Combined statistic:", format(x$combined_statistic), "\n")
  cat("Permutation p-value:", format(x$p_value), "\n")
  cat("Sample sizes:", x$n, "and", x$m, "\n")
  cat("Dimension:", x$dimension, "\n")
  cat("Projections:", x$nproj, "\n")
  cat("Permutations:", x$nperm, "\n")
  cat("Engine:", x$engine, "\n")

  invisible(x)
}

#' Summarize a random-projection MMD test result
#'
#' @param object An object of class \code{rpmmd_result}.
#' @param ... Additional arguments.
#'
#' @return A list of summary statistics.
#' @export
summary.rpmmd_result <- function(object, ...) {

  list(
    statistic = object$statistic,
    combined_statistic = object$combined_statistic,
    p_value = object$p_value,
    sample_sizes = c(object$n, object$m),
    dimension = object$dimension,
    nproj = object$nproj,
    nperm = object$nperm,
    engine = object$engine,
    projection_summary = summary(
      object$mmd_by_projection
    ),
    permutation_quantiles = stats::quantile(
      object$perm_statistics,
      probs = c(0.025, 0.5, 0.95, 0.975)
    )
  )
}

#' Plot a random-projection MMD test result
#'
#' @param x An object of class \code{rpmmd_result}.
#' @param main Plot title.
#' @param xlab Horizontal axis label.
#' @param col Histogram color.
#' @param border Histogram border color.
#' @param xlim Horizontal axis limits.
#' @param ... Additional arguments passed to the histogram plotting method.
#'
#' @return The original object, invisibly.
#' @export

plot.rpmmd_result <- function(
    x,
    main = "Permutation distribution",
    xlab = "Combined MMD statistic",
    col = "#6495B5",
    border = "white",
    xlim = NULL,
    ...
) {

  # Compute histogram without drawing it.
  h <- graphics::hist(
    x$perm_statistics,
    breaks = "Sturges",
    plot = FALSE
  )

  # Include both histogram boundaries and observed statistic.
  if (is.null(xlim)) {

    limits <- range(
      c(h$breaks, x$combined_statistic)
    )

    padding <- 0.06 * diff(limits)

    if (padding == 0) {
      padding <- 0.05 * max(1, abs(limits[1]))
    }

    xlim <- limits + c(-padding, padding)
  }

  # Preserve the user's graphical settings.
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::par(
    mar = c(5, 4.5, 4.5, 1.5),
    las = 1,
    bty = "l",
    cex.axis = 0.9,
    cex.lab = 1.05,
    cex.main = 1.25
  )

  graphics::plot(
    h,
    freq = TRUE,
    main = main,
    xlab = xlab,
    ylab = "Frequency",
    col = col,
    border = border,
    xlim = xlim,
    ...
  )

  graphics::abline(
    v = x$combined_statistic,
    col = "#B22222",
    lwd = 2.5,
    lty = 2
  )

  # Display the statistical results above the histogram.
  graphics::mtext(
    sprintf(
      "%s | Observed = %.4f | p-value = %.4f",
      toupper(x$statistic),
      x$combined_statistic,
      x$p_value
    ),
    side = 3,
    line = 0.4,
    cex = 0.85,
    col = "#444444"
  )

  invisible(x)
}

