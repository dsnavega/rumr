# This file is part of rumr
#
# Copyright (C) 2022, David Senhora Navega
#
# rumr is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# rumr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with rumr. If not, see <http://www.gnu.org/licenses/>.
#
# David Senhora Navega
# Laboratory of Forensic Anthropology
# Department of Life Sciences
# University of Coimbra
# Calçada Martim de Freitas, 3000-456, Coimbra
# Portugal

#' @author David Senhora Navega
#' @noRd
#'
compute_interval <- function(x, delta = 0.25) {

  n <- sum(!is.na(x))
  bw <- sd(x = x, na.rm = T)  * n ^ -0.25
  interval <- range(x = x, na.rm = T) + c(-bw, bw) * delta
  return(interval)

}

#' @author David Senhora Navega
#' @noRd
#'
clamp_value <- function(x, interval) {

  x <- pmin(pmax(x, interval[1], na.rm = TRUE), interval[2], na.rm = TRUE)
  return(x)

}

#' Weighted Quantile Estimation
#'
#' @author David Senhora Navega
#'
#' @noRd
#' @import stats
#'
weighted_quantile <- function (x, probs, weights) {

  if (!is.numeric(x))
    stop("\n(-) x must be a numeric vector.")

  if (missing(weights))
    weights <- rep(1/length(x), length(x))

  if (missing(probs)) {
    stop("\n(-) Please specifiy the target quantiles (probs argument).")
  } else {
    if(any(probs < 0) | any(probs > 1))
      stop("\n(-) Elements of probs must be between 0 and 1.")
  }

  df <- stats::na.omit(data.frame(x = x, weights = weights))
  x <- df$x
  weights <- df$weights

  compute_quantile <- Vectorize(function (x, weights, probs) {

    new.order <- order(x)
    ordered.weights <- weights[new.order]
    inequality <- cumsum(ordered.weights) < (probs * sum(ordered.weights))
    if(sum(inequality) == 0) {
      index <- 1
    } else {
      index <- max(which(inequality))
    }
    x[new.order][index]

  }, vectorize.args = "probs")

  return(compute_quantile(x, weights, probs))

}

#' @author David Senhora Navega
#' @noRd
#'
truncate_gaussian <- function(estimate, variance, interval, alpha) {

  confidence <- c(alpha / 2, 1 - (alpha / 2))
  interval <- sort(range(interval))

  fitted <- truncnorm::qtruncnorm(
    p = confidence, a = interval[1], b = interval[2],
    mean = estimate, sd = variance
  )

  if (estimate < fitted[1])
    fitted[1] <- estimate

  if (estimate > fitted[2])
    fitted[2] <- estimate

  predicted <- c(estimate, fitted)
  predicted <- clamp_value(predicted, interval)
  names(predicted) <- c("estimate", "lower", "upper")

  return(predicted)

}

#' @author David Senhora Navega
#' @noRd
#'
truncate_gaussian <- Vectorize(
  FUN = truncate_gaussian,
  vectorize.args = c("estimate", "variance"),
)

# Variance Estimation Model ----

#' @author David Senhora Navega
#' @noRd
#'
infer_variance <- function(known, predicted, interval, delta, exponent) {

  if (isFALSE(is.vector(known)))
    stop("\n(-) known must be a vector.")

  if (isFALSE(is.numeric(known)))
    stop("\n(-) known must be a numeric")

  if (isFALSE(is.vector(predicted)))
    stop("\n(-) predicted must be a vector.")

  if (isFALSE(is.numeric(predicted)))
    stop("\n(-) predicted must be a numeric")

  if (isTRUE(any(is.na(predicted)) | any(is.na(known))))
    stop("\n(-) NA values not allowed.")

  if (length(predicted) != length(known))
    stop("\n(-) predicted and known must have the same number of observations.")

  if (missing(interval))
    interval <- NULL

  if (missing(delta))
    delta <- 0.25

  if (isTRUE(is.null(interval)))
    interval = compute_interval(x = known, delta = delta)

  interval <- sort(range(interval))

  predicted <- clamp_value(x = predicted, interval = interval)

  A <- cbind(intercept = 1, predicted ^ exponent)
  b <- cbind(variance = abs(known - predicted))

  # Compute coefficients (Beta, B)
  B <- solve(t(A) %*% A, t(A) %*% b)

  # Hat Matrix (Diagonal)
  H <- diag(A %*% solve(t(A) %*% A) %*% t(A))

  # LOOCV
  fitted <- ((A %*% B) - (b * H)) / (1.0 - H)

  object <- list(
    coefficients = B,
    exponent = exponent,
    interval = interval,
    loocv = pmax(x = 1.2533 * as.vector(fitted), 0)
  )

  invisible(object)

}

#' @author David Senhora Navega
#' @noRd
#'
predict_variance <- function(object, predicted) {


  if (missing(predicted))
    return(object$loocv)

  if (isFALSE(is.vector(predicted)))
    stop("\n(-) predicted must be a vector.")

  if (isFALSE(is.numeric(predicted)))
    stop("\n(-) predicted must be a numeric")

  if (isTRUE(any(is.na(predicted))))
    stop("\n(-) NA values not allowed.")

  predicted <- clamp_value(predicted, object$interval)
  estimate <- cbind(1, predicted ^ object$exponent) %*% object$coefficients
  estimate <- pmax(as.vector(1.2533 * estimate), 0.0)

  return(estimate)

}

# Conformal Prediction Estimation Model ----

#' @author David Senhora Navega
#' @noRd
#'
infer_conformal <- function(known, predicted, interval, delta, signed, exponent) {

  if (isFALSE(is.vector(known)))
    stop("\n(-) known must be a vector.")

  if (isFALSE(is.numeric(known)))
    stop("\n(-) known must be a numeric")

  if (isFALSE(is.vector(predicted)))
    stop("\n(-) predicted must be a vector.")

  if (isFALSE(is.numeric(predicted)))
    stop("\n(-) predicted must be a numeric")

  if (isTRUE(any(is.na(predicted)) | any(is.na(known))))
    stop("\n(-) NA values not allowed.")

  if (length(predicted) != length(known))
    stop("\n(-) predicted and known must have the same number of observations.")

  if (missing(signed))
    signed <- TRUE

  if (isFALSE(is.logical(signed)))
    stop("\n(-) signed must be a logical (TRUE / FALSE).")

  if (missing(interval))
    interval <- NULL

  if (missing(delta))
    delta <- 0.25

  if (isTRUE(is.null(interval)))
    interval = compute_interval(x = known, delta = delta)

  interval <- sort(range(interval))

  predicted <- clamp_value(x = predicted, interval = interval)

  # Infer Variance
  variance_model <- infer_variance(
    known = known, predicted = predicted,
    interval = interval, delta = delta,
    exponent = exponent
  )

  variance <- predict_variance(variance_model)

  seq_alpha <- seq(0.01, 0.5, 0.01)
  confidence_list <- lapply(seq_alpha, function(alpha) {
    if (signed) {
      c(alpha / 2, 1 - (alpha / 2))
    } else {
      1 - alpha
    }
  })

  confidence_vector <- sort(unlist(confidence_list))
  names(confidence_vector) <- as.character(confidence_vector)

  if (signed) {
    residual <- known - predicted
  } else {
    residual <- abs(known - predicted)
  }

  conformal_vector <- sort(residual / variance)
  n_vector <- length(conformal_vector)
  conformal_score <- conformal_vector[ceiling(n_vector * confidence_vector)]
  names(conformal_score) <- as.character(confidence_vector)

  object <- list(
    model = variance_model,
    signed = signed,
    interval = interval,
    scale = conformal_score,
    predicted = predicted
  )

  invisible(object)

}

#' @author David Senhora Navega
#' @noRd
#'
predict_conformal <- function(object, predicted, alpha = 0.1) {

  if (missing(predicted)) {

    predicted <- object$predicted
    variance <- predict_variance(object$model)

  } else {

    if (isFALSE(is.vector(predicted)))
      stop("\n(-) predicted must be a vector.")

    if (isFALSE(is.numeric(predicted)))
      stop("\n(-) predicted must be a numeric")

    if (isTRUE(any(is.na(predicted))))
      stop("\n(-) NA values not allowed.")

    predicted <- clamp_value(predicted, object$interval)
    variance <- predict_variance(object$model, predicted)

  }

  if (alpha > 0.5 | alpha < 0.01)
    stop("\n(-) alpha must be a value between 0.01 and 0.5")

  alpha <- round(x = alpha, digits = 2)

  if (object$signed) {

    lwr_confidence <- as.character(alpha / 2)
    upr_confidence <- as.character(1 - (alpha / 2))
    lower <- predicted + (object$scale[lwr_confidence] * variance)
    upper <- predicted + (object$scale[upr_confidence] * variance)

  } else {

    confidence <- as.character(1 - alpha)
    lower <- predicted - (object$scale[confidence] * variance)
    upper <- predicted + (object$scale[confidence] * variance)

  }

  predicted <- cbind(predicted, lower, upper)
  colnames(predicted) <- c("estimate", "lower", "upper")

  predicted <- clamp_value(predicted, object$interval)
  rownames(predicted) <- NULL

  return(predicted)

}

# Local Uncertainty Estimation Model ----

#' David Senhora Navega
#' @noRd
#' @importFrom e1071 cmeans
#' @importFrom stats kmeans
#'
infer_local <- function(known, predicted, interval, delta, alpha, exponent) {

  if (isFALSE(is.vector(known)))
    stop("\n(-) known must be a vector.")

  if (isFALSE(is.numeric(known)))
    stop("\n(-) known must be a numeric")

  if (isFALSE(is.vector(predicted)))
    stop("\n(-) predicted must be a vector.")

  if (isFALSE(is.numeric(predicted)))
    stop("\n(-) predicted must be a numeric")

  if (isTRUE(any(is.na(predicted)) | any(is.na(known))))
    stop("\n(-) NA values not allowed.")

  if (length(predicted) != length(known))
    stop("\n(-) predicted and known must have the same number of observations.")


  if (alpha > 0.5 | alpha < 0.01)
    stop("\n(-) alpha must be a value between 0.01 and 0.5")
  alpha <- round(x = alpha, digits = 2)

  if (missing(interval))
    interval <- NULL

  if (missing(delta))
    delta <- 0.25

  if (isTRUE(is.null(interval)))
    interval = compute_interval(x = known, delta = delta)

  interval <- sort(range(interval))

  predicted <- clamp_value(x = predicted, interval = interval)

  # Residual
  residual <- known - predicted

  # Assess number of clusters from known values
  ncluster <- ceiling(log2(length(unique(known)) + 1)) + 1

  # Fuzzy clustering
  clst_df <- data.frame(known, predicted, error = 1.2533 * abs(residual))
  fcluster <- e1071::cmeans(
    x = scale(clst_df), centers = ncluster, iter.max = 100,
  )

  U <- fcluster$membership

  # Local Uncertainty Estimation Model
  Q <- t(sapply(seq_len(ncluster), function(ith) {
    confidence <- c(alpha / 2, 1 - (alpha / 2))
    weighted_quantile(x = residual, probs = confidence, weights = U[, ith])
  }))

  UQ <- U %*% Q
  colnames(UQ) <- c("lower", "upper")

  # Compute coefficients (Beta, B)
  A <- cbind(intercept = 1, predicted ^ exponent)
  B <- solve(t(A) %*% A, t(A) %*% UQ)

  # Hat Matrix (Diagonal)
  H <- diag(A %*% solve(t(A) %*% A) %*% t(A))

  # LOOCV
  endpoints <- ((A %*% B) - (UQ * H)) / (1.0 - H)
  estimate <- cbind(estimate = predicted, predicted + endpoints)
  estimate <- clamp_value(estimate, interval = interval)

  local_model <- list(
    coefficients = B,
    exponent = exponent,
    predicted = estimate,
    alpha = alpha
  )

  object <- list(
    model = local_model,
    interval = interval
  )

  invisible(object)

}

#' @author David Senhora Navega
#' @noRd
#'
predict_local <- function(object, predicted, alpha) {

  if (alpha > 0.5 | alpha < 0.01)
    stop("\n(-) alpha must be a value between 0.01 and 0.5")

  alpha <- round(x = object$model$alpha, digits = 2)

  if (missing(predicted)) {
    predicted <- object$model$predicted
  } else {

    if (isFALSE(is.vector(predicted)))
      stop("\n(-) predicted must be a vector.")

    if (isFALSE(is.numeric(predicted)))
      stop("\n(-) predicted must be a numeric")

    if (isTRUE(any(is.na(predicted))))
      stop("\n(-) NA values not allowed.")

    predicted <- clamp_value(x = predicted, interval = object$interval)
    local <- object$model
    endpoints <- cbind(intercept = 1, predicted ^ local[[2]]) %*% local[[1]]
    predicted <- cbind(estimate = predicted, predicted + endpoints)
    predicted <- clamp_value(x = predicted, interval = object$interval)

  }

  return(predicted)

}

# Truncated Gaussian Uncertainty Model ----

#' @author David Senhora Navega
#' @noRd
#' @import truncnorm
#'
infer_gaussian <- function(known, predicted, interval, delta, exponent) {

  if (isFALSE(is.vector(known)))
    stop("\n(-) known must be a vector.")

  if (isFALSE(is.numeric(known)))
    stop("\n(-) known must be a numeric")

  if (isFALSE(is.vector(predicted)))
    stop("\n(-) predicted must be a vector.")

  if (isFALSE(is.numeric(predicted)))
    stop("\n(-) predicted must be a numeric")

  if (isTRUE(any(is.na(predicted)) | any(is.na(known))))
    stop("\n(-) NA values not allowed.")

  if (length(predicted) != length(known))
    stop("\n(-) predicted and known must have the same number of observations.")

  if (missing(interval))
    interval <- NULL

  if (missing(delta))
    delta <- 0.25

  if (isTRUE(is.null(interval)))
    interval = compute_interval(x = known, delta = delta)

  interval <- sort(range(interval))

  predicted <- clamp_value(x = predicted, interval = interval)

  # Infer Variance
  variance_model <- infer_variance(
    known = known, predicted = predicted,
    interval = interval, delta = delta,
    exponent = exponent
  )

  object <- list(
    model = variance_model,
    interval = interval,
    predicted = predicted
  )

  invisible(object)

}

#' @author David Senhora Navega
#' @noRd
#'
predict_gaussian <- function(object, predicted, alpha = 0.1) {

  if (missing(predicted)) {

    predicted <- object$predicted
    variance <- predict_variance(object$model)

  } else {

    if (isFALSE(is.vector(predicted)))
      stop("\n(-) predicted must be a vector.")

    if (isFALSE(is.numeric(predicted)))
      stop("\n(-) predicted must be a numeric")

    if (isTRUE(any(is.na(predicted))))
      stop("\n(-) NA values not allowed.")

    predicted <- clamp_value(predicted, object$interval)
    variance <- predict_variance(object$model, predicted)

  }

  if (alpha > 0.5 | alpha < 0.01)
    stop("\n(-) alpha must be a value between 0.01 and 0.5")

  alpha <- round(x = alpha, digits = 2)

  fitted <- t(
    truncate_gaussian(predicted, variance, object$interval, alpha)
  )

  return(fitted)

}

#' @author David Senhora Navega
#' @noRd
#' @import ggplot2
#' @importFrom stats supsmu
#'
plot_gaussian <- function(
    object, predicted, alpha = 0.1, normalize = T, digits = 3, label = NULL
) {

  predicted <- clamp_value(predicted, object$interval)
  fitted <- predict_gaussian(object, predicted, alpha = alpha)
  variance <- predict_variance(object$model, predicted)

  # Prediction Domain
  x <- seq(
    from = object$interval[1], to = object$interval[2], length.out = 2048
  )

  x <- sort(c(fitted, x), decreasing = FALSE)
  x <- round(x, digits = digits)

  # Gaussian Density
  y <- truncnorm::dtruncnorm(
    x = x, a = object$interval[1], b = object$interval[2],
    mean = fitted[1], sd = variance
  )
  y <- round(y, digits = digits)

  # Super smooth data using Friedman's supsmu for prettier plot
  if (normalize) {
    df <- data.frame(x = c(min(x), x, max(x)), y = c(0, y / max(y), 0))
    supsmu_df <- stats::supsmu(df$x, df$y)
    df <- data.frame(x = supsmu_df$x, y = supsmu_df$y)
    df <- data.frame(
      x = c(min(df$x), df$x, max(df$x)), y = c(0, df$y / max(df$y), 0)
    )
    y_label <- "Density (Normalized)"
  } else {
    df <- data.frame(x = c(min(x), x, max(x)), y = c(0, y, 0))
    supsmu_df <- stats::supsmu(df$x, df$y)
    df <- data.frame(x = supsmu_df$x, y = supsmu_df$y)
    df <- data.frame(x = c(min(df$x), df$x, max(df$x)), y = c(0, df$y, 0))
    y_label <- "Density"
  }

  x_values <- round(fitted, digits = digits)

  subset_df <- df[df$x >= x_values[2] & df$x <= x_values[3],]

  lollipop_df <- data.frame(
    x = subset_df$x[subset_df$x %in% x_values],
    y = subset_df$y[subset_df$x %in% x_values]
  )

  x_values <- format(x_values, nsmall = digits)
  plt <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line(size = 0.75) +
    ggplot2::geom_area(data = subset_df, alpha = 0.25, linetype = "dashed") +
    ggplot2::geom_point(data = lollipop_df, size = 5, shape = 19) +
    ggplot2::geom_segment(
      data = lollipop_df, linetype = 2, size = 1,
      mapping = ggplot2::aes(x = x, xend = x, y = 0, yend = y)
    ) +
    ggplot2::scale_x_continuous(
      name = label,
      limits = range(x),
      breaks = seq(min(x) - (min(x) %% 10), max(x), 10),
      minor_breaks = seq(min(x) - (min(x) %% 5), max(x), 5),
    ) +
    ggplot2::ylab(label = y_label) +
    ggplot2::ggtitle(
      label = "Predictive Distribution (Truncated Gaussian)",
      subtitle = paste0(
        "Predicted: ",x_values[1]," [", x_values[2]," - ",x_values[3],"]\n",
        "Variance: ", round(variance, digits = digits),"\n",
        "Confidence: ", 1 - alpha
      )
    ) +
    ggplot2::theme_classic(base_line_size = 0.75) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 13),
      axis.ticks = ggplot2::element_line(
        size = 1, colour = "black", lineend = "round",linetype = 3
      ),
      axis.title = ggplot2::element_text(
        family = "sans", face = "plain", size = 14
      ),
      axis.text = ggplot2::element_text(
        family = "sans",face = "plain", size = 11, colour = "black"
      )
    )

  return(plt)

}
