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

#' Regression Uncertainty Modeler
#'
#' @author David Senhora Navega
#'
#' @export
#'
#' @param known a numeric vector. Known values used to train or fit a
#' statitical or machine learning regression model
#' @param predicted a numeric vector. Predicted values obtain from a regression
#' model. Ideally this vector should be obtained from cross-validation.
#' @param type a character. The type of uncertainty model to fit. See Details
#' @param alpha uncertainty level. Used when type = "local".
#' @param interval a numeric vector. This vector should have two elements with
#' the minimum and maximum values with the bounds for the regression function.
#' Default is NULL, which means that the values are inferred from known vector.
#' See Details
#' @param delta a numeric value. Controls the extrapolation outside of known
#' values. Must be between 0 and 1. Default is 0.25. See Details.

#' @param signed a logical. Defines if conformal prediction should be signed
#' (assymetrical) or unsigned (symmetrical) prediction intervals around the
#' predicted value. See Details.
#' @param exponent a numeric. Non-linear variance estimation
#'
#' @details
#' TODO
#' @references
#' TODO
#' @return a rumr class object.
#'
rumr <- function(
  known, predicted, type,
  alpha = 0.05, interval = NULL, delta = 0.25, signed = F, exponent = 1
) {

  type <- match.arg(type, c("gaussian", "conformal", "local", "variance"), F)

  model <- switch(type,

    gaussian = {
      infer_gaussian(known, predicted, interval, delta, exponent)
    },

    conformal = {
      infer_conformal(known, predicted, interval, delta, signed, exponent)
    },

    variance = {
      infer_variance(known, predicted, interval, delta, exponent)
    },

    local = {
      infer_local(known, predicted, interval, delta, alpha, exponent)
    }

  )

  object <- structure(
    .Data = list(
      type = type,
      model = model
    ),
    class = "rumr"
  )

  return(object)

}

#' Predict Method for rumr class.
#'
#' @author David Senhora Navega
#'
#' @export
#'
#' @param object a rumr class object
#' @param predicted a numeric vector. Predicted or fitted values from a
#' a regression model
#' @param alpha numeric value. Uncertainty level for the prediction interval.
#' For instance a 95% prediction level is given by alpha = 0.05. Confidence
#' level is defined as 1 - alpha. alpha must be between 0.01 and 0.5.
#' @param ... ...
#'
predict.rumr <- function(object, predicted, alpha = 0.05, ...) {

  if (isFALSE(inherits(object, "rumr")))
    stop("\n(-) object must a rumr model.")

  switch(object$type,

    gaussian = {
      predict_gaussian(object$model, predicted, alpha)
    },

    conformal = {
      predict_conformal(object$model, predicted, alpha)
    },

    local = {
      predict_local(object$model, predicted, alpha)
    },

    variance = {
      predict_variance(object$model, predicted)
    }

  )

}

#' @author David Senhora Navega
#' @export
#' @noRd
#'
print.rumr <- function(x, ...) {
  cat("rumr: regression uncertainty modeler")
  cat("\ntype:", x$type)
}

#' @author David Senhora Navega
#' @export
#' @noRd
#'
plot.rumr <- function(
  x, predicted, alpha = 0.05, normalize = F, digits = 3, label = NULL, ...
) {

  if (inherits(x, "rumr")) {
    if (x$type != "gaussian")
      stop("\n(-) rumr object type must be 'gaussian'.")

    if (!(is.numeric(predicted) & length(predicted) == 1))
      stop("predicted must a numeric vector of length 1.")

    plot_gaussian(x$model, predicted, alpha, normalize, digits, label)

  } else {
    stop("\n(-) x must be a 'rumr' class object.")
  }

}
