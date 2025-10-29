#' Dynamic Predictions
#'
#' @description Generic function for dynamic predictions.
#'
#' @param object Model object.
#' @param newdata Optional new data for prediction.
#' @param s Landmark time.
#' @param horizon Horizon time.
#' @param event Event indicator.
#' @param IC Confidence interval level.
#' @param nb.draws Number of Monte Carlo draws.
#' @rdname dynpred
#' @export
dynpred <- function(object, newdata, s, horizon, event, IC = 95,
                    nb.draws = 1000) {
  UseMethod("dynpred")
}
