#' @export
ranef <- function(object, ...) {
  UseMethod("ranef")
}

#' @export
plot <- function(object, ...) {
  UseMethod("plot")
}

#' @export
dynpred <- function(object, ...) {
  UseMethod("dynpred")
}
