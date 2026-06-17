#' @export
ranef <- function(object, ...) {
  UseMethod("ranef")
}


#' @export
dynpred <- function(object, ...) {
  UseMethod("dynpred")
}

#' @export
survmarg <- function(object, individual, time) {
  UseMethod("survmarg")
}
