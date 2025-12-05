#' @export
ranef <- function(object, ...) {
  UseMethod("ranef")
}


#' @export
dynpred <- function(object, ...) {
  UseMethod("dynpred")
}

#' @export
surv_marg <- function(object, individual, time) {
  UseMethod("surv_marg")
}
