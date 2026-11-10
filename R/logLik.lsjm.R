#' @export
logLik.lsjm <- function(object){
  if(object$result_step1$istop != 1){
    cat("--------------------------------------------------------------------------------------------------------- \n")
    cat("WARNING", "\n")
    cat("The first step estimation didn't reach convergence because :")
    if(x$result_step1$istop==2) cat("Maximum number of iteration reached without convergence.")
    if(x$result_step1$istop==4) {cat("The program stopped abnormally. No results can be displayed. \n")}
    cat("\n")
    cat("--------------------------------------------------------------------------------------------------------- \n")
    stop("\n")
  }
  ll <- obbject$result_step1$fn.value
  attr(ll, "df") <- length(object$result_step1$b)
  attr(ll, "nobs") <- object$control$Objectlsmm$control$Ind
  class(ll) <- "logLik"

  ll
}
