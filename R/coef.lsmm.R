#' @export

coef.lsmm <- function(object, ...){


  coef_fit <- object$table.res$Estimation[1:length(object$result_step1$b)]
  names(coef_fit) <- rownames(fit$table.res)[1:length(coef_fit)]

  coef_fit

}
