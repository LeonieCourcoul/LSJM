#' Management of data for residual part
#'
#' @param formGroup A formula which indicates the group variable
#' @param formFixed A formula which indicates the fixed effects for the longitudinal submodel
#' @param formRandom A formula which indicates the random effects for the longitudinal submodel
#' @param data.long1 A dataframe with the longitudinal data
#'
#' @importFrom stats model.frame model.matrix

data.manag.sigma <- function(formGroup, formFixed, formRandom, data.long1){

  data_long <- data.long1[unique(c(all.vars(formGroup), all.vars(formFixed), all.vars(formRandom)))]
  #y.new.prog <- data_long[all.vars(formFixed)][, 1]
  mfX <- model.frame(formFixed, data = data_long)
  X <- model.matrix(formFixed, mfX)
  mfU <- model.frame(formRandom, data = data_long)
  U <- model.matrix(formRandom, mfU)
  list.long <- list("X" = X, "U" = U)

  return(list.long)
}
