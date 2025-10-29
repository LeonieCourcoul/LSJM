#' Management of survival data
#'
#' @param formGroup A formula which indicates the group variable
#' @param formSurv A formula which indicates the survival submodel
#' @param data.long1 A dataframe with the longitudinal data
#'
#' @importFrom stats model.frame model.matrix
#'
data.manag.surv <- function(formGroup, formSurv, data.long1){
  tmp <- data.long1[unique(c(all.vars(formGroup),all.vars(formSurv)))]
  tmp <- unique(tmp)
  #Time <- tmp[all.vars(formSurv)][, 1]    # matrix of observed time such as Time=min(Tevent,Tcens)
  #event <- tmp[all.vars(formSurv)][, 2]  # vector of event indicator (delta)
  #event1 <- ifelse(event == 1, 1,0)
  #nTime <- length(Time)                   # number of subject having Time
  #zeros <- numeric(nTime)                 # for zero trick in Bayesian procedure
  # design matrice
  mfZ <- model.frame(formSurv, data = tmp)
  Z <- model.matrix(formSurv, mfZ)


  list("Z" = Z)

}
