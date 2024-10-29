#' lsmm : Estimation of linear mixed model for longitudinal data with a flexible subject-specific variability.
#'
#' This function fits linear mixed effects models in which
#' we suppose that the variance of the residual error is subject-specific. Three differents models can be estimated (see details below)
#' Parameters are estimated through a maximum likelihood method, using a Marquardt-Levenberg algorithm.
#'
#'
#'
#' @param formFixed A formula for the fixed effects of the longitudinal submodel
#' @param formRandom A formula for the random effects of the longitudinal submodel
#' @param formGroup A formula which indicates the group variable
#' @param timeVar A character which indicates the time variable
#' @param formVar A character : type of variability either 'classic' for a classical linear mixed model or 'cov-dependent' for a subject-specific and covariate-dependent covariable variability or 'inter-intra' for distinguishing inter-visit from intra-visit variability
#' @param formFixedVar A formula for the fixed effects of the variance predictor if formVar == 'cov-dependent'
#' @param formRandomVar A formula for the random effects of the variance predictor if formVar == 'cov-dependent'
#' @param random_inter A logical indicating if the inter-visits variability is subject-specific when formVar = 'inter-intra'
#' @param random_intra A logical indicating if the intra-visit variability is subject-specific when formVar = 'inter-intra'
#' @param formGroupVisit A formula which indicates the visit indicator variable  when formVar = 'inter-intra'
#' @param correlated_re A logical indicating if the random effects of the trend are correlated to the random effects of the variability when formVar is in c('cov-dependent', 'inter-intra')
#' @param data.long A dataframe with the longitudinal data
#' @param S1 An integer : the number of QMC draws for the first step
#' @param S2 An integer : the number of QMC draws for the second step
#' @param nproc An integer : the number of processors for parallel computing
#' @param clustertype one of the supported types from \code{makeCluster} function
#' @param maxiter optional maximum number of iterations for the marqLevAlg iterative algorithm.
#' @param print.info logical indicating if the outputs of each iteration should be written
#' @param file optional character giving the name of the file where the outputs of each iteration should be written (if print.info=TRUE)
#' @param epsa optional threshold for the convergence criterion based on the parameter stability.
#' @param epsb optional threshold for the convergence criterion based on the objective function stability.
#' @param epsd  optional threshold for the relative distance to maximum. This criterion has the nice interpretation of estimating the ratio of the approximation error over the statistical error, thus it can be used for stopping the iterative process whatever the problem.
#' @param binit optional initials parameters.
#'
#' @return A lsmm object which contains the following elements :
#' \describe{
#' \item{\code{table.res}}{The table of results : Estimation and SE}
#' \item{\code{result_step1}}{A marqLevAlg object with the results of the first step estimation.}
#' \item{\code{result_step1}}{A marqLevAlg object with the results of the second step estimation.}#'
#' \item{\code{info_conv_step2}}{Information about convergence (criteria and time)}
#' \item{\code{control}}{A list of control elements}
#'
#' }
#'
#' @export
#'
#' @examples
#'
#' \donttest{
#'
#' First example : a linear mixed model with subject-specific time-dependent and covariate-dependent variability
#'
#' Second example : a linear mixed model with subject-specific inter-visits and intra-visits variabilities
#'
#' Third example : a linear mixed model with subject-specific inter-visits variability and constant intra-visit variability
#'
#' }
#'
#'
#'
#'
#'
#'
lsmm <- function(formFixed, formRandom, formGroup, timeVar,
                 formVar = "classic", formFixedVar = NULL, formRandomVar = NULL,
                 random_inter = F, random_intra = F, formGroupVisit = NULL, correlated_re = F,
                 data.long,
                 S1 = 500, S2= 5000,
                 nproc = 1, clustertype = "SOCK", maxiter = 100, print.info = FALSE,
                 file = "", epsa = 1e-04, epsb = 1e-04, epsd = 1e-04, binit = NULL){

  if(missing(formFixed)) stop("The argument formFixed must be specified")
  if(missing(formRandom)) stop("The argument formRandom must be specified")
  if(!inherits((formFixed),"formula")) stop("The argument formFixed must be a formula")
  if(!inherits((formRandom),"formula")) stop("The argument formRandom must be a formula")
  if(missing(formGroup)) stop("The argument formGroup must be specified")
  if(!inherits((formGroup),"formula")) stop("The argument formGroup must be a formula")
  if(missing(data.long)) stop("The argument data.long must be specified")
  if(!inherits((data.long),"data.frame")) stop("The argument data.long must be a data frame")
  if(nrow(data.long) == 0) stop("Data should not be empty")

  if(missing(timeVar)) stop("The argument timeVar must be specified")
  if(!inherits((timeVar),"character")) stop("The argument timeVar must be a character")
  if(length(timeVar) != 1) stop("The argument timeVar must be of length 1")
  if(!(timeVar %in% colnames(data.long))) stop("Unable to find variable 'timeVar' in 'data.long'")


  if(missing(formVar)) stop("The argument formVar must be specified")
  if(!inherits((formVar),"character")) stop("The argument formVar must be a character")
  if(length(formVar) != 1 ||!(formVar %in% c("classic", "cov-dependent", "inter-intra"))) stop ("The argument formVar must be of lenght 1 and must be 'classic' or 'cov-dependent' or 'inter-intra'")
  if(formVar == "cov-dependent" && missing(formFixedVar)) stop("The argument formFixedVar must be specified")
  if(formVar == "cov-dependent" && missing(formRandomVar)) stop("The argument formRandomVar must be specified")
  if(formVar == "cov-dependent" && !inherits((formFixedVar),"formula")) stop("The argument formFixedVar must be a formula")
  if(formVar == "cov-dependent" && !inherits((formRandomVar),"formula")) stop("The argument formRandomVar must be a formula")
  if(formVar == "inter-intra" && missing(random_inter)) stop("The argument random_inter must be specified")
  if(formVar == "inter-intra" && missing(random_intra)) stop("The argument random_intra must be specified")
  if(formVar == "inter-intra" && !inherits((random_inter),"logical")) stop("The argument random_inter must be a logical")
  if(formVar == "inter-intra" && !inherits((random_intra),"logical")) stop("The argument random_intra must be a logical")
  if(formVar == "inter-intra" && missing(formGroupVisit)) stop("The argument formGroupVisit must be specified")
  if(formVar == "inter-intra" && !inherits((formGroupVisit),"formula")) stop("The argument formGroupVisit must be a formula")
  if(formVar %in% c("cov-dependent", "inter-intra") && !inherits((correlated_re),"logical")) stop("The argument correlated_re must be a logical")

  if(!inherits((S1),"numeric")) stop("The argument S1 must be a numeric")
  if(!is.null(S2) && !inherits((S2),"numeric")) stop("The argument S2 must be a numeric")

  #if(!(all.vars(formFixed) %in% colnames(data.long))) stop("All variables used in the argument formFixed must be in data.long")
  #if(!(all.vars(formRandom) %in% colnames(data.long))) stop("All variables used in the argument formRandom must be in data.long")
  #if(!(all.vars(formGroup) %in% colnames(data.long))) stop("All variables used in the argument formGroup must be in data.long")
  #if(formVar == "cov-dependent" &&!(all.vars(formFixedVar) %in% colnames(data.long))) stop("All variables used in the argument formFixedVar must be in data.long")
  #if(formVar == "cov-dependent" &&!(all.vars(formRandomVar) %in% colnames(data.long))) stop("All variables used in the argument formRandomVar must be in data.long")
  #if(formVar == "inter-intra" &&!(all.vars(formGroupVisit) %in% colnames(data.long))) stop("All variables used in the argument formGroupVisit must be in data.long")
  #browser()
  time.prog1 <- Sys.time()
  data.long <- as.data.frame(data.long)
  id <- as.integer(data.long[all.vars(formGroup)][,1])
  if(!("id" %in% colnames(data.long))){
    data.long <- cbind(data.long, id = id)
  }
  else{
    data.long$id <- as.integer(data.long$id)
  }
  idVar <- "id"
  list.long <- data.manag.long(formGroup,formFixed, formRandom,data.long)
  y.new_glob <- list.long$y.new
  data.long <- as.data.frame(data.long)
  data.long$y.new <-  y.new_glob
  data.long <- as.data.frame(data.long)

  browser()
  if(formVar == "classic"){
    lsmm.result <- lsmm_classic(formFixed, formRandom, formGroup, data.long, idVar, list.long, time.prog1,S1 , S2, nproc , clustertype, maxiter, print.info ,file, epsa, epsb, epsd, binit)
  }
  else{
    if(formVar == "cov-dependent"){
      lsmm.result <- lsmm_covDep(formFixed, formRandom, formGroup, formFixedVar, formRandomVar,correlated_re,data.long, idVar,list.long,time.prog1,S1 , S2,nproc , clustertype, maxiter, print.info ,file, epsa, epsb, epsd, binit)
    }
    else{
      lsmm.result <- lsmm_interintra(formFixed, formRandom, formGroup, random_inter, random_intra, formGroupVisit,correlated_re , data.long, idVar, list.long, time.prog1,S1 , S2,nproc , clustertype, maxiter, print.info ,file, epsa, epsb, epsd, binit)
    }
  }
  lsmm.result$control$timeVar <- timeVar

  lsmm.result

}
