#' lsmm : Estimation of a linear mixed model for longitudinal data with a flexible subject-specific variability.
#'
#' This function fits linear mixed effects models in which
#' we can suppose that the variance of the residual error is subject-specific. Three differents models can be estimated (see details below).
#' Parameters are estimated through a maximum likelihood method, using a Marquardt-Levenberg algorithm.
#'
#' @details
#'
#'
#' The model is defined by:
#' \eqn{Y_{ij} = Y_{i}(t_{ij}) = \widetilde{Y}_i(t_{ij}) + \epsilon_{ij} = X_{ij}^{\top} \beta+Z_{ij}^{\top} b_{i}+\epsilon_{ij}},
#' where $X_{ij}$ and $Z_{ij}$ vectors of explanatory variables for subject $i$ at visit $j$, respectively associated with the fixed-effect vector \eqn{\beta} and the subject-specific random-effect vector \eqn{b_i}.
#'
#' 1. Standard linear mixed model:
#' In this case, \eqn{b_i \sim \mathcal{N}(0,B)}, with \eqn{B} an unspecified matrix and the measurements error \eqn{\epsilon_{ij}} are independent Gaussian errors with variance \eqn{\sigma^2_{\epsilon}}.
#'
#' 2. Location-scale mixed model with time and/or covariate-dependent variability:
#' In this model we assume the following specification for the residual error:
#' \eqn{\epsilon_{ij} \sim \mathcal{N}(0,\sigma_i^2)} with \eqn{ \log(\sigma_i(t_{ij}))  = O_{ij}^{\top} \mu+M_{ij}^{\top} \tau_{i}}.
#'  where $O_{ij}$ and $M_{ij}$ vectors of explanatory variables for subject $i$ at visit $j$, respectively associated with the fixed-effect vector \eqn{\mu} and the subject-specific random-effect vector \eqn{\tau_i}.
#' For the random effects we assume :
#' \eqn{\quad\left(\begin{array}{c}
#'              b_{i} \\
#'              \tau_i
#'              \end{array}\right) \sim N\left(\left(\begin{array}{c}
#'                                                   0 \\
#'                                                   0
#'                                                   \end{array}\right),\left(\begin{array}{cc}
#'                                                                            \Sigma_{b} & \Sigma_{\tau b} \\
#'                                                                            \Sigma_{\tau b}' & \Sigma_{\tau}
#' \end{array}\right)\right)}
#'
#' As conventionally assumed in practice, the random effects \eqn{b_i} can be considered independent of the errors by setting \eqn{\Sigma_{\tau b} =0}.
#'
#' 3. Location-scale mixed model distinguishing within and between visits variabilities:
#' In some studies, multiple measurements of a marker are collected at each time point.
#' For example, in medical research, 2 or 3 blood pressure readings are typically taken per visit,
#' and intra-visit variability can be informative. To capture this, we propose an LSMM that distinguishes within- and between-visit variabilities.
#' We introduce an additional level in longitudinal data, grouping repeated measurements by time.
#' For each subject \eqn{i} \eqn{(i=1,...,N)}, \eqn{Y_{ijl}} represents the \eqn{l}-th \eqn{(l=1,...,n_{ij})} measurement at visit \eqn{j} \eqn{(j=1,...,n_i)} and time \eqn{t_{ij}}.
#' We then define the following LSMM to decompose individual residual variance within and between visits:
#' \eqn{
#' \left\{
#'   \begin{array}{ll}
#'   Y_{ijl} =  \widetilde{Y}_i(t_{ij}) + \epsilon_{ij} + \nu_{ijl} = X_{ij}^{\top} \beta+Z_{ij}^{\top} b_{i}+\epsilon_{ij} + \nu_{ijl}, \\
#'   \epsilon_{ij} \sim \mathcal{N}(0,\sigma_i^2) \hspace{4mm} \text{with} \hspace{3mm} \log(\sigma_i)  = \mu_\sigma + \tau_{\sigma i},\\
#'   \nu_{ijl} \sim \mathcal{N}(0,\kappa_i^2) \hspace{3mm} \text{with} \hspace{3mm} \log(\kappa_i)  = \mu_\kappa + \tau_{\kappa i},\\
#'   \end{array}
#'   \right.}
#'
#' The parameter \eqn{\mu_\sigma} and \eqn{\mu_\kappa} are the fixed intercepts for the between-visits and within-visit variances respectively.
#' The subject-specific random-effect \eqn{b_i} and \eqn{\tau_i = (\tau_{\sigma i},\tau_{\kappa i})^\top} are assumed to be Gaussian as presented previously.
#'
#'
#'
#'
#' @param formFixed A formula for the fixed effects of the longitudinal submodel
#' @param formRandom A formula for the random effects of the longitudinal submodel
#' @param formGroup A formula which indicates the group variable
#' @param timeVar A character which indicates the time variable
#' @param formVar A character : type of variability either 'standard' for a standard linear mixed model or 'cov-dependent' for a subject-specific and covariate-dependent covariable variability or 'inter-intra' for distinguishing inter-visit from intra-visit variability
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
#' First example : a standard linear mixed model (constant residual variance, in time and between subjects)
#'
#'
#' Second example : a linear mixed model with subject-specific time-dependent and covariate-dependent variability
#'
#' Third example : a linear mixed model with subject-specific inter-visits and intra-visits variabilities
#'
#' Fourth example : a linear mixed model with subject-specific inter-visits variability and constant intra-visit variability
#'
#' }
#'
#'
#'
#'
#'
#'
lsmm <- function(formFixed, formRandom, formGroup, timeVar,
                 formVar = "standard", formFixedVar = NULL, formRandomVar = NULL,
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
  if(length(formVar) != 1 ||!(formVar %in% c("standard", "cov-dependent", "inter-intra"))) stop ("The argument formVar must be of lenght 1 and must be 'standard' or 'cov-dependent' or 'inter-intra'")
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

  if(formVar == "standard"){
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
