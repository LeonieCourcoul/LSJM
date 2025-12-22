#' lsmm : Estimation of a linear mixed model for longitudinal data with flexible subject-specific variability.
#'
#' This function fits linear mixed effects models for longitudinal data, allowing
#' the residual variance to be subject-specific. Three different model types can
#' be estimated (see *Details*).Parameters are estimated by maximum likelihood
#' using a Marquardt-Levenberg algorithm.
#'
#'
#'
#' @details
#'
#'
#' The model is defined as:
#'
#' \eqn{Y_{ij} = Y_{i}(t_{ij}) = \widetilde{Y}_i(t_{ij}) + \epsilon_{ij} = X_{ij}^{\top} \beta+Z_{ij}^{\top} b_{i}+\epsilon_{ij}},
#'
#' where \eqn{X_{ij}} and \eqn{Z_{ij}} are vectors of explanatory variables for subject \eqn{i}
#' at time \eqn{t_{ij}}, associated with the fixed-effect vector \eqn{\beta} and
#' the subject-specific random-effect vector \eqn{b_i}, respectively.
#'
#' **A. Standard linear mixed model**
#'
#' In this case, \eqn{b_i \sim \mathcal{N}(0,B)}, with \eqn{B} an unstructured
#' covariance matrix, and the measurement errors \eqn{\epsilon_{ij}} are independent
#' Gaussian errors with variance \eqn{\sigma^2_{\epsilon}}.
#'
#' **B. Location-scale mixed model with time and/or covariate-dependent variability**
#'
#' In this model, the residual error variance can vary across subjects and over time:
#'  we assume the following specification for the residual error:
#'
#' \eqn{\epsilon_{ij} \sim \mathcal{N}(0,\sigma_i^2)} with \eqn{ \log(\sigma_i(t_{ij}))  = O_{ij}^{\top} \mu+M_{ij}^{\top} \tau_{i}}.
#'
#' where \eqn{O_{ij}} and \eqn{M_{ij}} are vectors of explanatory variables for subject \eqn{i}
#' at visit \eqn{j}, associated with the fixed-effect vector \eqn{\mu} and
#' the subject-specific random-effect vector \eqn{\tau_i}, respectively.
#'
#' The random effects are assumed jointly Gaussian:
#'
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
#' By convention, random effects for the mean (\eqn{b_i}) can be assumed
#' independent from those for the variance (\eqn{\tau_i}) by setting \eqn{\Sigma_{\tau b} =0}.
#'
#' **C. Location-scale mixed model distinguishing within and between visits variabilities**
#'
#' In some studies, multiple measurements of the same marker are collected during each visit.
#' For instance, in clinical research, two or three blood pressure readings are
#' typically recorded per visit, and within-visit variability may carry information.
#'
#' To account for this, we introduce an LSMM that distinguishes within- and
#' between-visit variabilities.
#'
#' We introduce an additional level in the data hierarchy, grouping repeated
#' measures by visit.
#'
#' For each subject \eqn{i} \eqn{(i=1,...,N)}, \eqn{Y_{ijl}} represents the
#' \eqn{l}-th measurement \eqn{(l=1,...,n_{ij})} at visit \eqn{j} \eqn{(j=1,...,n_i)}
#' and time \eqn{t_{ij}}. We then define
#'
#' \eqn{
#' \left\{
#'   \begin{array}{ll}
#'   Y_{ijl} =  \widetilde{Y}_i(t_{ij}) + \epsilon_{ij} + \nu_{ijl} = X_{ij}^{\top} \beta+Z_{ij}^{\top} b_{i}+\epsilon_{ij} + \nu_{ijl}, \\
#'   \epsilon_{ij} \sim \mathcal{N}(0,\sigma_i^2) \hspace{4mm} \text{with} \hspace{3mm} \log(\sigma_i)  = \mu_\sigma + \tau_{\sigma i},\\
#'   \nu_{ijl} \sim \mathcal{N}(0,\kappa_i^2) \hspace{3mm} \text{with} \hspace{3mm} \log(\kappa_i)  = \mu_\kappa + \tau_{\kappa i},\\
#'   \end{array}
#'   \right.}
#'
#' where \eqn{\mu_\sigma} and \eqn{\mu_\kappa} are fixed intercepts for the
#' between-visits and within-visit variances, respectively. The subject-specific
#' random-effect \eqn{b_i} and \eqn{\tau_i = (\tau_{\sigma i},\tau_{\kappa i})^\top}
#' are assumed to be Gaussian as in model (B).
#'
#'
#'
#'
#' @param formFixed Formula specifying the fixed effects of the longitudinal model.
#' @param formRandom Formula specifying the random effects of the longitudinal model.
#' @param formGroup Formula specifying the grouping variable (typically the subject ID).
#' @param timeVar Character string specifying the time variable.
#' @param formVar Character string indicating the type of variability: **"standard"** for a standard LMM,
#' **"cov-dependent"** for a covariate-dependent residual variance, or **"inter-intra"**
#' to distinguish inter- and intra-visit variability.
#' @param formFixedVar Formula specifying the fixed effects for the variance predictor (if \code{formVar == "cov-dependent"}).
#' @param formRandomVar Formula specifying the random effects for the variance predictor (if \code{formVar == "cov-dependent"}).
#' @param random_inter Logical indicating whether the between-visit variability is subject-specific (used when \code{formVar = "inter-intra"}).
#' @param random_intra Logical indicating whether the within-visit variability is subject-specific (used when \code{formVar = "inter-intra"}).
#' @param formGroupVisit Formula specifying the visit indicator variable (used when \code{formVar = "inter-intra"}).
#' @param correlated_re Logical indicating whether the random effects for the mean and variance submodels are correlated (used when \code{formVar} is in \code{c("cov-dependent", "inter-intra")}).
#' @param data.long Data frame containing the longitudinal data.
#' @param S1 Integer specifying the number of QMC draws for the first step.
#' @param S2 Integer specifying the number of QMC draws for the second step.
#' @param nproc Integer specifying the number of processors for parallel computing.
#' @param clustertype Character string indicating the cluster type supported by \code{makeCluster}.
#' @param maxiter Optional integer specifying the maximum number of iterations for the Marquardt–Levenberg algorithm (default to 100).
#' @param print.info Logical indicating whether iteration details should be printed (False by default).
#' @param file Optional character string giving the name of the file where iteration outputs are written (if \code{print.info = TRUE}).
#' @param epsa Optional numeric threshold for convergence based on parameter stability.
#' @param epsb Optional numeric threshold for convergence based on objective function stability.
#' @param epsd  Optional numeric threshold for the relative distance to the maximum. This criterion has the nice interpretation of estimating the ratio of the approximation error over the statistical error, thus it can be used for stopping the iterative process whatever the problem.
#' @param binit Optional vector of initial parameters values.
#'
#' @return An object of class \code{lsmm} containing:
#' \describe{
#' \item{\code{table.res}}{Table of parameter estimates and standard errors.}
#' \item{\code{result_step1}}{A \code{marqLevAlg} object with first-step estimation results.}
#' \item{\code{result_step2}}{A \code{marqLevAlg} object with second-step estimation results.}
#' \item{\code{info_conv_step1}}{Information on first-step convergence (criteria and computation time).}
#' \item{\code{info_conv_step2}}{Information on second-step convergence (criteria and computation time).}
#' \item{\code{control}}{List of control parameters used during estimation.}
#'
#' }
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' data(threeC)
#' threeC$age.visit65 <- (threeC$age.visit-65)/10
#' threeC$SBP <- threeC$SBP/10
#' threeC <- threeC
#' threeC <- dplyr::group_by(threeC, ID, num.visit)
#' threeC <- dplyr::mutate(threeC, SBPvisit = mean(SBP))
#' threeC_ex1 <- threeC[!duplicated(threeC[, c("ID", "num.visit")]),
#'                                  c("ID", "SBPvisit", "age.visit65", "sex")]
#'
#' #First example : a standard linear mixed model (constant residual variance,
#' # in time and between subjects, case A in details)
#'
#' m1 <- lsmm(formFixed = SBPvisit ~ age.visit65+I(age.visit65^2),
#'                formRandom = ~ age.visit65+I(age.visit65^2),
#'                formGroup = ~ ID,
#'                timeVar = 'age.visit65',
#'                data.long = threeC_ex1,
#'                formVar = "standard",
#'                S1 = 500,
#'                S2 = 1000,
#'                nproc = 1)
#'
#' summary(m1)
#'
#'
#' #Second example : a linear mixed model with subject-specific time-dependent
#' # and covariate-dependent variability (case B in details)
#'
#' #We adjust the individual residual variability on age and the sex.
#'
#' m2 <- lsmm(formFixed = SBPvisit ~ age.visit65,
#'                formRandom = ~ age.visit65,
#'                formGroup = ~ ID,
#'                timeVar = 'age.visit65',
#'                data.long = threeC_ex1,
#'                formVar = "cov-dependent",
#'                formFixedVar = ~ age.visit65+sex,
#'                formRandomVar = ~ age.visit65,
#'                correlated_re = FALSE,
#'                S1 = 500,
#'                S2 = 1000,
#'                nproc = 1)
#'
#' summary(m2)
#'
#' #Third example : a linear mixed model with subject-specific inter-visits and
#' # intra-visits variabilities
#'
#' threeC_ex2 <- threeC[, c("ID", "SBP", "age.visit65", "sex", "num.visit")]
#'
#' m3 <- lsmm(formFixed = SBP ~ age.visit65,
#'                formRandom = ~ age.visit65,
#'                formGroup = ~ ID,
#'                timeVar = 'age.visit65',
#'                data.long = threeC_ex2,
#'                formVar = "inter-intra",
#'                random_inter = TRUE,
#'                random_intra = TRUE,
#'                formGroupVisit = ~num.visit,
#'                correlated_re = FALSE,
#'                S1 = 500,
#'                S2 = 1000,
#'                nproc = 1)
#'
#' summary(m3)
#'
#'
#'
#'
#' #Fourth example : a linear mixed model with subject-specific inter-visits
#' # variability and constant intra-visit variability
#'
#' m4 <- lsmm(formFixed = SBP ~ age.visit65+sex,
#'                formRandom = ~ age.visit65,
#'                formGroup = ~ ID,
#'                timeVar = 'age.visit65',
#'                data.long = threeC_ex2,
#'                formVar = "inter-intra",
#'                random_inter = TRUE,
#'                random_intra = FALSE,
#'                formGroupVisit = ~num.visit,
#'                correlated_re = FALSE,
#'                S1 = 500,
#'                S2 = 1000,
#'                nproc = 1)
#'
#' summary(m4)
#'
#' }
#'
#'
#'
#' @importFrom dplyr %>% group_by mutate
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
