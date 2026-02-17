#' lsjm : Estimation of a location-scale joint model for longitudinal data with flexible subject-specific variability and complex survival structures
#'
#' This function fits joint models in which the residual error variance is allowed
#' to be subject-specific and the survival submodel may include either a single
#' event, competing risks, or an illness-death process.
#' Up to nine different model structures can be estimated (see *Details*).
#' Parameters are estimated by maximum likelihood using a Marquardt-Levenberg algorithm.
#'
#' @details
#'
#' A joint model is composed of two submodels:
#' (1) a linear mixed model or location-scale mixed model(see \code{?LSJM::lsmm}), and
#' (2) a survival model.
#'
#' Three types of survival processes are supported and are described below.
#'
#' **A. A single event: proportional hazards model**
#'
#' In this case, we consider only a single event. Let \eqn{T_i = \min(T^*_{i1}, C_i)} denote the observed time,
#' where \eqn{T^*_{i1}} is the true event time and \eqn{C_i} the censoring time for subject \eqn{i}.
#' Let \eqn{\delta_i \in \{0,1\}} be the event indicator.
#' The hazard function is given by:
#' \deqn{
#' \lambda_i^{01}(t|r_i) = \lambda_0(t) \exp\left(
#' W_i^{01\top}\gamma^{01} +
#' g_y^{01}(b_i, t)^\top \alpha_b^{01} +
#' g_\tau^{01}(\tau_i, t)^\top \alpha_\tau^{01}
#' \right)
#' }
#'
#' where:
#' * \eqn{\lambda_0(t)} is the baseline hazard function,
#' * \eqn{W_i^{01}} is a vector of baseline covariates with associated coefficients \eqn{\gamma^{01}},
#' * \eqn{\alpha_b^{01}} are regression coefficients for the function \eqn{g_y}, representing the association between the event risk and the mean trajectory of \eqn{Y},
#' * \eqn{\alpha_\tau^{01}} are regression coefficients for the function \eqn{g_\tau}, representing the association between the event risk and the residual variance.
#'
#' The association function \eqn{g_y(b_i,t)} can be defined as:
#' \itemize{
#'   \item \eqn{g_y(b_i,t) = \tilde{y}_i(t)} — current value,
#'   \item \eqn{g_y(b_i,t) = \tilde{y}'_i(t) = \frac{\partial \tilde{y}_i(t)}{\partial t}} — current slope,
#'   \item \eqn{g_y(b_i,t) = (\tilde{y}_i(t), \tilde{y}'_i(t))} — both value and slope,
#'   \item \eqn{g_y(b_i,t) = b_i} — random effects.
#' }
#'
#' The association function \eqn{g_\tau(\tau_i, t)} is defined according to the longitudinal model type:
#' \itemize{
#'   \item If a standard linear mixed model with homogeneous residual variance is used: \eqn{g_\tau(\tau_i, t) = 0}.
#'   \item If a location–scale mixed model with time- or covariate-dependent variance is used: \eqn{g_\tau(\tau_i, t) = \sigma_i(t)}.
#'   \item If both within- and between-visit variances are modelled: \eqn{g_\tau(\tau_i, t) = (\sigma_i, \kappa_i)^\top},
#'   where \eqn{\alpha_\tau = (\alpha_\sigma, \alpha_\kappa)} correspond to between-visit and within-visit variabilities, respectively.
#' }
#'
#' The baseline hazard function \eqn{\lambda_0(t)} can follow different parametric forms:
#' \itemize{
#'   \item Exponential: \eqn{\lambda_0(t) = \exp(\alpha_0)},
#'   \item Weibull: \eqn{\lambda_0(t) = \zeta^2 t^{\zeta^2 - 1} \exp(\alpha_0)},
#'   \item Gompertz: \eqn{\lambda_0(t) = \kappa_1^2 \exp(\kappa_2 t)},
#'   \item Cubic B-splines with \eqn{Q} knots:
#'   \eqn{\lambda_0(t) = \exp\left( \sum_{q=1}^{Q+4} \eta_q B_q(t, \nu) \right)},
#'   where \eqn{B_q(t, \nu)} denotes the q-th B-spline basis function with knot vector \eqn{\nu}.
#' }
#'
#' **B. Competing events: cause-specific model**
#'
#' When multiple types of events are possible, two competing causes can be modeled.
#' Let \eqn{T_i = \min(T^*_{i1}, T^*_{i2}, C_i)} be the observed time, where \eqn{T^*_{ik}} is the true event time for cause \eqn{k} (with \eqn{k \in \{1, 2\}}),
#' and \eqn{C_i} the censoring time.
#' The event indicator is \eqn{\delta_i \in \{0, 1, 2\}}, where \eqn{\delta_i = k} if cause \eqn{k} occurred and \eqn{0} otherwise.
#' The cause-specific hazard is defined as:
#'
#' \deqn{
#' \lambda_i^{0k}(t) = \lambda_0^{0k}(t)
#' \exp\left(
#' W_i^{0k\top}\gamma^{0k} +
#' g_y^{0k}(b_i,t)^\top \alpha_b^{0k} +
#' g_\tau^{0k}(\tau_i,t)^\top \alpha_\tau^{0k}
#' \right)
#' }
#'
#' The definitions of \eqn{\lambda_0^{0k}}, \eqn{W_i^{0k}}, \eqn{g_y^{0k}}, and \eqn{g_\tau^{0k}} follow the same principles as in the single-event model.

#' **C. Semi-competing events: illness–death model**
#' In this setting, transition to state (1) may be interval-censored, whereas transition to state (2) is observed exactly.
#' The observed event data for subject \eqn{i} are given by:
#' \eqn{D_i = (T_{0i}, L_i, R_i, \delta_i^{(1)}, T_i, \delta_i^{(2)})^\top},
#' where:
#' * \eqn{T_{0i}} — entry time (in case of delayed entry),
#' * \eqn{L_i} — time of the last visit where the subject was in state (0),
#' * \eqn{R_i} — time of the first visit where the subject was observed in state (1),
#' * \eqn{T_i} — minimum between the transition time to state (2) and the censoring time,
#' * \eqn{\delta_i^{(1)} = \mathbb{1}_{R_i < T_i}} — indicator of transition to state (1),
#' * \eqn{\delta_i^{(2)}} — indicator of transition to state (2).
#'
#' Transition intensities from state \eqn{k \in \{0,1\}} to state \eqn{l \in \{1,2\}} follow a proportional hazards model under the Markov assumption:
#'
#' \deqn{
#' \lambda_i^{kl}(t|b_i, \tau_i) =
#' \lambda_0^{kl}(t) \exp\left(
#' W_i^{kl\top}\gamma^{kl} +
#' g_y^{kl}(b_i,t)^\top \alpha_b^{kl} +
#' g_\tau^{kl}(\tau_i,t)^\top \alpha_\tau^{kl}
#' \right)
#' }
#'
#' where the terms are defined analogously to the single-event model.
#'
#' @param Objectlsmm The result from the \code{lsmm()} function.
#' @param survival_type Character string specifying the survival scheme: \code{"Single"}, \code{"CR"} (competing risks), or \code{"IDM"} (illness–death model).
#' @param formSurv_01 One-sided formula specifying covariates for the 0–1 transition.
#' @param formSurv_02 One-sided formula specifying covariates for the 0–2 transition.
#' @param formSurv_12 One-sided formula specifying covariates for the 1-2 transition.
#' @param sharedtype_01 Character vectors defining the dependence structure between the longitudinal and the survival submodel for the 0-1 transition. For standard mixed models, valid options include \code{"value"}, \code{"slope"}, and \code{"random effects"}.
#' For location–scale mixed models, additional terms such as \code{"variability"}, \code{"variability inter"}, or \code{"variability intra"} can be included when subject-specific variances are modelled.
#' @param sharedtype_02 Character vectors defining the dependence structure between the longitudinal and the survival submodel for the 0-2 transition. For standard mixed models, valid options include \code{"value"}, \code{"slope"}, and \code{"random effects"}.
#' For location–scale mixed models, additional terms such as \code{"variability"}, \code{"variability inter"}, or \code{"variability intra"} can be included when subject-specific variances are modelled.
#' @param sharedtype_12 Character vectors defining the dependence structure between the longitudinal and the survival submodel for the 1-2 transition. For standard mixed models, valid options include \code{"value"}, \code{"slope"}, and \code{"random effects"}.
#' For location–scale mixed models, additional terms such as \code{"variability"}, \code{"variability inter"}, or \code{"variability intra"} can be included when subject-specific variances are modelled.
#' @param hazardBase_01 Character strings specifying the baseline hazard function for transition 0-1: one of \code{"Exponential"}, \code{"Weibull"}, \code{"Gompertz"}, or \code{"Splines"}.
#' @param hazardBase_02 Character strings specifying the baseline hazard function for transition 0-2: one of \code{"Exponential"}, \code{"Weibull"}, \code{"Gompertz"}, or \code{"Splines"}.
#' @param hazardBase_12 Character strings specifying the baseline hazard function for transition 1-2: one of \code{"Exponential"}, \code{"Weibull"}, \code{"Gompertz"}, or \code{"Splines"}.
#' @param delta1 One-sided formula defining the event indicator for the first event (\eqn{1} for event, \eqn{0} otherwise).
#' @param delta2 One-sided formula defining the indicator for the second event (\eqn{1} for event, \eqn{0} otherwise).
#' @param Time_T One-sided formula specifying the event time variable.
#' @param Time_L For IDM models, one-sided formula giving the time of the last observation in state (0).
#' @param Time_R For IDM models, one-sided formula giving the first observed time in state (1).
#' @param Time_T0 One-sided formula specifying the time of entry (for delayed entry).
#' @param formSlopeFixed One-sided formula for the time derivative of the fixed effects (if "slope" is included in \code{sharedtype}).
#' @param formSlopeRandom One-sided formula for the time derivative of the random effects (if "slope" is included in \code{sharedtype}).
#' @param index_beta_slope Vector indicating the indices of fixed-effect parameters used in the slope association.
#' @param index_b_slope Vector indicating the indices of random-effect parameters used in the slope association.
#' @param nb.knots.splines Integer giving the number of internal knots for spline-based baseline hazards.
#' @param nb_pointsGK Integer specifying the number of Gauss–Kronrod quadrature points (between 7 and 15; default is 15).
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
#'
#' @return An object of class \code{lsjm} containing:
#' \describe{
#' \item{\code{table.res}}{Table of parameter estimates and standard errors.}
#' \item{\code{result_step1}}{A \code{marqLevAlg} object with first-step estimation results.}
#' \item{\code{result_step2}}{A \code{marqLevAlg} object with second-step estimation results.}
#' \item{\code{info_conv_step1}}{Information on first-step convergence (criteria and computation time).}
#' \item{\code{info_conv_step2}}{Information on second-step convergence (criteria and computation time).}
#' \item{\code{control}}{List of control parameters used during estimation.}
#' }
#'
#' @examples
#'
#' \donttest{
#'
#' # First, run the examples from the lsmm function (see ?lsmm).
#'
#' # Example 1: Illness-death model with time-dependent subject-specific
#' # variability
#' data(threeC)
#' threeC$age.visit65 <- (threeC$age.visit-65)/10
#' threeC$SBP <- threeC$SBP/10
#' threeC <- threeC
#' threeC <- dplyr::group_by(threeC, ID, num.visit)
#' threeC <- dplyr::mutate(threeC, SBPvisit = mean(SBP))
#' threeC_ex1 <- threeC[!duplicated(threeC[, c("ID", "num.visit")]),
#'                                   c("ID", "SBPvisit", "age.visit65", "sex")]
#'
#' m1 <- lsmm(formFixed = SBPvisit ~ age.visit65,
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
#' l1 <- lsjm(m1,
#'            survival_type = 'IDM',
#'            formSurv_01=~1,
#'            formSurv_02=~sex,
#'            formSurv_12=~sex,
#'            sharedtype_01 = c("value", "variability"),
#'            sharedtype_02 = c("value", "slope", "variability"),
#'            sharedtype_12 = c("value", "variability"),
#'            hazardBase_01 = "Weibull",
#'            hazardBase_02 = "Weibull",
#'            hazardBase_12 = "Splines",
#'            delta1=~dem,
#'            delta2=~death,
#'            Time_T =~age.final65,
#'            Time_L =~age.last65,
#'            Time_R =~age.first65,
#'            Time_T0 =~age0_65,
#'            formSlopeFixed =~1,
#'            formSlopeRandom =~1,
#'            index_beta_slope = c(2),
#'            index_b_slope = c(2),
#'            nb.knots.splines = c(0,0,1),
#'            S1 = 1000,
#'            S2 = 2000,
#'            nproc = 10)
#'
#' summary(l1)
#'
#'
#' # Example 2: Competing risks model with between- and within-visit
#' # subject-specific variabilities
#'
#' m2 <- lsmm(formFixed = SBP ~ age.visit65,
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
#'
#' l2 <- lsjm(Objectlsmm = m2,
#'            survival_type = 'CR',
#'            formSurv_01 = ~ sex,
#'            formSurv_02 = ~ sex,
#'            sharedtype_01 = c("value", "variability inter"),
#'            sharedtype_02 = c("value", "variability inter",
#'                              "variability intra"),
#'            hazardBase_01 = "Weibull",
#'            hazardBase_02 = "Weibull",
#'            delta1 = ~ demCR,
#'            delta2 = ~ deathCR,
#'            Time_T0 = ~ age0_65,
#'            Time_T = ~ age65_CR,
#'            nproc = 5,
#'            S1 = 1000,
#'            S2 = 2000)
#'
#' summary(l2)
#' }
#'
#'
#' @import survival
#' @import flexsurv
#' @import marqLevAlg
#' @import splines
#' @importFrom survival Surv
#' @export
#'
#'





lsjm <- function(Objectlsmm,  survival_type = c('Single', 'CR', 'IDM'),
                 formSurv_01, formSurv_02 = NULL, formSurv_12 = NULL,
                 sharedtype_01, sharedtype_02 = NULL, sharedtype_12 = NULL,
                 hazardBase_01, hazardBase_02 = NULL, hazardBase_12 = NULL,
                 delta1, delta2 = NULL, Time_T, Time_L = NULL, Time_R = NULL, Time_T0 = NULL,
                 formSlopeFixed = NULL, formSlopeRandom = NULL, index_beta_slope = NULL, index_b_slope = NULL, nb.knots.splines = c(1,1,1),
                 nb_pointsGK = 15, S1 = 1000, S2 = 5000, binit = NULL,
                 nproc = 1, clustertype = "SOCK", maxiter = 500,
                 print.info = FALSE, file = NULL, epsa = 1e-04, epsb = 1e-04, epsd = 1e-04
){
  if(missing(Objectlsmm)) stop("The argument Objectlsmm must be specified")
  if(!inherits(Objectlsmm, "lsmm_classic") && !inherits(Objectlsmm, "lsmm_covDep") && !inherits(Objectlsmm, "lsmm_interintra")) stop("use only \"lsmm\" objects")


  timeVar <- Objectlsmm$control$timeVar
  formVar <- Objectlsmm$control$formVar

  if(missing(survival_type)) stop("The argument survival_type must be specified")
  if(!inherits((survival_type),"character")) stop("The argument survival_type must be a character")
  if(length(survival_type) != 1) stop("The argument survival_type must be of length 1")
  if(!(survival_type %in% c('Single', 'CR', 'IDM'))) stop("The argument survival_type must be equal to 'Single' when there is only one event or 'CR' for competing risk model or 'IDM' for an illness-death model with possibly interval censoring")
  if(missing(formSurv_01)) stop("The argument formSurv_01 must be specified")
  if(!inherits((formSurv_01),"formula")) stop("The argument formSurv_01 must be a formula")
  if(survival_type %in% c('CR', 'IDM') && missing(formSurv_02)) stop("The argument formSurv_02 must be specified")
  if(survival_type %in% c('CR', 'IDM') && !inherits((formSurv_02),"formula")) stop("The argument formSurv_01 must be a formula")
  if(survival_type %in% c('IDM') && missing(formSurv_12)) stop("The argument formSurv_12 must be specified")
  if(survival_type %in% c('IDM') && !inherits((formSurv_12),"formula")) stop("The argument formSurv_12 must be a formula")

  if(!inherits((delta1),"formula")) stop("The argument delta1 must be a formula")
  if(survival_type %in% c('CR','IDM') && missing(delta2)) stop("The argument delta2 must be specified")
  if(survival_type %in% c('CR','IDM') && !inherits((delta2),"formula")) stop("The argument delta2 must be a formula")

  deltas <- list('delta1' = delta1,
                 'delta2' = delta2)

  if(missing(sharedtype_01)) stop("The argument sharedtype_01 must be specified")
  if(!inherits((sharedtype_01),"character")) stop("The argument sharedtype_01 must be a character")
  if(survival_type %in% c('CR', 'IDM') && missing(sharedtype_02)) stop("The argument sharedtype_02 must be specified")
  if(survival_type %in% c('CR', 'IDM') && !inherits((sharedtype_02),"character")) stop("The argument sharedtype_02 must be a character")
  if(survival_type %in% c('IDM') && missing(sharedtype_12)) stop("The argument sharedtype_12 must be specified")
  if(survival_type %in% c('IDM') && !inherits((sharedtype_12),"character")) stop("The argument sharedtype_12 must be a character")

  if(missing(hazardBase_01)) stop("The argument hazardBase_01 must be specified")
  if(!inherits((hazardBase_01),"character")) stop("The argument hazardBase_01 must be a character")
  if(survival_type %in% c('CR', 'IDM') && missing(hazardBase_02)) stop("The argument hazardBase_02 must be specified")
  if(survival_type %in% c('CR', 'IDM') && !inherits((hazardBase_02),"character")) stop("The argument hazardBase_02 must be a character")
  if(survival_type %in% c('IDM') && missing(hazardBase_12)) stop("The argument hazardBase_12 must be specified")
  if(survival_type %in% c('IDM') && !inherits((hazardBase_12),"character")) stop("The argument hazardBase_12 must be a character")
  if('slope' %in% c(sharedtype_01, sharedtype_02, sharedtype_12) && missing(formSlopeFixed)) stop("The argument formSlopeFixed must be specified")
  if('slope' %in% c(sharedtype_01, sharedtype_02, sharedtype_12) && !inherits((formSlopeFixed),"formula")) stop("The argument formSlopeFixed must be specified")
  if('slope' %in% c(sharedtype_01, sharedtype_02, sharedtype_12) && missing(formSlopeRandom)) stop("The argument formSlopeRandom must be specified")
  if('slope' %in% c(sharedtype_01, sharedtype_02, sharedtype_12) && !inherits((formSlopeRandom),"formula")) stop("The argument formSlopeRandom must be specified")
  if('slope' %in% c(sharedtype_01, sharedtype_02, sharedtype_12) && missing(index_beta_slope)) stop("The argument index_beta_slope must be specified")
  if('slope' %in% c(sharedtype_01, sharedtype_02, sharedtype_12) && missing(index_b_slope)) stop("The argument index_b_slope must be specified")

  Time <- list(Time_L = Time_L,
               Time_R = Time_R,
               Time_T = Time_T,
               Time_T0 = Time_T0)


  if(formVar == 'classic'){
    if(!all(sharedtype_01 %in% c("value", "slope", "random effects"))) stop("The argument of sharedtype_01 must be in c('value', 'slope', 'random effects')")
    if(survival_type %in% c('CR', 'IDM') && !all(sharedtype_02 %in% c("value", "slope", "random effects"))) stop("The argument of sharedtype_02 must be in c('value', 'slope', 'random effects')")
    if(survival_type %in% c('IDM') && !all(sharedtype_12 %in% c("value", "slope", "random effects"))) stop("The argument of sharedtype_12 must be in c('value', 'slope', 'random effects')")
    if(survival_type == 'Single'){
      result <- lsjm_classicSingle(Objectlsmm, Time, deltas, hazardBase_01,  nb.knots.splines,
                                   formSurv_01,   nb_pointsGK, sharedtype_01, formSlopeFixed, formSlopeRandom,
                                   index_beta_slope , index_b_slope, timeVar,
                                   S1, S2, binit, nproc , clustertype, maxiter,
                                   print.info, file , epsa , epsb , epsd )
    }
    else{
      if(survival_type == 'CR'){
        result <- lsjm_classicCR(Objectlsmm, Time, deltas, hazardBase_01, hazardBase_02, nb.knots.splines,
                                 formSurv_01, formSurv_02,  nb_pointsGK, sharedtype_01, sharedtype_02,
                                 formSlopeFixed, formSlopeRandom,
                                 index_beta_slope , index_b_slope, timeVar,
                                 S1, S2, binit, nproc , clustertype, maxiter,
                                 print.info, file , epsa , epsb , epsd)
      }
      else{
        result <- lsjm_classicIDM(Objectlsmm, Time, deltas, hazardBase_01, hazardBase_02, hazardBase_12, nb.knots.splines,
                                  formSurv_01, formSurv_02, formSurv_12, nb_pointsGK, sharedtype_01, sharedtype_02, sharedtype_12,
                                  formSlopeFixed, formSlopeRandom,
                                  index_beta_slope , index_b_slope, timeVar,
                                  S1, S2, binit, nproc , clustertype, maxiter,
                                  print.info, file , epsa , epsb , epsd)
      }
    }
  }
  else{
    if(formVar == 'cov-dependent'){
      if(!all(sharedtype_01 %in% c("value", "slope", "variability", "random effects"))) stop("The argument of sharedtype_01 must be in c('value', 'slope', 'variability', 'random effects')")
      if(survival_type %in% c('CR', 'IDM') && !all(sharedtype_02 %in% c("value", "slope", "variability", "random effects"))) stop("The argument of sharedtype_02 must be in c('value', 'slope', 'variability', 'random effects')")
      if(survival_type %in% c('IDM') && !all(sharedtype_12 %in% c("value", "slope", "variability", "random effects"))) stop("The argument of sharedtype_12 must be in c('value', 'slope', 'variability', 'random effects')")
      if(survival_type == 'Single'){
        result <- lsjm_covDepSingle(Objectlsmm, Time, deltas, hazardBase_01,  nb.knots.splines,
                                    formSurv_01,  nb_pointsGK, sharedtype_01,
                                    formSlopeFixed, formSlopeRandom,
                                    index_beta_slope , index_b_slope, timeVar,
                                    S1, S2, binit, nproc , clustertype, maxiter,
                                    print.info, file , epsa , epsb , epsd)
      }
      else{
        if(survival_type == 'CR'){
          result <- lsjm_covDepCR(Objectlsmm, Time, deltas, hazardBase_01, hazardBase_02, nb.knots.splines,
                                  formSurv_01, formSurv_02,  nb_pointsGK, sharedtype_01, sharedtype_02,
                                  formSlopeFixed, formSlopeRandom,
                                  index_beta_slope , index_b_slope,timeVar,
                                  S1, S2, binit, nproc , clustertype, maxiter,
                                  print.info, file , epsa , epsb , epsd)
        }
        else{
          result <- lsjm_covDepIDM(Objectlsmm, Time, deltas, hazardBase_01, hazardBase_02, hazardBase_12, nb.knots.splines,
                                   formSurv_01, formSurv_02, formSurv_12, nb_pointsGK, sharedtype_01, sharedtype_02, sharedtype_12,
                                   formSlopeFixed, formSlopeRandom,
                                   index_beta_slope , index_b_slope,timeVar,
                                   S1, S2, binit, nproc , clustertype, maxiter,
                                   print.info, file , epsa , epsb , epsd)
        }
      }
    }
    else{
      if(Objectlsmm$control$var_inter && Objectlsmm$control$var_intra){
        if(!all(sharedtype_01 %in% c("value", "slope", "variability inter", "variability intra", "random effects"))) stop("The argument of sharedtype_01 must be in c('value', 'slope', 'variability inter', 'variability intra', 'random effects')")
        if(survival_type %in% c('CR', 'IDM') && !all(sharedtype_02 %in% c("value", "slope", "variability inter", "variability intra", "random effects"))) stop("The argument of sharedtype_02 must be in c('value', 'slope', 'variability inter', 'variability intra', 'random effects')")
        if(survival_type %in% c('IDM') && !all(sharedtype_12 %in% c("value", "slope", "variability inter", "variability intra", "random effects"))) stop("The argument of sharedtype_12 must be in c('value', 'slope', 'variability inter', 'variability intra', 'random effects')")
      }
      else{
        if(Objectlsmm$control$var_inter){
          if(!all(sharedtype_01 %in% c("value", "slope", "variability inter", "random effects"))) stop("The argument of sharedtype_01 must be in c('value', 'slope', 'variability inter', 'random effects')")
          if(survival_type %in% c('CR', 'IDM') && !all(sharedtype_02 %in% c("value", "slope", "variability inter", "random effects"))) stop("The argument of sharedtype_02 must be in c('value', 'slope', 'variability inter', 'random effects')")
          if(survival_type %in% c('IDM') && !all(sharedtype_12 %in% c("value", "slope", "variability inter", "random effects"))) stop("The argument of sharedtype_12 must be in c('value', 'slope', 'variability inter', 'random effects')")
        }
        else{
          if(Objectlsmm$control$var_intra){
            if(!all(sharedtype_01 %in% c("value", "slope", "variability intra", "random effects"))) stop("The argument of sharedtype_01 must be in c('value', 'slope', 'variability intra', 'random effects')")
            if(survival_type %in% c('CR', 'IDM') && !all(sharedtype_02 %in% c("value", "slope", "variability intra", "random effects"))) stop("The argument of sharedtype_02 must be in c('value', 'slope', 'variability intra', 'random effects')")
            if(survival_type %in% c('IDM') && !all(sharedtype_12 %in% c("value", "slope", "variability intra", "random effects"))) stop("The argument of sharedtype_12 must be in c('value', 'slope', 'variability intra', 'random effects')")
          }
          else{
            if(Objectlsmm$control$var_intra){
              if(!all(sharedtype_01 %in% c("value", "slope", "random effects"))) stop("The argument of sharedtype_01 must be in c('value', 'slope', 'random effects')")
              if(survival_type %in% c('CR', 'IDM') && !all(sharedtype_02 %in% c("value", "slope", "random effects"))) stop("The argument of sharedtype_02 must be in c('value', 'slope', 'random effects')")
              if(survival_type %in% c('IDM') && !all(sharedtype_12 %in% c("value", "slope", "random effects"))) stop("The argument of sharedtype_12 must be in c('value', 'slope', 'random effects')")
            }
          }
        }
      }
      if(survival_type == 'Single'){
        result <- lsjm_interintraSingle(Objectlsmm, Time, deltas, hazardBase_01, nb.knots.splines,
                                        formSurv_01, nb_pointsGK, sharedtype_01,
                                        formSlopeFixed, formSlopeRandom,
                                        index_beta_slope , index_b_slope,timeVar,
                                        S1, S2, binit,nproc , clustertype, maxiter,
                                        print.info, file , epsa , epsb , epsd)
      }
      else{
        if(survival_type == 'CR'){
          result <- lsjm_interintraCR(Objectlsmm, Time, deltas, hazardBase_01, hazardBase_02, nb.knots.splines,
                                      formSurv_01, formSurv_02,  nb_pointsGK, sharedtype_01, sharedtype_02,
                                      formSlopeFixed, formSlopeRandom,
                                      index_beta_slope , index_b_slope,timeVar,
                                      S1, S2, binit,nproc , clustertype, maxiter,
                                      print.info, file , epsa , epsb , epsd)
        }
        else{
          result <- lsjm_interintraIDM(Objectlsmm, Time, deltas, hazardBase_01, hazardBase_02, hazardBase_12, nb.knots.splines,
                                       formSurv_01, formSurv_02, formSurv_12, nb_pointsGK, sharedtype_01, sharedtype_02, sharedtype_12,
                                       formSlopeFixed, formSlopeRandom,
                                       index_beta_slope , index_b_slope,timeVar,
                                       S1, S2, binit,nproc , clustertype, maxiter,
                                       print.info, file , epsa , epsb , epsd)
        }
      }

    }

  }
  result
}
