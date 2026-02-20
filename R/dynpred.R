#' dynpred: compute the dynamic predictions of an event.
#'
#' @description
#'
#' This function can be used only with an \code{lsjm} object.
#' Individual dynamic predictions can be computed and plotted.
#' They are defined as the predicted probability of having event \eqn{k} between time \eqn{s} and \eqn{s+t}
#' given that the subject \eqn{i} has not experienced any event before time \eqn{s},
#' and knowing all marker measures collected until time \eqn{s},
#' denoted by \eqn{Y_i(s)}, and the set of estimated parameters. The prediction is defined for subject \eqn{i} by:
#' \deqn{
#' \pi_i^{0k}(s,t;\widehat{\theta}) = P(s<T_i<s+t, \delta_i = k \mid T_i > s, \mathcal{Y}_i(s), \widehat{\theta})
#' = \frac{\int \left[\int_s^{s+t} \exp{(-\sum_{c=1}^K \Lambda_{i}^{0c}(u|r_i,\widehat{\theta}))}\lambda_{i}^{0k}(u|r_i,\widehat{\theta})du \right]f(\mathcal{Y}_i(s)|r_i,\widehat{\theta})f(r_i|\widehat{\theta})dr_i}
#' {\int\exp{(-\sum_{c=1}^K \Lambda_{i}^{0c}(s|r_i,\widehat{\theta}))}f(\mathcal{Y}_i(s)|r_i,\widehat{\theta})f(r_i|\widehat{\theta})dr_i}
#' }
#'
#' For illness-death and competing risks models, \eqn{K = 2}, whereas for a model with a single event, \eqn{K = 1}.
#' For illness-death and single-event models, only the prediction from state (0) to state (1) can be computed,
#' while for competing risks, transitions (0→1) and (0→2) can both be predicted.
#'
#' As in the estimation procedure, the integral over the random effects is computed by QMC approximation and the integral over time by the Gauss-Kronrod quadrature.
#'
#' The 95\% confidence interval of predictions is obtained by the following Monte Carlo algorithm, which can be computationally intensive. For \eqn{L} large enough and \eqn{l = 1,...,L} (e.g. \eqn{L = 1000}):
#' \itemize{
#'   \item Generate \eqn{\widetilde{\theta}^{(l)} \sim \mathcal{N}(\widehat{\theta},V(\widehat{\theta}))}, where \eqn{V(\widehat{\theta})} is the inverse of the Hessian matrix at \eqn{\widehat{\theta}};
#'   \item Compute  \eqn{\widetilde{\pi}^{(l)}_i(s,t;\widetilde{\theta}^{(l)})} from the equation above;
#'   \item Compute the 95\% confidence interval from the 2.5th and 97.5th percentiles of the  \eqn{L}-sample of  \eqn{\widetilde{\pi}^{(l)}_i(s,t;\widetilde{\theta}^{(l)})}.
#' }
#'
#' @param object An \code{lsjm} object.
#' @param newdata A dataset containing the individuals for whom predictions are to be computed.
#' @param s The landmark time (a single numeric value).
#' @param horizon The horizon time of prediction. The function computes the probability of experiencing the event between \code{s} and \code{horizon}. This can be a vector or a single value.
#' @param event An integer indicating for which event the prediction is computed. In the case of a competing risks model, it can be 1 or 2, whereas for other models it can only be 1.
#' @param CI An integer between 0 and 100 indicating the confidence level (default: \code{95} for a 95% CI). If \code{CI = NULL}, no confidence interval is computed.
#' @param nb.draws An integer giving the number of Monte Carlo draws used to compute the confidence interval (default: \code{1000}).
#'
#' @return A dataframe containing the table of predictions (\code{table.pred}).
#'
#' @examples
#'
#' \donttest{
#'
#' # Begin by running the examples from the lsjm function (see ?lsjm).
#'
#' # Example 1: prediction of the risk to be diagnosed with dementia between 85
#' # and 95 years old given the blood pressure measurements before 85 years old
#' # for individual "10003", with an IDM model s = (85-65)/10
#'
#' data(threeC)
#' threeC$age.visit65 <- (threeC$age.visit-65)/10
#' threeC$SBP <- threeC$SBP/10
#' threeC <- threeC
#' threeC <- dplyr::group_by(threeC, ID, num.visit)
#' threeC <- dplyr::mutate(threeC, SBPvisit = mean(SBP))
#' threeC_ex1 <- threeC[!duplicated(threeC[, c("ID", "num.visit")]),
#'                                   c("ID", "SBPvisit", "age.visit65", "dem","death", "sex","age.first65",
#'                                   "age.last65","age.final65","age0_65")]
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
#' # Example 2: prediction of the risk to be diagnosed with dementia
#' # (dynpDementia) and to die (dynpDeath) between 85 and 95 years old
#' # given the blood pressure measurements before 85 years old for
#' # individual "10003", with a competing risk model
#' ind2 <- threeC_ex2[which(threeC_ex2$ID == 10003),]
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
#'
#' dynpDementia <- dynpred(l2, ind2,  s = 2, horizon = seq(2.1,3,0.1),
#'                   event = 1, nb.draws = 1000)
#'
#' dynpDeath <- dynpred(l2, ind2,  s = 2, horizon = seq(2.1,3,0.1),
#'                   event = 2, nb.draws = 1000)
#'}
#'
#'
#' @rdname dynpred
#' @export
dynpred <- function(object, newdata, s, horizon, event, CI = 95,
                    nb.draws = 1000) {
  UseMethod("dynpred")
}
