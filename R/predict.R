#' predict: Prediction of some quantities for each subjects
#'
#' This function computes different predicted quantities. For a \code{lsmm} object, it is possible to compute, for each subject, their random effects, their value of the marker and the residual variability (for each measurement time)).
#' For a \code{lsjm} object we could also compute cumulative hazard function for each risk transition.
#'
#'
#' @param object Either a \code{lsmm} object or a \code{lsjm} object
#' @param which A vector of characters to indicate which predictions are computed. \code{"RE"} corresponds to the random effects, \code{"Y"} to the marker and \code{"Cum"} to the cumulative risk function(s) (only in the case of a \code{lsjm} object)
#' @param data.long A dataframe containing the longitudinal data for making predictions.
#' @param Objectranef Optional \code{ranef} object containing the predicted random effects for each individual. The Default is NULL and the predict functions should compute the random effects.
#'
#' @return A table for each type of prediction (RE/Y/Cum)
#' \describe{
#' \item{\code{predictRE}}{A table with predicted random effects for each subjects.}
#' \item{\code{predictY}}{A table with predicted marker and residual variability for each subjects and at each measurement time.}
#' \item{\code{predictCum_01}}{A table for predicted cumulative hazard function for transition 0-1.}
#' \item{\code{predictCum_02}}{A table for predicted cumulative hazard function for transition 0-1.}
#' \item{\code{predictCum_12}}{A table for predicted cumulative hazard function for transition 0-1.}
#' \item{\code{grid.time.Cum}}{A vector of times for which the cumulative hazard functions are computed.}
#' }
#'
#' @examples
#' \donttest{
#'
#' data(threeC)
#' threeC$age.visit65 <- (threeC$age.visit-65)/10
#' threeC$SBP <- threeC$SBP/10
#' threeC <- threeC
#' threeC <- dplyr::group_by(threeC, ID, num.visit)
#' threeC <- dplyr::mutate(threeC, SBPvisit = mean(SBP))
#' threeC_ex1 <- threeC[!duplicated(threeC[, c("ID", "num.visit")]), c("ID", "SBPvisit", "age.visit65", "sex")]
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
#' pred.m1 <- predict(m1, which = c("RE","Y"), data.long = threeC_ex1)
#'
#' pred.l1 <- predict(l1, which = c("RE","Y","Cum"), data.long = threeC_ex1)
#'
#'
#' }
#'
#'
#'
#'
#' @name predict
#' @rdname predict
#'
NULL
