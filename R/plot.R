#' plot: Plot methods for \code{lsjm} or \code{lsmm} objects
#'
#' This function produces different plots (longitudinal and survival goodness-of-fit,
#' individual trajectory) of a fitted object of class \code{lsmm} or \code{lsjm}.
#'
#' With \code{which="long.fit"}, this function allows to assess the fit of the longitudinal submodel comparing the mean of marker predictions collected in some windows of times (defined by \code{break.times} or using percentiles) to the mean of the observed measurements and its 95\% confidence interval.
#'
#' With \code{which="traj.ind"}, represents the individual trajectory with its prediction interval.
#'
#' For a \code{lsjm} object only, with \code{which="link"}, the function
#' allows to assess the fit of the survival submodel. For each transition, the predicted cumulative hazard function at each event time is computed given the predicted random effects. Then the mean of the predicted cumulative hazard functions are compared with there Nelson-Aalen estimator in the case of a single event or with competing risks. For an illness-death model, the predicted cumulative hazard function for each transition is compared to an illness-death model estimated by penalized likelihood accounting for interval censoring with the \code{SmoothHazard} package
#'
#' @param x A \code{lsmm} or a \code{lsjm} object
#' @param which A character indicating which plot must be display,
#' @param Objectpredict An object of the \code{predict} function
#' @param break.times A vector of breaking times to create windows if \code{which = 'long.fit'}
#' @param ID.ind A vector providing the id of subjects for whom we we wish to trace the individual trajectory.
#' @param ObjectSmoothHazard A SmoothHazard object (only for IDM survival model)
#' @param xlim the x limits (x1, x2) of the plot. The default value, NULL, indicates that the range of the finite values to be plotted should be used.
#' @param ylim the y limits (y1, y2) of the plot. The default value, NULL, indicates that the range of the finite values to be plotted should be used.
#' @param ... Further arguments passed to [graphics::plot()] or other methods.
#'
#' @examples
#'
#' \dontrun{
#' data(threeC)
#' threeC$age.visit65 <- (threeC$age.visit-65)/10
#' threeC$SBP <- threeC$SBP/10
#' threeC <- threeC %>% group_by(ID, num.visit) %>% mutate(SBPvisit = mean(SBP))
#' threeC$age65_CR <- NA
#' threeC$age65_CR[which(threeC$dem == 1)] <- (threeC$age.last65[which(threeC$dem == 1)] +threeC$age.first65[which(threeC$dem == 1)])/2
#' threeC$age65_CR[which(threeC$dem == 0)] <- threeC$age.final65[which(threeC$dem == 0)]
#' threeC$demCR <- threeC$dem
#' threeC$deathCR <- NA
#' threeC$deathCR[which(threeC$dem == 1)] <- 0
#' threeC$deathCR[which(threeC$dem == 0)] <- threeC$death[which(threeC$dem == 0)]
#' threeC <- threeC %>% group_by(ID) %>% filter(age.visit65 <= age65_CR)
#'
#'
#' threeC_ex1 <- threeC[!duplicated(threeC[, c("ID", "num.visit")]), c("ID", "SBPvisit", "age.visit65", "sex","age0_65", "demCR","deathCR","age65_CR")]
#'
#'
#' # Estimation of a lsmm with time-dependent variability model:
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
#' # Predictions:
#' pred.m2 <- predict(m2, which = c("RE", "Y"))
#'
#' #Plot:
#' plot(m2, which = "traj.fit", predictObject = pred.m2,break.times = (seq(65,95,by = 2.5)-65)/10)
#' plot(m2, which = "traj.ind", predictObject = pred.m2,  ID.ind = c(10003,120010))
#'
#' # Estimation of a lsjm with time-dependent variability model and competing events:
#' l2 <- lsjm(Objectlsmm = m2,
#'            survival_type = 'CR',
#'            formSurv_01 = ~ sex,
#'            formSurv_02 = ~ sex,
#'            sharedtype_01 = c("value", "variability"),
#'            sharedtype_02 = c("value", "variability"),
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
#' # Predictions:
#' pred.l2 <- predict(l2, which = c("RE", "Y", "Cum"))
#'
#' # Plot:
#' plot(l2, which = "traj.fit", predictObject = pred.l2,break.times = (seq(65,95,by = 2.5)-65)/10)
#' plot(l2, which = "traj.ind", predictObject = pred.l2, ID.ind = c(10003,120010))
#' plot(l2, which = "survival.fit", predictObject = pred.l2)
#' }
#'
#'
#' @name plot.lsjm
#' @rdname plot.lsjm
#'
#'
#'
NULL



