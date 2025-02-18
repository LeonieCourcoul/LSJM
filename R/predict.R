#' Prediction of some quantities for each subjects
#'
#' This function computes different predicted quantities. For a \code{lsmm} object, it is possible to compute, for each subject, their random effects, their value of the marker and the residual variability (for each measurement time)).
#' For a \code{lsjm} object we could also compute cumulative hazard function for each risk transition.
#'
#'
#' @param object Either a lsmm object or a lsjm object
#' @param which A vector of characters to indicate which predictions are computed. "RE" corresponds to the random effects, "Y" to the marker and "Cum" to the cumulative risk function(s) (only in the case of a lsjm object)
#' @param data.long A dataframe which contains the longitudinal data for making predictions.
#'
#' @return A table for each type of prediction (RE/Y/Cum)
#' \describe{
#' \item{\code{predictRE}}{A table with predicted random effects for each subjects}
#' \item{\code{predictY}}{A table with predicted marker and residual variability for each subjects and at each measurement time}
#' \item{\code{predictCum...}}{to be completed}
#' }
#'
#' @examples
#' \donttest{#'
#' data(threeC)
#' threeC$age.visit65 <- (threeC$age.visit-65)/10
#' threeC$SBP <- threeC$SBP/10
#' threeC <- threeC
#' threeC <- threeC %>% group_by(ID, num.visit) %>% mutate(SBPvisit = mean(SBP))
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
#' pred.m1 <- predict(m1, which = c("RE","Y"), data.long = threeC_ex1)
#'
#' head(pred.m1$predictRE)
#' head(pred.m1$predictY)
#'
#' l1 <- ... *to be completed*
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
