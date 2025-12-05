#' ranef: the predicted random effects for all subjects
#'
#' @param object Either a lsmm object or a lsjm object
#' @param ... Further arguments
#'
#' @examples
#'
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
#' r1 <- ranef(m1)
#' r1.1 <- ranef(l1)
#'
#'}
#'
#' @name ranef
#' @rdname ranef
#'
#' @export
#'
#'
NULL
