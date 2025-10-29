#' The population-average probability of being free of any events at a given time $t$ and covariates
#'
#'
#' @param object A lsjm object
#' @param individual Subject profile for which we want to compute the probability
#' @param time a numeric, the time $t$ at which we want to compute the probability
#' @return A numeric, the population-average probability of being free of any events at a given time $t$ and chosen covariates
#'
#' @examples
#' \donttest{#'
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
#'
#' #l1 <- ... *to be completed*
#'
#' We compute the average-probability of being alive and without dementia at age 85 for a Women
#' # individual <- threeC_ex1[1, c("Sex")]
#' # surv_marg(l1, individual, time = 2)
#'
#'
#' }
#'
#'
#'
#'
#' @name surv_marg
#' @rdname surv_marg
#'
NULL
