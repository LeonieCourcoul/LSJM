#' ThreeC dataset
#'
#' A random subsample of 500 subjects from the French Three-Cities cohort,
#' aimed at assessing the relation between vascular factors and dementia
#' in the elderly, for the example of the LSJM package.
#'
#' The sample includes participants without dementia at baseline and aged
#' 65 years old or older. Repeated measures of systolic blood pressure were
#' collected over a maximum period of 20 years. At each visit, systolic
#' blood pressure was measured two or three times.
#'
#' @format A data frame with 5248 rows and 13 variables:
#' \describe{
#'   \item{ID}{id of each subject}
#'   \item{Apoe4}{indicator of the apoe4 gene}
#'   \item{Edu}{Educational level (0 = less than 10 years, 1 = more than 10 years)}
#'   \item{SBP}{systolic blood pressure in mmHg}
#'   \item{age.visit}{age at measurement of SBP}
#'   \item{age.final}{age at death or censoring (last news)}
#'   \item{age0}{age at entry in the cohort}
#'   \item{age.last}{age of last visit without dementia}
#'   \item{age.first}{age at visit diagnosis for dementia}
#'   \item{sex}{sex of the subject}
#'   \item{dem}{indicator of dementia}
#'   \item{death}{indicator of death}
#'   \item{num.visit}{visit identifier}
#'   \item{age.visit65}{(age.visit-65)/10}
#'   \item{age.final65}{(age.final-65)/10}
#'   \item{age.last65}{(age.last-65)/10}
#'   \item{age.first65}{(age.first-65)/10}
#'   \item{age0_65}{(age0-65)/10}
#' }
#'
#' @docType data
#' @usage data(threeC)
#' @keywords datasets
#' @name threeC
"threeC"
