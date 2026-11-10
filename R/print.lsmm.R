#' @export

print.lsmm <- function(object, ...){

  model_name <- switch(
    class(object)[1],
    "lsmm_classic" = "Linear mixed-effect model",
    "lsmm_covDep" = "Location-scale linear mixed model",
    "lsmm_interintra" = "Location-scale linear mixed model",
    "lsjm_classicSingle" = "Joint model for a single event",
    "lsjm_classicCR" = "Joint model for competing risks",
    "lsjm_classicIDM" = "Joint model for an illness-death model",
    "lsjm_covDepSingle" = "Location-scale joint model for a single event",
    "lsjm_covDepCR" = "Location-scale joint model for competing risks",
    "lsjm_covDepIDM" = "Location-scale joint model for an illness-death model",
    "lsjm_interintraSingle" = "Location-scale joint model for a single event",
    "lsjm_interintraCR" = "Location-scale joint model for competing risks",
    "lsjm_interintraIDM" = "Location-scale joint model for an illness-death model",
    "Unknown LSJM model"
    )

  cat("\n", model_name, "\n", sep = "")
  cat(strrep("-", nchar(model_name)), "\n", sep = "")

  cat("Number of subjects", object$control$Objectlsmm$control$Ind, "\n")
  cat("Number of parameters", length(object$result_step1$b), "\n")
  cat("Log-likelihood", round(object$result_step1$fn.value,2),"\n")

  invisible(object)
}
