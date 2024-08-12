data.time <- function(data.id, Time, formFixed, formRandom, timeVar){
  if (!timeVar %in% names(data.id))
    stop("\n'timeVar' does not correspond to one of the columns in formulas")
  data.id[[timeVar]] <- Time
  mfX.id <- model.frame(formFixed, data = data.id)
  Xtime <- model.matrix(formFixed, mfX.id)
  mfU.id <- model.frame(formRandom, data = data.id)
  Utime <- model.matrix(formRandom, mfU.id)

  list("Xtime" = Xtime, "Utime" = Utime)
}
