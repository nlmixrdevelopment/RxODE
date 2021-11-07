.thetamodelVars <- rex::rex(or("tv", "t", "pop", "POP", "Pop", "TV", "T", "cov", "err", "eff"))
.thetaModelReg <- rex::rex(or(
  group(start, .thetamodelVars),
  group(.thetamodelVars, end)))

.etaParts <- c(
  "eta", "ETA", "Eta", "ppv", "PPV", "Ppv", "iiv", "Iiv", "bsv", "Bsv", "BSV",
  "bpv", "Bpv", "BPV", "psv", "PSV", "Psv")

.etaModelReg <- rex::rex(or(group(start, or(.etaParts)), group(or(.etaParts), end)))


#' @export
#' @rdname model
model.function <- function(x, ..., envir=parent.frame()) {
  .ui <- RxODE(x)
  model(x=.ui, ..., envir=envir)
}

#' @export
#' @rdname model
model.rxUi <- function(x, ..., envir=parent.frame()) {
  .ret <- .copyUi(x) # copy so (as expected) old UI isn't affected by the call
}
