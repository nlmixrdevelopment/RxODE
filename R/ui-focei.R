#' Get the THETA/ETA lines from RxODE UI
#'
#' @param rxui This is the RxODE ui object
#' @return The theta/eta lines
#' @author Matthew L. Fidler
#' @noRd
.uiGetThetaEta <- function(rxui) {
  .iniDf <- rxui$iniDf
  .w <- which(!is.na(.iniDf$ntheta))
  .thetas <- lapply(.w, function(i) {
    eval(parse(text=paste0("quote(THETA[", .iniDf$ntheta[i],"] <- ", .iniDf$name[i], ")")))
  })
  .etas <- NULL
  .i2 <- .iniDf[-.w, ]
  if (length(.i2$name) > 0) {
    .i2 <- .i2[.i2$neta1 == .i2$neta2, ]
    .etas <- lapply(seq_along(.i2$name), function(i) {
      eval(parse(text=paste0("quote(ETA[", .i2$neta1[i],"] <- ", .i2$name[i], ")")))
    })
  }
  c(.thetas, .etas)
}
