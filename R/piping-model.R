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
  .modelLines <- .quoteCallInfoLines(match.call(expand.dots = TRUE)[-(1:2)])
}

#' Find the line of the original expression
#'
#' @param expr Expression
#' @param rxui RxODE UI
#' @param errorLine When TRUE,Should the error line be considered, or consider
#'   the non-error line
#' @return The return can be:
#'
#'  - Line number of from the original model, if the model matches one
#'  line uniquely.
#'
#'  - Negative number; This indicates the model property (ie
#'  `f(depot)`) is not found in the model, but the differential
#'  equation `d/dt(depot)` is defined uniquely i the model and is
#'  found at the -(line number) returned.
#'
#'  - `NA` which means the line (or related line) was not found in the
#'  model
#'
#'  - `NULL` which means a line is duplicated in the model
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.getModelLineFromExpression <- function(expr, rxui, errorLine=FALSE) {
  .origLines <- rxui$lstExpr
  .errLines <- rxui$predDf$line
  .ret <- NA_integer_
  for (.i in seq_along(.origLines)) {
    .isErrorLine <- .i %in% .errLines
    if ((errorLine && .isErrorLine) ||
          (!errorLine && !.isErrorLine)) {
      .expr <- .origLines[[.i]]
      if (identical(.expr[[2]], expr)) {
        if (is.na(.ret)) {
          .ret <- .i
        } else {
          return(NULL)
        }
      }
    }
  }
  if (!is.na(.ret)) {
    return(.ret)
  }
  if (!errorLine && length(expr) == 2L) {
    .expr1 <- expr[[1]]
    .expr2 <- expr[[2]]
    if (identical(.expr1, quote(`f`)) ||
          identical(.expr1, quote(`F`)) ||
          identical(.expr1, quote(`alag`)) ||
          identical(.expr1, quote(`lag`)) ||
          identical(.expr1, quote(`rate`)) ||
          identical(.expr1, quote(`dur`))) {
      .expr3 <- eval(parse(text=paste0("quote(d/dt(",as.character(.expr2),"))")))
      for (.i in seq_along(.origLines)) {
        .expr <- .origLines[[.i]]
        if (identical(.expr[[2]], .expr3)) {
          return(-.i)
        }
      }
    }
  }
  return(NA_integer_)
}

