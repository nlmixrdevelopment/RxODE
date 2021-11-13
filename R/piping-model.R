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

#' This gives a equivalent left handed expression
#'
#' @param expr This exchanges expressions ie f() to F()
#' @return NULL if there isn't an equaivlent expression, or the
#'   equivalent expression.
#' @author Matthew L. Fidler
#' @noRd
.getModelLineEquivalentLhsExpression <- function(expr) {
  .expr3 <- NULL
  if (length(expr) == 2L) {
    .expr1 <- expr[[1]]
    .expr2 <- expr[[2]]
    if (identical(.expr1, quote(`f`))) {
      .expr3 <- eval(parse(text=paste0("quote(F(",as.character(.expr2),"))")))
    }
    if (identical(.expr1, quote(`F`))) {
      .expr3 <- eval(parse(text=paste0("quote(f(",as.character(.expr2),"))")))
    }
    if (identical(.expr1, quote(`lag`))) {
      .expr3 <- eval(parse(text=paste0("quote(alag(",as.character(.expr2),"))")))
    }
    if (identical(.expr1, quote(`alag`))) {
      .expr3 <- eval(parse(text=paste0("quote(lag(",as.character(.expr2),"))")))
    }
  }
 .expr3
}

#' Get the model line number from the expression
#'
#' @param expr The Left handed expression to find
#' @param altExpr The alternate expression to find (or `NULL`), this
#'   is for expressions that are equivalent like `f(center)` and
#'   `F(center)`.
#' @param useErrorLine This is a boolean indicating if error lines
#'   should be considered (`TRUE`).
#' @param errLines This is an integer vector of the lines where an
#'   error is defined in the model.
#' @param origLines This is a list of lines in the `model({})` block
#'   of the equation.
#' @return `NULL` for duplicated lines, otherwise the line number (if
#'   it is found) and `NA` if it is not found.
#' @author Matthew L. Fidler
#' @noRd
.getModelineFromExperssionsAndOriginalLines <- function(expr, altExpr, useErrorLine,
                                                        errLines, origLines) {
  .ret <- NA_integer_
  for (.i in seq_along(origLines)) {
    .isErrorLine <- .i %in% errLines
    if ((useErrorLine && .isErrorLine) ||
          (!useErrorLine && !.isErrorLine)) {
      .expr <- origLines[[.i]]
      if (identical(.expr[[2]], expr)) {
        if (is.na(.ret)) {
          .ret <- .i
        } else {
          return(NULL)
        }
      } else if (!is.null(altExpr)) {
        if (identical(.expr[[2]], altExpr)) {
          if (is.na(.ret)) {
            .ret <- .i
          } else {
            return(NULL)
          }
        }
      }
    }
  }
  .ret
}
#' Get the negative model line for ode based on ode property
#'
#' @param expr Left handed expression
#' @param origLines List of original lines (expressions)
#' @param errorLine Consider error lines
#' @return `NA` if the ODE isn't found. Negative number for the lowest
#'   line of the ODE compartment.
#' @author Matthew L. Fidler
#' @noRd
.getNegativeModelLineForDiffFromProperty <- function(expr, origLines, errorLine) {
  .ret <- NA_integer_
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
      for (.i in seq_along(origLines)) {
        .expr <- origLines[[.i]]
        if (identical(.expr[[2]], .expr3)) {
          if (is.na(.ret)) {
            .ret <- -.i
          } else {
            .ret <- min(-.i, .ret)
          }
        }
      }
    }
  }
  .ret
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
.getModelLineFromExpression <- function(lhsExpr, rxui, errorLine=FALSE) {
  .origLines <- rxui$lstExpr
  .errLines <- rxui$predDf$line
  .expr3 <- .getModelLineEquivalentLhsExpression(lhsExpr)
  .ret <- .getModelineFromExperssionsAndOriginalLines(expr, .expr3, errorLine, .errLines, .origLines)
  if (is.null(.ret)) {
    return(NULL)
  } else if (!is.na(.ret)) {
    return(.ret)
  }
  .getNegativeModelLineForDiffFromProperty(expr, .origLines, errorLine)
}

.thetamodelVars <- rex::rex(or("tv", "t", "pop", "POP", "Pop", "TV", "T", "cov", "err", "eff"))
.thetaModelReg <- rex::rex(or(
  group(start, .thetamodelVars),
  group(.thetamodelVars, end)))

.etaParts <- c(
  "eta", "ETA", "Eta", "ppv", "PPV", "Ppv", "iiv", "Iiv", "bsv", "Bsv", "BSV",
  "bpv", "Bpv", "BPV", "psv", "PSV", "Psv")

.etaModelReg <- rex::rex(or(group(start, or(.etaParts)), group(or(.etaParts), end)))

#' Find the variables in the expression
#'
#' @param x Expression
#' @return Character vector of variables in the expression
#' @author Matthew L. Fidler
#' @noRd
.findVariablesInExpression <- function(x) {
  if (is.atomic(x)) {
    character()
  } else if (is.name(x)) {
    return(as.character(x))
  } else {
    if (is.call(x)) {
      .x1 <- x[-1]
    } else {
      .x1 <- x
    }
    unique(unlist(lapply(.x1, .findVariablesInExpression)))
  }
}

#' @export
rxUiGet.mvFromExpression <- function(x, ...) {
  .x <- x[[1]]
  .exact <- x[[2]]
  eval(call("rxModelVars",as.call(c(list(quote(`{`)), .x$lstExpr[-.x$predDf$line]))))
}
attr(rxUiGet.mvFromExpression, "desc") <- "Calculate model variables from stored (possibly changed) expression"
