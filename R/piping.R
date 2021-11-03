#'  This copies the RxODE UI object so it can be modifed
#'
#' @param ui Original UI object
#' @return Copied UI object
#' @author Matthew L. Fidler
#' @noRd
.copyUi <- function(ui) {
  .ret <- new.env(parent=emptyenv())
  lapply(ls(ui, envir=ui, all.names=TRUE), function(item){
    assign(item, get(item, envir=ui), envir=.ret)
  })
  class(.ret) <- class(ui)
  .ret
}

#' @export
#' @rdname ini
ini.function <- function(x, ..., envir=parent.frame()) {
  .ui <- RxODE(x)
  ini(x=.ui, ..., envir=envir)
}

#'  Returns quoted call information
#'
#' @param callInfo Call information
#' @return Quote call information.  for `name=expression`, change to
#'   `name<-expression` in quoted call list
#' @author Matthew L. Fidler
#' @noRd
.iniQuoteLines <- function(callInfo) {
  lapply(seq_along(callInfo), function(i) {
    .name <- names(callInfo)[i]
    if (!is.null(.name)) {
      if (.name != "") {
        # Changed named items to
        # quote(name <- expression)
        return(as.call(list(quote(`<-`), .enQuote(.name),
                            eval(call("quote", callInfo[[i]])))))
      }
    }
    eval(call("quote", callInfo[i]))
  })
}

.iniModifyThetaDf <- function(ini, lhs, rhs, doFix, doUnfix) {
  .w <- which(ini$name == lhs)
  if (length(.w) != 1) {
    stop("Error updating parameter name '", lhs, "'", call.=FALSE)
  }
  .curFix <- ini$fix[.w]
  if (doFix) {
    if (.curFix) {
      warning("trying to fix '", lhs, "', but already fixed")
    } else {
      ini$fix[.w] <- TRUE
    }
  } else if (doUnfix) {
    if (.curFix) {
      ini$fix[.w] <- FALSE
    } else {
      warning("trying to unfix '", lhs, "', but already unfixed")
    }
  }
  if (is.null(rhs)) {
  } else if (length(rhs) == 1)  {
    ini$est[.w] <- rhs
  } else if (length(rhs) == 2) {
    ini$lower[.w] <- rhs[1]
    ini$est[.w] <- rhs[2]
  } else if (length(rhs) == 3) {
    ini$lower[.w] <- rhs[1]
    ini$est[.w] <- rhs[2]
    ini$upper[.w] <- rhs[3]
  } else {
    stop("piping for '", lhs, "' failed, the estimates should be between 1-3 values")
  }
  ini
}

.iniHandleFixOrUnfix <- function(expr, rxui, envir=parent.frame()) {
  .assignOp <- expr[[1]]
  if (identical(.assignOp, quote(`<-`)) ||
        identical(.assignOp, quote(`=`))) {
    .lhs <- as.character(expr[[2]])
    .rhs <- expr[[3]]
    .doFix <- .doUnfix <- FALSE
    if (is.name(.rhs)) {
      if (identical(.rhs, quote(`fix`)) ||
          identical(.rhs, quote(`fixed`)) ||
          identical(.rhs, quote(`FIX`)) ||
          identical(.rhs, quote(`FIXED`))) {
        .doFix <- TRUE
      } else if (identical(.rhs, quote(`unfix`)) ||
                   identical(.rhs, quote(`unfixed`)) ||
                   identical(.rhs, quote(`UNFIX`)) ||
                   identical(.rhs, quote(`UNFIXED`))) {
        .doUnfix <- TRUE
      }
      .rhs <- NULL
    } else if (identical(.rhs[[1]], quote(`fix`)) ||
          identical(.rhs[[1]], quote(`fixed`)) ||
          identical(.rhs[[1]], quote(`FIX`)) ||
          identical(.rhs[[1]], quote(`FIXED`))) {
      .doFix <- TRUE
      .rhs[[1]] <- quote(`c`)
    } else if (identical(.rhs[[1]], quote(`unfix`)) ||
          identical(.rhs[[1]], quote(`unfixed`)) ||
          identical(.rhs[[1]], quote(`UNFIX`)) ||
          identical(.rhs[[1]], quote(`UNFIXED`))) {
      .doUnfix <- TRUE
      .rhs[[1]] <- quote(`c`)
    }
    if (!is.null(.rhs)) {
      .rhs <- eval(.rhs, envir=envir)
    }
    assign("iniDf", .iniModifyThetaDf(rxui$ini, .lhs, .rhs, .doFix, .doUnfix),
           envir=rxui)
    return(invisible())
  }
}

#' @export
#' @rdname ini
ini.rxUi <- function(x, ..., envir=parent.frame()) {
  .ret <- .copyUi(x) # copy so (as expected) old UI isn't affected by the call
  .iniLines <- .iniQuoteLines(match.call(expand.dots = TRUE)[-(1:2)])
  lapply(.iniLines, function(line){
    .iniHandleFixOrUnfix(line, .ret, envir=envir)
  })
  .ret
}
