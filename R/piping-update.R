#'
#'
#' @param object RxODE UI object
#' @param ... Lines to update
#' @param envir Environment for evaluating ini({}) style calls
#' @export
update.rxUi <- function(object, ..., envir=parent.frame()) {
  # x, ..., envir=parent.frame()
  x <- object
  .ret <- .copyUi(x) # copy so (as expected) old UI isn't affected by the call
  .modelLines <- .quoteCallInfoLines(match.call(expand.dots = TRUE)[-(1:2)], envir=envir)
  .modelHandleModelLines(.modelLines, .ret, modifyIni=TRUE, envir)
}
