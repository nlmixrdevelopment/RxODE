#' @export
update.rxUi <- function(object, ..., envir=parent.frame()) {
  # x, ..., envir=parent.frame()
  x <- object
  .ret <- .copyUi(x) # copy so (as expected) old UI isn't affected by the call
  .lines <- .quoteCallInfoLines(match.call(expand.dots = TRUE)[-(1:2)])
  .modelLines <- vapply(.lines, function(line){
    .cur <- try(.collectWarnings(.iniHandleFixOrUnfix(line, .ret, envir=envir),
                                 lst = TRUE), silent=TRUE)
    if (inherits(.cur, "try-error")) {
      return(TRUE)
    } else {
      lapply(.cur[[2]], function(w){
        warning(w, call.=FALSE)
      })
      return(FALSE)
    }
  }, logical(1))
  .modelHandleModelLines(.lines[which(.modelLines)], .ret)
}
