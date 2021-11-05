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
