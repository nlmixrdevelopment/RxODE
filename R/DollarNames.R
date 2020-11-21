#' @importFrom utils .DollarNames
#' @export
.DollarNames.rxEt <- function(x, pattern) {
  grep(pattern, .Call(`_RxODE_etDollarNames`, x), value = TRUE)
}

#' @export
.DollarNames.rxSolve <- function(x, pattern) {
  grep(pattern, .Call(`_RxODE_rxSolveDollarNames`, x), value = TRUE)
}
