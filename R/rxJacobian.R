## Refactor Jacobian calculation
##' Faster expand.grid
##'
##' Only support x and y as characters right now
##'
##' @param x first element (must be character)
##' @param y second element (must be character)
##' @return Expand grid
##' @author Matthew Fidler
##' @keywords internal
##' @examples
##'
##' ##
##' rxExpandGrid(letters, letters)
##'
##' ## Another fast method; See
##' ## https://stackoverflow.com/questions/10405637/use-outer-instead-of-expand-grid
##'
##' expand.grid.jc <- function(seq1,seq2) {
##'  cbind(Var1 = rep.int(seq1, length(seq2)),
##'   Var2 = rep.int(seq2, rep.int(length(seq1),length(seq2))))
##' }
##'
##' \dontrun{
##'  microbenchmark::microbenchmark(rxExpandGrid(letters, letters), expand.grid.jc(letters, letters))
##' }
##' @export
rxExpandGrid <- function(x, y){
    rxExpandGrid_(x, y)
}

