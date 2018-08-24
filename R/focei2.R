.foceiSetup <- function(obj, data, theta, thetaFixed = NULL,
                        skipCov=NULL, rxInv = NULL,
                        lower = NULL, upper = NULL, etaMat = NULL,
                        control = RxODE::foceiControl()){
    loadNamespace("n1qn1");
    .Call(`_RxODE_foceiSetup_`, obj, data, theta, thetaFixed, skipCov, rxInv, lower, upper, etaMat, control); # nolint
}


.nearPd <- function(mat){
    if (any(is.na(mat))){
        ## cat("Bad matrix:\n");
        ## print(mat);
        return(mat)
    } else {
        return(as.matrix(Matrix::nearPD(mat)$mat));
    }
}
##' Cox Box transformation
##'
##' @param x data to transform
##' @param lambda Cox-box lambda parameter
##' @return Cox-Box Transformed Data
##' @author Matthew L. Fidler
##' @export
coxBox <- function(x, lambda=1){
    .Call(`_RxODE_coxBox_`, x, lambda, 0L)
}

##' Yeo-Johnson Transformation
##'
##' @param x data to transform
##' @param lambda Cox-box lambda parameter
##' @return Yeo-Johnson  Transformed Data
##' @author Matthew L. Fidler
##' @export
yeoJohnson <- function(x, lambda=1){
    .Call(`_RxODE_coxBox_`, x, lambda, 1L)
}
