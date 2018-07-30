.foceiSetup <- function(obj, data, theta, thetaFixed = NULL,
                        skipCov=NULL, rxInv = NULL,
                        lower = NULL, upper = NULL, etaMat = NULL,
                        control = RxODE::foceiControl()){
    loadNamespace("n1qn1");
    .Call(`_RxODE_foceiSetup_`, obj, data, theta, thetaFixed, skipCov, rxInv, lower, upper, etaMat, control); # nolint
}


.nearPd <- function(mat){
    if (any(is.na(mat))){
        cat("Bad matrix:\n");
        print(mat);
        return(mat)
    } else {
        return(as.matrix(Matrix::nearPD(mat)$mat));
    }
}
