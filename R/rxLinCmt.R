##' Interface to solved linear compartments implemented in RxODE
##'
##' @param t Times to Solve
##' @param params Parameters for solved linear compartments
##' @param events Event Table for solving
##' @param cmt The compartment for the solving, by default this is one
##' @return A matrix with the solved equations and the derivatives of
##'     the equation based on the parameters you provide.
##'
##' @author Matthew L. Fidler
##' @export
rxLinCmt <- function(t=NULL,
                     params=NULL,
                     events=NULL,
                     cmt=1){
    names(params) <- toupper(names(params));
    lst <- as.list(params);
    lst$cmt      <- as.integer(cmt);
    lst$dosing   <- as.matrix(events$get.dosing());
    if (is.null(t)){
        lst$t <- as.matrix(events$get.sampling()[, 1]);
    } else {
        lst$t <- as.matrix(sort(t), ncol=1);
    }
    env <- list2env(lst);
    linCmtEnv(env);
    return(env$ret);
}
