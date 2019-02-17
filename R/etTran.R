##'@export
print.rxEtTran <- function(x,...){
    print(as.data.frame(x));
    .cls <- class(x);
    .lst <- attr(.cls, ".RxODE.lst")
    cat("\n Covariates (non time-varying):\n")
    print(.lst$cov1)
}
