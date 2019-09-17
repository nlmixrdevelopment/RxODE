##' Validate RxODE
##'
##' This allows easy validation/qualification of nlmixr by running the
##' testing suite on your system.
##' @param full Should a full validation be performed?  (By default
##'     \code{TRUE})
##' @author Matthew L. Fidler
##' @export
rxValidate <- function(full=TRUE){
    ## rxVersion(" Validation", TRUE);
    if (is.character(full)){
        if (full=="covr"){
            Sys.setenv("NOT_CRAN"="true", "covr"="true");
            covr::report()
        } else {
            old.wd <- getwd();
            on.exit({setwd(old.wd); Sys.setenv(RxODE_VALIDATION_FULL="false", NOT_CRAN="")});
            Sys.setenv(RxODE_VALIDATION_FULL="false", "NOT_CRAN"="true")
            path <- file.path(system.file("tests", package = "RxODE"),"testthat")
            setwd(path)
            testthat::test_dir(path, filter=full);
            Sys.setenv(RxODE_VALIDATION_FULL="true")
            Sys.getenv("RxODE_VALIDATION_FULL")
            testthat::test_dir(path, filter=full);
        }
    } else {
        old.wd <- getwd();
        on.exit({setwd(old.wd)});
        path <- file.path(system.file("tests", package = "RxODE"),"testthat")
        setwd(path)
        Sys.setenv("NOT_CRAN"="funny")
        message("CRAN only tests")
        message("================================================================================")
        pt <- proc.time();
        testthat::test_dir(path);
        message("================================================================================")
        message("Timing of CRAN tests (should be under 60 seconds)")
        message("================================================================================")
        print(proc.time() - pt);
        message("================================================================================")
        message("Normal tests")
        message("================================================================================")
        Sys.setenv("NOT_CRAN"="true")
        testthat::test_dir(path);
        if (full){
            message("================================================================================")
            message("Full Validation tests")
            message("================================================================================")
            Sys.setenv(RxODE_VALIDATION_FULL="true")
            testthat::test_dir(path);
        }
    }
}

##' @rdname rxValidate
##' @export
rxTest <- rxValidate
