##' Validate RxODE
##'
##' This allows easy vaildation/qualification of nlmixr by running the
##' testing suite on your system.
##' @param type Type of test or fitler of test type
##' @author Matthew L. Fidler
##' @export
rxValidate <- function(type=NULL){
    ## rxVersion(" Validation", TRUE);
    if (is.character(type)){
        if (type=="covr"){
            Sys.setenv("NOT_CRAN"="true", "covr"="true");
            covr::report()
        } else {
            old.wd <- getwd();
            on.exit({setwd(old.wd); Sys.setenv(NOT_CRAN="")});
            Sys.setenv("NOT_CRAN"="true")
            path <- file.path(system.file("tests", package = "RxODE"),"testthat")
            setwd(path)
            testthat::test_dir(path, filter=type);
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
    }
}

##' @rdname rxValidate
##' @export
rxTest <- rxValidate
