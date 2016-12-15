rxReq <- function(pkg){
    ## nocov start
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package \"%s\" needed for this function to work. Please install it.", pkg), call. = FALSE);
    }
    ## nocov end
}
asTbl <- function(obj){
    if (getOption("RxODE.prefer.tbl", TRUE) && requireNamespace("dplyr", quietly = TRUE)){
        return(dplyr::as.tbl(as.data.frame(obj)));
    } else {
        return(as.data.frame(obj))
    }
}

rxTbl <- function(x, msg){
    if (getOption("RxODE.prefer.tbl", TRUE) && class(x) == "data.frame" && requireNamespace("dplyr", quietly = TRUE)){
        if (!missing(msg)){
            if (getOption("RxODE.verbose", TRUE)){ ## nocov start
                cat(sprintf("Change solved object to dplyr's tbl for %s\n", msg));
            } ## nocov end
        }
        return(dplyr::as.tbl(x))
    } else {
        if (!missing(msg)){
            if (getOption("RxODE.verbose", TRUE)){ ## nocov start
                cat(sprintf("Change solved object to data.frame for %s\n", msg))
            } ## nocov end
        }
        return(x)
    }
}
##' Cleanup anonymous dlls
##'
##' This cleans up any dlls created by text files
##'
##' @param wd What directory should be cleand
##'
##' This cleans up all files named rx-*.dll and associated files as
##' well as call_dvode.o and associated files
##'
##' @return TRUE if successful
##'
##' @author Matthew L. Fidler
##' @export
rxClean <- function(wd = getwd()){
    owd <- getwd();
    setwd(wd);
    on.exit(setwd(owd));
    pat <- "^(Makevars|(rx.*)[.](o|dll|s[ol]|c|rx))$"
    files <- list.files(pattern = pat);
    for (f in files){
        if (f == "Makevars" && file.exists("tran.c")){
            warning("Ignoring Makevars since 'tran.c' is in the same directory.")
        } else {
            try(dyn.unload(f), silent = TRUE);
            unlink(f);
        }
    }
    return(length(list.files(pattern = pat)) == 0);
}
##' RxODE load directory
##'
##' @title RxODE load directory
##' @param ... File name under load directory, or nothing for load
##'     directory alone.
##' @return RxODE load Directory
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxLoadDir <- function(...){
    tmp <- getLoadedDLLs()$RxODE;
    class(tmp) <- "list"
    loadDir <- dirname(tmp$path);
    return(file.path(loadDir, ...))
}

##' Return the include directory
##'
##' The include directory has the headers that may be needed to build
##' functions against the RxODE library.
##'
##' @title RxODE C headers include directory
##' @param ... Additional parameters sent to file.path
##' @return RxODE include directory
##' @author Matthew L. Fidler
##' @export
rxIncludeDir <- function(...){
    incl <- system.file("include", package = "RxODE");
    if (file.exists(file.path(incl, "d.h"))){
        return(file.path(incl, ...));
    } else {
        ## nocov start
        incl <- system.file("src", package = "RxODE");
        if (file.exists(file.path(incl, "d.h"))){
            return(file.path(incl, ...));
        } else {
            stop("Cannot find d.h in a include directory.  RxODE installation may be corrupt.")
        }
        ## nocov end
    }
}

refresh <- function(){
    ## nocov start
    cat("Dparser Version\n");
    print(dparser::dpVersion())
    source(devtools::package_file("build/refresh.R"))
    ## nocov end
}

