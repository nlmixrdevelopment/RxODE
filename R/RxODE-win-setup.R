##' Returns a list of physical drives that have been or currently are
##' mounted to the computer.
##'
##' This excludes network drives.  See
##' \url{https://www.forensicmag.com/article/2012/06/windows-7-registry-forensics-part-5}
##'
##' @param duplicates Return drives with duplicate entries in
##'     \code{SYSTEM\\MountedDevices}; These are likely removable media.  By default this is \code{FALSE}
##' @return Drives with letters
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxPhysicalDrives <- memoise::memoise(function(duplicates=FALSE){
                                 if (.Platform$OS.type == "unix"){
        return(NULL)
    } else {
        ## This lists all the drive letters (and volume
        ## information) of drives mounted to your computer.
        n <- names(utils::readRegistry("SYSTEM\\MountedDevices"))
        reg <- rex::rex(start, "\\DosDevices\\", capture(or("A":"Z", "a":"z"), ":"), end)
        ns <- n[regexpr(reg, n) != -1];
        if (length(n) > 0){
            ns <- toupper(gsub(reg, "\\1", ns));
            dups <- unique(ns[duplicated(ns)]);
            if (length(dups) > 1){
                if (duplicates){
                    ## Duplicate drive names are more likely to be removable media letters (like usb/cd/etc.)
                    d <- paste0(sort(unique(dups)), "\\")
                    w <- which(!sapply(d, removableDrive))
                    if (length(d) >= 1){
                        return(d)
                    } else {
                        return("C:\\")
                    }
                } else {
                    d <- paste0(sort(unique(ns[!(ns %in% dups)])), "\\");
                    w <- which(!sapply(d, removableDrive))
                    d <- d[w]
                    if (length(d) >= 1){
                        return(d)
                    } else {
                        return("C:\\")
                    }
                }
            } else {
                d <- paste0(sort(ns), "\\");
                w <- which(!sapply(d, removableDrive));
                d <- d[w]
                if (length(d) >= 1){
                    return(d)
                } else {
                    return("C:\\")
                }
            }
            ret <- ns;
        } else {
            ret <- "C:\\";
        }
        return(ret)
    }
})
.rxPythonBaseWin <- function(){
    if (.Platform$OS.type == "unix"){
    } else if (file.exists(R.home("rxPython"))){
        R.home("rxPython")
    } else {
        .keys <- .rxNlmixr();
        if (!is.null(.keys$pythonBase)){
            .pythonBase <- .keys$pythonBase;
        } else {
            .keys <- NULL
            .keys <- try(utils::readRegistry("SOFTWARE\\Python\\PythonCore", hive = "HCU", ## view = "32-bit",
                                            maxdepth = 3), silent=TRUE);
            if (is.null(.keys) || length(.keys) == 0 || inherits(.keys, "try-error")){
                try(.keys <- utils::readRegistry("SOFTWARE\\Python\\PythonCore", hive = "HLM", ## view = "32-bit",
                                                maxdepth = 3), silent = TRUE);
            }
            .pythonBase <- NULL
            for (i in seq_along(.keys)){
                try(.pythonBase <- .keys[[i]]$InstallPath$`(Default)`, silent=TRUE)
                if (!is.null(.pythonBase)){
                    if (file.exists(file.path(.pythonBase, "python.exe"))){
                        break;
                    } else {
                        .pythonBase <- NULL;
                    }
                }
            }
        }
        if (!is.null(.pythonBase)){
            .pythonBase <- gsub("\\", "/", .normalizePath(gsub(rex::rex(any_of("/", "\\"), end), "", .pythonBase)), fixed=TRUE);
        }
        return(.pythonBase)
    }
}

.rxNlmixr <- memoise::memoise(function(){
    .keys <- NULL
    .keys <- try(utils::readRegistry(sprintf("SOFTWARE\\nlmixr%s", ifelse(.Platform$r_arch == "i386", "32", "")),
                                    hive = "HCU", maxdepth = 2), silent = TRUE);
    if (!inherits(.keys, "try-error")){
        .ca1 <- commandArgs()[[1]];
        if (.ca1 != "RStudio"){
            .ca1 <- .normalizePath(.ca1);
            if (regexpr(rex::rex(substring(.normalizePath(.keys[[1]]), 2)), .ca1) == -1){

                return(list())
            }
        }
        .lst <- list(
            rtoolsBase = .normalizePath(file.path(.keys[[1]], "rtools"), winslash="/", mustWork=FALSE),
            pythonBase=.normalizePath(file.path(.keys[[1]], "python"), winslash="/", mustWork=FALSE)
        );
        if (dir.exists(.lst$rtoolsBase) && dir.exists(.lst$pythonBase)){
            return(.lst);
        }
    }
    return(list());
})

.rxRtoolsBaseWin <- memoise::memoise(function(retry=FALSE){
    if (.Platform$OS.type == "unix" || getOption("RxODE.rtools.setup", FALSE)){
        return("");
    } else if (file.exists(R.home("rtools"))){
        return(R.home("rtools"))
    } else {
        return(NULL)
    }
})
##' Setup Rtools path
##'
##' @param rm.rtools Remove the Rtools from the current path specs.
##'
##' @param rm.python Remove Python from the current path specs.
##'
##' @param retry Should you retry to find Rtools?  If NA, don't throw
##'     an error if it isn't found.
##'
##' @author Matthew L. Fidler
.rxWinRtoolsPath <- function(rm.rtools=TRUE, rm.python=TRUE, retry=FALSE){
    ## Note that devtools seems to assume that rtools/bin is setup
    ## appropriately, and figures out the c compiler from there.
    if (.Platform$OS.type == "unix" || getOption("RxODE.rtools.setup", FALSE)){
        return(TRUE)
    } else if (Sys.which("make") != ""){
        return(TRUE)
    } else {
        .path <- unique(sapply(sub(rex::rex('"', end), "", sub(rex::rex(start, '"'), "",
                                   gsub("/", "\\\\", strsplit(Sys.getenv("PATH"), ";")[[1]]))),
                              function(x){
            if (file.exists(x)){
                return(.normalizePath(x));
            } else {
                return("");
            }
        }))
        .path <- .path[.path != ""];
        if (!inherits(rm.rtools, "logical")) rm.rtools <- FALSE
        if (!inherits(rm.python, "logical")) rm.python <- FALSE
        if (rm.rtools){
            .path <- .path[regexpr(rex::rex(or("Rtools", "RTOOLS", "rtools")), .path) == -1]
        }
        if (rm.python){
            .path <- .path[regexpr(rex::rex(or("Python", "python", "PYTHON")), .path) == -1]
        }
        .rPath <- .normalizePath(file.path(Sys.getenv("R_HOME"), paste0("bin", Sys.getenv("R_ARCH"))));
        .path <- c(.rPath, .path);

        ## Look in the registry...
        ## This is taken from devtools and adapted.
        .rtoolsBase <- .rxRtoolsBaseWin(retry=retry);
        if (is.null(.rtoolsBase)) return("");
        .x <- file.path(.rtoolsBase, ifelse(.Platform$r_arch == "i386","mingw_32/bin", "mingw_64/bin"));
        if (file.exists(.x)){
            Sys.setenv(rxBINPREF=gsub("([^/])$", "\\1/", gsub("\\\\", "/", .normalizePath(.x))));
        }
        .exists <- try(file.exists(.rtoolsBase), silent=TRUE);
        if (inherits(.exists, "try-error")) .exists <- FALSE
        if (.exists){
            .gcc <- list.files(.rtoolsBase, "gcc",full.names=TRUE)[1]
            if (is.na(.gcc)){
                .gcc <- "";
            }
            ## This allows both toolchains to be present, but RxODE should still work...
            for (.x in rev(c(file.path(.rtoolsBase, "bin"),
                            file.path(.rtoolsBase, "pandoc"),
                            file.path(.rtoolsBase, "qpdf/bin"),
                            ## file.path(.rtoolsBase, "mingw_32/bin") ## Rtools sets up the mingw_32/bin first (even if x64)
                            file.path(.rtoolsBase, ifelse(.Platform$r_arch == "i386","mingw_32/bin", "mingw_64/bin")),
                            file.path(.rtoolsBase, ifelse(.Platform$r_arch == "i386","mingw_32/bin", "mingw_64/bin")),
                            file.path(.rtoolsBase, ifelse(.Platform$r_arch == "i386","mingw_32/opt/bin", "mingw_64/opt/bin"))
                            ## ifelse(.gcc == "", "", file.path(.gcc, "bin")),
                            ## ifelse(.gcc == "", "", ifelse(.Platform$r_arch == "i386",file.path(.gcc, "bin32"), file.path(.gcc, "bin64"))
                            ## )
                            ))){
                if (file.exists(.x)){
                    .path <- c(.normalizePath(.x), .path);
                }
            }
            .pythonBase <- .rxPythonBaseWin();
            if (!is.null(.pythonBase)){
                python <- .normalizePath(file.path(.pythonBase, "python.exe"));
                if (file.exists(python)){
                    ## Sometimes there are 2 competing python
                    ## installs.  Make sure to setup everything, just
                    ## in case... Otherwise it can crash R :(
                    Sys.setenv(PYTHON_EXE=.normalizePath(python)); ## For PythonInR
                    .path <- c(.normalizePath(.pythonBase), .path);
                    Sys.setenv(PYTHONHOME=.normalizePath(.pythonBase));
                    Sys.setenv(PYTHON_INCLUDE=.normalizePath(file.path(.pythonBase, "include")));
                    .pythonLibBase <- .normalizePath(file.path(.pythonBase, "libs"))
                    .pythonLibName <- list.files(.pythonLibBase, "[0-9][0-9]+\\.lib$");
                    if (length(.pythonLibName) == 1){
                        Sys.setenv(PYTHON_LIB=.normalizePath(file.path(.pythonBase, "libs", .pythonLibName)));
                    } else {
                        Sys.unsetenv("PYTHON_LIB")
                    }
                    Sys.setenv(PYTHONPATH=paste(c(.normalizePath(file.path(.pythonBase, "DLLs")),
                                                  .normalizePath(file.path(.pythonBase, "Lib")),
                                                  .normalizePath(file.path(.pythonBase, "Lib", "site-packages"))),
                                                collapse=";"));
                    Sys.unsetenv("PYTHONSTARTUP");
                }
            }
            ## java <- as.vector(Sys.which("java"));
            ## if (java != ""){
            ##     java <- sub(rex::rex(one_of("/", "\\"), except_any_of("/", "\\", "\n"), end), "", java)
            ## }
            .keys <- NULL;
            ## Is there a 64 bit aspell that should be checked for...?
            .keys <- try(utils::readRegistry("SOFTWARE\\Aspell", hive="HLM", view="32-bit", maxdepth=3), silent=TRUE);
            ## Add aspell for cran checks...
            if (!is.null(.keys)){
                if (any(names(.keys) == "Path")){
                    if (file.exists(.keys$Path)){
                        .path <- c(.normalizePath(.keys$Path), .path);
                    }
                }
            }
            ## Last CRAN check for Rtools is qpdf
            .qpdf <- c(paste0(.rtoolsBase, "/qpdf/bin"), paste0(rxPhysicalDrives(), "/qpdf/bin"))
            for (.p in .qpdf){
                if (file.exists(.p)){
                    .path <- c(.normalizePath(.p), .path);
                    break;
                }
            }
            .path <- .path[.path != ""];
            .path <- paste(unique(.path), collapse=";");
            Sys.setenv(PATH=.path);
            return(TRUE);
        } else {
            return(FALSE)
        }
    }
}
.installRxPython <- function(dir=R.home("rxPython"),
                             opt="InstallAllUsers=0 AssociateFiles=0 Shortcuts=0"){
    dir <- .normalizePath(dir);
    opt <- sprintf("%s DefaultJustForMeTargetDir=%s", opt, dir, dir)
    .x64 <- (regexpr("x64", Sys.info()["release"]) != -1)
    ## if (x64) opt <- sprintf("/passive %s", opt)
    .tmp <- try({installr::install.python(installer_option=opt, x64=.x64)}, silent=TRUE)
    if (inherits(.tmp, "try-error")){
        ## Installr is a bit old in MRAN support that possibility
        installr::install.python(installer_option=opt)
    }
    system(sprintf("%s/python -m pip install --upgrade pip", R.home("rxPython")))
    system(sprintf("%s/python -m pip install --upgrade sympy", R.home("rxPython")))
    system(sprintf("%s/python -m pip install --upgrade numpy", R.home("rxPython")))
}
##' Setup Python and SymPy for windows
##'
##' @author Matthew L. Fidler
##' @export
rxWinPythonSetup <- function(){
    .base <- .rxPythonBaseWin()
    if (is.null(.base)){
        .installRxPython()
        .base <- .rxPythonBaseWin()
        if (is.null(.base)){
            stop("RxODE requires Python. Please install an appropriate version and add it to the system path.")
        }
    }
    if (file.access(paste(.base, "/Lib/site-packages", sep=""),2)==-1){
      stop("The Python library path does not appear to be writeable. Please rectify this situation, restart R, and try again.")
    }
    .tmp <- try({rxSymPyVersion()})
    if (inherits(.tmp, "try-error")){
        system(sprintf("%s/python -m pip install --upgrade pip", .rxPythonBaseWin()))
        system(sprintf("%s/python -m pip install --upgrade sympy", .rxPythonBaseWin()))
        system(sprintf("%s/python -m pip install --upgrade numpy", .rxPythonBaseWin()))
    }
    utils::install.packages("reticulate");
    .tmp <- try({rxSymPyVersion()});
    if (inherits(.tmp, "try-error")){
        rxReq("remotes")
        remotes::install_github("nlmixrdevelopment/SnakeCharmR")
        .tmp <- try({rxSymPyVersion()});
        if (inherits(.tmp, "try-errror")){
            stop("Cannot setup RxODE<->sympy link");
        }
    }
}

##' Setup Windows components for RxODE
##'
##' @inheritParams .rxWinRtoolsPath
##' @author Matthew L. Fidler
##' @export
rxWinSetup <- function(rm.rtools=TRUE, rm.python=TRUE){
    if (!.rxWinRtoolsPath(rm.rtools, rm.python)){
        message("This package requires Rtools!\nPlease download from http://cran.r-project.org/bin/windows/Rtools/,\ninstall and restart your R session before proceeding.\n");
    }
    rxWinPythonSetup();
}
