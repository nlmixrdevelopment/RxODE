##' Returns a list of physical drives that have been or currently are
##' mounted to the computer.
##'
##' This excludes network drives.  See
##' \url{https://www.forensicmag.com/article/2012/06/windows-7-registry-forensics-part-5}
##'
##' @param duplicates Return drives with duplicate entires in
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

##' Return Rtools base
##'
##' @return Rtools base path, or "" on unix-style platforms.
##' @author Matthew L. Fidler
.rxRtoolsBaseWin <- memoise::memoise(function(){
    if (.Platform$OS.type == "unix"){
        return("");
    } else {
        ## Prefer nlmixr rtools over everything
        .keys <- .rxNlmixr()
        if (!is.null(.keys$rtoolsBase)){
            .rtoolsBase <- .keys$rtoolsBase;
        } else {
            ## The grep solution assumes that the path is setup correctly;
            .gcc <- Sys.which("gcc.exe")
            .rtools <- sub("[/\\](mingw).*", "", .gcc);
            if (file.exists(file.path(.rtools, "Rtools.txt"))){
                return(.rtools)
            } else {
                ## Rtools doesn't add itself to the path by default.  To
                ## remove install headaches, fish for the path a bit.

                ## The general solution also corrects the problem of
                ## having msys or cygwin compilers on top of the Rtools
                ## compiler, and will adjust the path (just because which
                ## shows a different path doesn't mean Rtools isn't
                ## there.)
                ## This is what Rtools installer is supposed to do.
                ## There is some discussion on devtools if this really occurs...
                .rtoolsBase <- "C:/Rtools";
                if (!file.exists(.rtoolsBase)) {
                    .keys <- try(utils::readRegistry("SOFTWARE\\R-core\\Rtools", hive = "HCU",
                                                    view = "32-bit", maxdepth = 2), silent = TRUE)
                    if (is.null(.keys) || length(.keys) == 0)
                        .keys <- try(utils::readRegistry("SOFTWARE\\R-core\\Rtools", hive = "HLM",
                                                        view = "32-bit", maxdepth = 2), silent = TRUE)
                    if (!inherits(.keys, "try-error")){
                        for (i in seq_along(.keys)) {
                            .version <- names(.keys)[[i]]
                            .key <- .keys[[.version]]
                            if (!is.list(.key) || is.null(.key$InstallPath)) next;
                            install_path <- .normalizePath(.key$InstallPath, mustWork = FALSE,
                                                          winslash = "/");
                            if (file.exists(install_path)){
                                .rtoolsBase <- install_path;
                            }
                        }
                    }
                }
                .ver <- R.Version();
                .ver <- paste0(.ver$major, ".", gsub(rex::rex(start, capture(except_any_of(".")), ".",
                                                            anything, end), "\\1", .ver$minor))
                if (!file.exists(.rtoolsBase)){## Based on Issue #2, Rtools may also be installed to RBuildTools;  This is also reflected on the R-stan website.
                    .rtoolslist <- apply(expand.grid(c("Rtools", paste0("Rtools/", .ver),
                                                      "RBuildTools", paste0("RBuildTools/", .ver)), rxPhysicalDrives()), 1,
                                        function(x){ paste0(x[2], x[1])});
                    for (.path in .rtoolslist){
                        if (file.exists(.path)){
                            return(.path)
                        }
                    }
                }
                if (file.exists(.rtoolsBase)){
                    return(.rtoolsBase)
                } else if (file.exists(.rtools)) {
                    message("gcc available, assuming it comes from rtools...\nRxODE may not work with other compilers.\n")
                    return(.rtools)
                } else {
                    message("This package requires Rtools!\nPlease download from http://cran.r-project.org/bin/windows/Rtools/,\ninstall and restart your R session before proceeding.")
                    return("c:/Rtools")
                }
            }
        }
    }
})
##' Setup Rtools path
##'
##' @param rm.rtools Remove the Rtools from the current path specs.
##' @param rm.python Remove Python from the current path specs.
##'
##' @author Matthew L. Fidler
.rxWinRtoolsPath <- function(rm.rtools=TRUE, rm.python=TRUE){
    ## Note that devtools seems to assume that rtools/bin is setup
    ## appropriately, and figures out the c compiler from there.
    if (.Platform$OS.type == "unix"){
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
        .rtoolsBase <- .rxRtoolsBaseWin();
        .x <- file.path(.rtoolsBase, ifelse(.Platform$r_arch == "i386","mingw_32/bin", "mingw_64/bin"));
        if (file.exists(.x)){
            Sys.setenv(BINPREF=gsub("([^/])$", "\\1/", gsub("\\\\", "/", .normalizePath(.x))));
        }
        if (file.exists(.rtoolsBase)){
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
            ## Last Cran check for Rtools is qpdf
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
##' Setup Python and SymPy for windows
##'
##' @author Matthew L. Fidler
##' @export
rxWinPythonSetup <- function(){
    .base <- .rxPythonBaseWin()
    if (is.null(.base)){
        stop("RxODE requires Python. Please install an appropriate version and add it to the system path.")
    }
    if (file.access(paste(.base, "/Lib/site-packages", sep=""),2)==-1){
      stop("The Python library path does not appear to be writeable. Please rectify this situation, restart R, and try again.")
    }
    message("Attempting to install SymPy. This may take a few seconds...")
    try(system("python -m pip install sympy"))

    message("Please restart your R session before using RxODE.")
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
