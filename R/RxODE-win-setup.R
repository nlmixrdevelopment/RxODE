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
            rtoolsBase = .normalizePath(file.path(.keys[[1]], "rtools"), winslash="/", mustWork=FALSE));
        if (dir.exists(.lst$rtoolsBase)){
            return(.lst);
        }
    }
    return(list());
})

.rxRtoolsBaseWin <- memoise::memoise(function(retry=FALSE){
    if (.Platform$OS.type == "unix" || getOption("RxODE.rtools.setup", FALSE)){
        return("");
    } else {
        for (.i in rxPhysicalDrives()){
            if (file.exists(paste0(.i, "Rtools"))){
                return(paste0(.i, "Rtools"))
            }
        }
        if (file.exists(R.home("rtools"))){
            return(R.home("rtools"))
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
                    .exists <- try(file.exists(.rtoolsBase), silent=TRUE);
                    if (inherits(.exists, "try-error")) .exists <- FALSE
                    if (!.exists) {
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
                    .exists <- try(file.exists(.rtoolsBase), silent=TRUE);
                    if (inherits(.exists, "try-error")) .exists <- FALSE
                    if (!.exists){## Based on Issue #2, Rtools may also be installed to RBuildTools;  This is also reflected on the R-stan website.
                        .rtoolslist <- apply(expand.grid(c("Rtools", paste0("Rtools/", .ver),
                                                           "RBuildTools", paste0("RBuildTools/", .ver)), rxPhysicalDrives()), 1,
                                             function(x){ paste0(x[2], x[1])});
                        for (.path in .rtoolslist){
                            if (file.exists(.path)){
                                return(.path)
                            }
                        }
                    }
                    .exists <- try(file.exists(.rtoolsBase), silent=TRUE);
                    if (inherits(.exists, "try-error")) .exists <- FALSE
                    if (.exists){
                        return(.rtoolsBase)
                    } else if (file.exists(.rtools)) {
                        message("'gcc' available, assuming it comes from 'rtools'");
                        message("RxODE may not work with other compilers")
                        return(.rtools)
                    } else {
                        if (is.na(retry)) return(FALSE)
                        if (file.exists(file.path(Sys.getenv("BINPREF"), "gcc"))){
                            return(NULL)
                        }
                        if (retry){
                            stop("RxODE requires 'Rtools'\nPlease download from http://cran.r-project.org/bin/windows/Rtools/,\ninstall and restart your R session before proceeding")
                        }
                        try({installr::install.Rtools()});
                        return(.rxRtoolsBaseWin(retry=TRUE))
                    }
                }
            }
        }
    }

})
##' Setup Rtools path
##'
##' @param rm.rtools Remove the Rtools from the current path specs.
##'
##' @param retry Should you retry to find Rtools?  If NA, don't throw
##'     an error if it isn't found.
##'
##' @author Matthew L. Fidler
.rxWinRtoolsPath <- function(rm.rtools=TRUE, retry=FALSE){
    ## Note that devtools seems to assume that rtools/bin is setup
    ## appropriately, and figures out the c compiler from there.
    if (.Platform$OS.type == "unix" || getOption("RxODE.rtools.setup", FALSE)){
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
        if (rm.rtools){
            .path <- .path[regexpr(rex::rex(or("Rtools", "RTOOLS", "rtools")), .path) == -1]
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
##' Setup Windows components for RxODE
##'
##' @inheritParams .rxWinRtoolsPath
##' @author Matthew L. Fidler
##' @export
rxWinSetup <- function(rm.rtools=TRUE){
    if (!.rxWinRtoolsPath(rm.rtools)){
        message("RxODE requires 'rtools'\nPlease download from http://cran.r-project.org/bin/windows/Rtools/,\ninstall and restart your R session before proceeding\n");
    }
}
