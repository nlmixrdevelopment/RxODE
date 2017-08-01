rxWget <- function(url, to){
    if (Sys.which("wget") == ""){
        if (.Platform$r_arch == "i386"){
            if (file.exists("c:/RTOOLS/mingw_32/bin")){
                download.file("https://eternallybored.org/misc/wget/current/wget.exe", "c:/RTOOLS/mingw_32/bin/wget.exe");
                if (!file.exists("c:/RTOOLS/mingw_32/bin/wget.exe")){
                    stop("Cannot get wget...");
                }
            } else {
                download.file("https://eternallybored.org/misc/wget/current/wget.exe", "wget.exe");

            }
        } else {
            if (file.exists("c:/RTOOLS/mingw_64/bin")){
                download.file("https://eternallybored.org/misc/wget/current/wget64.exe", "c:/RTOOLS/mingw_64/bin/wget.exe");
                if (!file.exists("c:/RTOOLS/mingw_64/bin/wget.exe")){
                    stop("Cannot get wget...");
                }
            } else {
                download.file("https://eternallybored.org/misc/wget/current/wget64.exe", "wget.exe");
            }
        }
    }
    if (Sys.which("wget") == ""){
        stop("Wget not working...");
    }
    if (Sys.which("wget") != ""){
        download.file(url, to, method="wget", extra="--progress=dot --no-check-certificate");
    }
}

rxPythonBaseWin <- function(){
    if(.Platform$OS.type == "unix"){
    } else {
        keys <- NULL;
        keys <- try(utils::readRegistry("SOFTWARE\\Python\\PythonCore", hive = "HCU", ## view = "32-bit",
                                        maxdepth = 3), silent=TRUE);
        if (is.null(keys) || length(keys) == 0 || inherits(keys, "try-error")){
            try(keys <- utils::readRegistry("SOFTWARE\\Python\\PythonCore", hive = "HLM", ## view = "32-bit",
                                            maxdepth = 3), silent = TRUE);
        }
        python.base <- NULL
        for (i in seq_along(keys)){
            try(python.base <- keys[[i]]$InstallPath$`(Default)`, silent=TRUE)
            if (!is.null(python.base)){
                if (file.exists(file.path(python.base, "python.exe"))){
                    break;
                } else {
                    python.base <- NULL;
                }
            }
        }
        if (!is.null(python.base)){
            python.base <- gsub("\\", "/", utils::shortPathName(gsub(rex::rex(any_of("/", "\\"), end), "", python.base)), fixed=TRUE);
        }
        return(python.base)
    }
}
##' Return Rtools base
##'
##' @return Rtools base path, or "" on unix-style platforms.
##' @author Matthew L. Fidler
rxRtoolsBaseWin <- function(){
    if(.Platform$OS.type == "unix"){
        return("");
    } else {
        rtools.base <- "C:/Rtools";
        if (!file.exists(rtools.base)){
            keys <- NULL;
            keys <- utils::readRegistry("SOFTWARE\\R-core\\Rtools", hive = "HCU", view = "32-bit", maxdepth = 2);
            if (is.null(keys) || length(keys) == 0)
                try(keys <- utils::readRegistry("SOFTWARE\\R-core\\Rtools", hive = "HLM", view = "32-bit", maxdepth = 2), silent = TRUE);
            if (is.null(keys) || length(keys) == 0){
                stop("Cannot use this package because Rtools isn't setup appropriately...")
            }

            for(i in seq_along(keys)) {
                version <- names(keys)[[i]]
                key <- keys[[version]]
                if (!is.list(key) || is.null(key$InstallPath)) next;
                install_path <- normalizePath(key$InstallPath, mustWork = FALSE, winslash = "/");
                if (file.exists(install_path)){
                    rtools.base <- install_path;
                }
            }
        }
        return(rtools.base)
    }
}
##' Setup Rtools path
##'
##' @param rm.rtools Remove the Rtools from the current path specs.
##'
##' @author Matthew L. Fidler
rxWinRtoolsPath <- function(rm.rtools=TRUE){
    ## Note that devtools seems to assume that rtools/bin is setup appropriately, and figures out the c compiler from there.
    if(.Platform$OS.type == "unix"){
        return(TRUE)
    } else {
        path <- unique(sapply(gsub("/", "\\\\", strsplit(Sys.getenv("PATH"), ";")[[1]]), function(x){
            if (file.exists(x)){
                return(normalizePath(x));
            } else {
                return("");
            }
        }))
        path <- path[path != ""];
        if (rm.rtools){
            path <- path[regexpr(rex::rex(or("Rtools", "RTOOLS", "rtools")), path) == -1]
        }
        r.path <- normalizePath(file.path(Sys.getenv("R_HOME"),paste0("bin",Sys.getenv("R_ARCH"))));
        path <- c(r.path, path);
        ## Look in the registry...
        ## This is taken from devtools and adapted.
        rtools.base <- rxRtoolsBaseWin();
        if (file.exists(rtools.base)){
            gcc <- list.files(rtools.base, "gcc",full.names=TRUE)[1]
            if (is.na(gcc)){
                gcc <- "";
            }
            for (x in rev(c(file.path(rtools.base, "bin")##,
                            ## file.path(rtools.base, "mingw_32/bin") ## Rtools sets up the mingw_32/bin first (even if x64)
                            ## file.path(rtools.base, ifelse(.Platform$r_arch == "i386","mingw_32/bin", "mingw_64/bin")),
                            ## file.path(rtools.base, ifelse(.Platform$r_arch == "i386","mingw_32/bin", "mingw_64/bin")),
                            ## file.path(rtools.base, ifelse(.Platform$r_arch == "i386","mingw_32/opt/bin", "mingw_64/opt/bin")),
                            ## ifelse(gcc == "", "", file.path(gcc, "bin")),
                            ## ifelse(gcc == "", "", ifelse(.Platform$r_arch == "i386",file.path(gcc, "bin32"), file.path(gcc, "bin64"))
                            ## )
                            ))){
                if (file.exists(x)){
                    path <- c(normalizePath(x), path);
                }
            }
            python.base <- rxPythonBaseWin();
            if (!is.null(python.base)){
                python <- normalizePath(file.path(python.base, "python.exe"));
                if (file.exists(python)){
                    Sys.setenv(PYTHON_EXE=python); ## For PythonInR
                    path <- c(normalizePath(python.base), path);
                    Sys.setenv(PYTHONHOME=python.base);
                    lib.path <- file.path(python.base, "Lib");
                    if (length(list.files(lib.path)) > 0){
                        Sys.setenv(PYTHONPATH=paste(python.base, normalizePath(lib.path), collapse=";"));
                    }
                }
            }
            ## java <- as.vector(Sys.which("java"));
            ## if (java != ""){
            ##     java <- sub(rex::rex(one_of("/", "\\"), except_any_of("/", "\\", "\n"), end), "", java)
            ## }
            keys <- NULL;
            ## Is there a 64 bit aspell that should be checked for...?
            keys <- try(utils::readRegistry("SOFTWARE\\Aspell", hive="HLM", view="32-bit", maxdepth=3), silent=TRUE);
            ## Add aspell for cran checks...
            if (!is.null(keys)){
                if (any(names(keys) == "Path")){
                    if (file.exists(keys$Path)){
                        path <- c(normalizePath(keys$Path), path);
                    }
                }
            }

            path <- path[path != ""];
            path <- paste(unique(path), collapse=";");
            Sys.setenv(PATH=path);
            return(TRUE);
        } else {
            return(FALSE)
        }
    }
}
##' Setup python and sympy for windows
##'
##' @author Matthew L. Fidler
##' @export
rxWinPythonSetup <- function(){
    base <- rxPythonBaseWin()
    if (is.null(base)){
        stop("This requires python.  Please setup and add to path.")
    }
    shell("pip install sympy")
    if (!requireNamespace("SnakeCharmR", quietly = TRUE)){
        devtools::install_github("asieira/SnakeCharmR");
    }
    message("To be safe, please restart R before using RxODE with SymPy")
}
##' Setup Windows components for RxODE
##'
##' @param rm.rtools Remove Rtools from path?
##' @author Matthew L. Fidler
##' @export
rxWinSetup <- function(rm.rtools=TRUE){
    if (!rxWinRtoolsPath(rm.rtools)){
        cat("This package will not work without Rtools being installed!\n");
        cat("Currently downloading https://cran.r-project.org/bin/windows/Rtools/Rtools33.exe...\n");
        rtools <- tempfile(fileext="-Rtools33.exe");
        cat("Downloading to", rtools, "\n");
        rxWget("http://cran.r-project.org/bin/windows/Rtools/Rtools33.exe", rtools);
        system(rtools);
        unlink(rtools);
        if (!rxWinRtoolsPath(rm.rtools)){
            stop("Rtools not setup correctly.");
        }
    }
    rxWinPythonSetup();
}
