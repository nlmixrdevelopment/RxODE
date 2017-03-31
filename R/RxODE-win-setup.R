rxWget <- function(url, to){
    if (.Platform$GUI == "RTerm"){
        if (Sys.which("wget") != ""){
            download.file(url, to, method="wget", extra="--progress=dot --no-check-certificate");
        } else if (Sys.which("curl" != "")){
            download.file(url, to, method="curl");
        } else {
            download.file(url, to);
        }
    } else{
        download.file(url, to);
    }
}
##' Setup Rtools path
##'
##' @param rm.rtools Remove the Rtools from the current path specs.
##'
##' @author Matthew L. Fidler
rxWinRtoolsPath <- function(rm.rtools=TRUE){
    if(.Platform$OS.type == "unix"){
        return(TRUE)
    } else {
        path <- unique(sapply(gsub("/", "\\\\", strsplit(Sys.getenv("PATH"), ";")[[1]]), function(x){
            if (file.exists(x)){
                return(normalizePath(x));
            } else {
                return(x);
            }
        }))
        if (rm.rtools){
            path <- path[regexpr(rex::rex(or("Rtools", "RTOOLS", "rtools")), path) == -1]
        }
        r.path <- normalizePath(file.path(Sys.getenv("R_HOME"),paste0("bin",Sys.getenv("R_ARCH"))));
        path <- c(r.path, path);
        if (file.exists("C:/Rtools")){
            gcc <- list.files("c:/Rtools", "gcc",full.names=TRUE)[1]
            for (x in c("c:/Rtools/bin", ifelse(.Platform$r_arch == "i386","C:/Rtools/mingw_32/bin", "C:/Rtools/mingw_64/bin"),
                        ifelse(.Platform$r_arch == "i386","C:/Rtools/mingw_32/opt/bin", "C:/Rtools/mingw_64/opt/bin"),
                        file.path(gcc, "bin"),
                        ifelse(.Platform$r_arch == "i386",file.path(gcc, "bin32"), file.path(gcc, "bin64")))){
                path <- c(normalizePath(x), path);
            }
            if (file.exists("C:/Python27/python.exe")){
                path <- c(normalizePath("c:/Python27"), path);
            }
            path <- paste(path, collapse=";");
            Sys.setenv(PATH=unique(path));
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
    if (!file.exists("C:/Python27/lib")){
        python <- tempfile(fileext="-python-2.7.13.msi");
        rxWget("https://www.python.org/ftp/python/2.7.13/python-2.7.13.msi", python);
        cat("Install python to the default location (c:/Python27)\n");
        system(sprintf("explorer \"%s\"",python))
        readline(prompt="Press [enter] to continue");
        unlink(python)
        if (!file.exists("C:/Python27/lib")){
            stop("Python installation unsuccessful.");
        }
    }
    if (!file.exists("C:/Python27/Lib/site-packages/sympy")){
        system("c:/Python27/Scripts/pip install sympy");
    }
    if (!requireNamespace("SnakeCharmR", quietly = TRUE)){
        devtools::install_github("mattfidler/SnakeCharmR");
    }
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
