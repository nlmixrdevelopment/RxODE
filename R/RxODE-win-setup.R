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
                return("");
            }
        }))
        path <- path[path != ""];
        if (rm.rtools){
            path <- path[regexpr(rex::rex(or("Rtools", "RTOOLS", "rtools")), path) == -1]
        }
        r.path <- normalizePath(file.path(Sys.getenv("R_HOME"),paste0("bin",Sys.getenv("R_ARCH"))));
        path <- c(r.path, path);
        if (file.exists("C:/Rtools")){
            gcc <- list.files("c:/Rtools", "gcc",full.names=TRUE)[1]
            if (is.na(gcc)){
                gcc <- "";
            }
            for (x in c("c:/Rtools/bin", ifelse(.Platform$r_arch == "i386","C:/Rtools/mingw_32/bin", "C:/Rtools/mingw_64/bin"),
                        ifelse(.Platform$r_arch == "i386","C:/Rtools/mingw_32/opt/bin", "C:/Rtools/mingw_64/opt/bin"),
                        ifelse(gcc == "", "", file.path(gcc, "bin")),
                        ifelse(gcc == "", "", ifelse(.Platform$r_arch == "i386",file.path(gcc, "bin32"), file.path(gcc, "bin64"))))){
                if (file.exists(x)){
                    path <- c(normalizePath(x), path);
                }
            }
            if (file.exists("C:/Python27/python.exe")){
                path <- c(normalizePath("c:/Python27"), path);
                Sys.setenv(PYTHONHOME="c:/Python27");
            }
            path <- path[path != ""];
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
        ## unlink("python-2.7.13.msi")
        if (!file.exists("python-2.7.13.msi")){
            rxWget("https://www.python.org/ftp/python/2.7.13/python-2.7.13.msi", "python-2.7.13.msi");
        }
        cat("Install python to the default location (c:/Python27)\n");
        shell("start python-2.7.13.msi")
        readline(prompt="Press [enter] to continue");
        if (!file.exists("C:/Python27/lib")){
            stop("Python installation unsuccessful.");
        } else {
            unlink("python-2.7.13.msi");
        }
    }
    if (!requireNamespace("SnakeCharmR", quietly = TRUE)){
        devtools::install_github("asieira/SnakeCharmR");
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
