rxWget <- function(url, to){
  cat("Checking for wget.exe...\n")
  
  if (Sys.which("wget") == ""){
        if (.Platform$r_arch == "i386"){
          
          rtoolslist <- paste(letters, ":/RTOOLS/mingw_32/bin", sep="")
          rtoolspath <- rtoolslist[which(file.exists(rtoolslist, sep=""))]
          
            if (any(file.exists(rtoolspath))){
              try(download.file("https://eternallybored.org/misc/wget/current/wget.exe", paste(rtoolspath, "/wget.exe", sep="")))
              if (!file.exists(paste(rtoolspath, "/wget.exe", sep=""))){
                    stop(paste("Cannot install wget. Please download from https://eternallybored.org/misc/wget/current/wget.exe and install into ",
                               rtoolspath, " before continuing.\n", sep=""))
                }
            } else {
                try(download.file("https://eternallybored.org/misc/wget/current/wget.exe", "wget.exe"))

            }
        } else {
          
          rtoolslist <- paste(letters, ":/RTOOLS/mingw_64/bin", sep="")
          rtoolspath <- rtoolslist[which(file.exists(rtoolslist, sep=""))]
          
          if (any(file.exists(rtoolspath))){
                try(download.file("https://eternallybored.org/misc/wget/current/wget64.exe", paste(rtoolspath, "/wget.exe", sep="")))
                if (!file.exists(paste(rtoolspath, "/wget.exe", sep=""))){
                  stop(paste("Cannot install wget. Please download from https://eternallybored.org/misc/wget/current/wget64.exe and install as wget.exe into ",
                             rtoolspath, " before continuing.\n", sep=""))
                }
            } else {
                try(download.file("https://eternallybored.org/misc/wget/current/wget64.exe", "wget.exe"))
            }
        }
    }
    if (Sys.which("wget") == ""){
        stop("wget is not available. Please install manually before continuing.");
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
        
        if (length(grep("rtools", tolower(Sys.which("gcc.exe"))))==0) {
          stop("RxODE cannot be installed, since Rtools isn't set up appropriately. Please (re)install it and try again.\n")
        } 
        
        # if (!file.exists(rtools.base)){
        #     keys <- NULL
        #     try(keys <- utils::readRegistry("SOFTWARE\\R-core\\Rtools", hive = "HCU", view = "32-bit", maxdepth = 2), silent = TRUE)
        #     if (is.null(keys) || length(keys) == 0)
        #         try(keys <- utils::readRegistry("SOFTWARE\\R-core\\Rtools", hive = "HLM", view = "32-bit", maxdepth = 2), silent = TRUE)
        #     if (is.null(keys) || length(keys) == 0){
        #         stop("Cannot use this package because Rtools isn't setup appropriately...")
        #     }
        # 
            # for(i in seq_along(keys)) {
            #     version <- names(keys)[[i]]
            #     key <- keys[[version]]
            #     if (!is.list(key) || is.null(key$InstallPath)) next;
            #     install_path <- normalizePath(key$InstallPath, mustWork = FALSE, winslash = "/");
            #     if (file.exists(install_path)){
            #         rtools.base <- install_path;
            #     }
            # }
        # }
      
      rtoolslist  <- paste(letters, ":/Rtools", sep="")
      rtools.base <- rtoolslist[which(file.exists(rtoolslist, sep=""))]
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
            for (x in rev(c(file.path(rtools.base, "bin"),
                            ## file.path(rtools.base, "mingw_32/bin") ## Rtools sets up the mingw_32/bin first (even if x64)
                            file.path(rtools.base, ifelse(.Platform$r_arch == "i386","mingw_32/bin", "mingw_64/bin")),
                            file.path(rtools.base, ifelse(.Platform$r_arch == "i386","mingw_32/bin", "mingw_64/bin")),
                            file.path(rtools.base, ifelse(.Platform$r_arch == "i386","mingw_32/opt/bin", "mingw_64/opt/bin"))
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
                        Sys.setenv(PYTHONPATH=paste(python.base, normalizePath(lib.path), sep=";"));
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
        stop("RxODE requires Python. Please install an appropriate version and add it to the system path.")
    }
    if (file.access(paste(base, "/Lib/site-packages", sep=""),2)==-1){
      stop("The Python library path does not appear to be writeable. Please rectify this situation, restart R, and try again.")
    } 
    message("Attempting to install simpy. This may take a few seconds...")
    try(system("python -m pip install sympy"))
    
    if (!requireNamespace("SnakeCharmR", quietly = TRUE)){
        message("Attempting to install SnakeCharmR. This may take a few seconds...")
        devtools::install_github("asieira/SnakeCharmR");
    }
    message("Please restart your R session before using RxODE with sympy.")
}
##' Setup Windows components for RxODE
##'
##' @param rm.rtools Remove Rtools from path?
##' @author Matthew L. Fidler
##' @export
rxWinSetup <- function(rm.rtools=TRUE){
    if (!rxWinRtoolsPath(rm.rtools)){
        cat("This package requires Rtools! Please download from http://cran.r-project.org/bin/windows/Rtools/, install and restart your R session before proceeding.\n");
        # cat("Currently downloading https://cran.r-project.org/bin/windows/Rtools/Rtools33.exe...\n");
        # rtools <- tempfile(fileext="-Rtools33.exe");
        # cat("Downloading to", rtools, "\n");
        # rxWget("http://cran.r-project.org/bin/windows/Rtools/Rtools33.exe", rtools);
        # system(rtools);
        # unlink(rtools);
        # if (!rxWinRtoolsPath(rm.rtools)){
        #     stop("Rtools not setup correctly.");
        # }
    }
    rxWinPythonSetup();
}
