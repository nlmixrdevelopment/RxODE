.pkgModelCurrent  <- TRUE

.setPkgModels <- function(value){ ## For testing
    assignInMyNamespace(".pkgModelCurrent",value);
}

.rxPkgInst <- function(obj){
    .wd <- getwd()
    if (regexpr(obj$package, .wd) != -1){
        .inst <- gsub(paste0("(", obj$package, ").*"), "\\1", .wd);
    } else {
        .inst <- system.file(package=obj$package)
    }
    if (assertthat::is.writeable(.inst)){
        if (regexpr("inst$", .inst) != -1) return(.inst)
        .inst2 <- file.path(.inst, "inst")
        if (file.exists(.inst2)) return(.inst2)
        .html <- file.path(.inst, "html");
        if (file.exists(.html)) return(.inst);
        return(.inst2)
    } else {
        .inst  <- "~/.rxCache/";
        if (assertthat::is.writeable(.inst)){
            return(.inst)
        }
        return(rxTempDir())
    }
}

.rxPkgDir <- function(obj){
    return(file.path(.rxPkgInst(obj), "rx"));
}

.rxPkgDll <- function(obj){
    obj$mdir <- .rxPkgDir(obj);
    .pkgInfo  <- getLoadedDLLs()[[obj$package]]
    if (!all(is.null(.pkgInfo))){
        if (obj$isValid()){
            .tmp <- .pkgInfo
            class(.tmp)  <- "list"
            return(.tmp$path)
        } else {
            return(file.path(obj$mdir, basename(obj$rxDll$dll)))
        }
    } else {
        return(file.path(obj$mdir, basename(obj$rxDll$dll)))
    }
}

.rxNewMvStr  <- function(obj){
    gsub("[.].*",paste0("_new_",.Platform$r_arch,"_model_vars"),basename(obj$rxDll$dll))
}

.rxPkgLoaded  <- function(pkg){
    .si  <- sessionInfo();
    return(length(intersect(pkg, c(names(.si$otherPkgs)## ,names(.si$loadedOnly)
                                   ))) != 0)
}

##' Use model object in your package
##' @param obj model to save.
##' @param internal If this is run internally.  By default this is FALSE
##' @inheritParams usethis::use_data
##' @export
rxUse <- function(obj, overwrite = TRUE, compress = "bzip2",
                  internal=FALSE){
    rxReq("usethis")
    rxReq("devtools")
    internal  <-  internal;
    if (missing(obj)){
        .env <- new.env();
        assign("internal", internal, .env)
        assign("overwrite", overwrite, .env)
        assign("compress", compress, .env)
        sapply(list.files(devtools::package_file("inst/rx"), full.names=TRUE),
               unlink,force=TRUE,recursive=TRUE)
        .models <- c();
        for (.f in list.files(path=devtools::package_file("data"),
                              pattern="\\.rda$",full.names=TRUE)){
            load(.f, envir=.env)
            .f2 <- basename(.f);
            .f2 <- substr(.f2, 0, nchar(.f2) - 4)
            if (is(.env[[.f2]], "RxODE")){
                .env[[.f2]]$package <- NULL;
                message(sprintf("Recompile %s (if needed)", .f2))
                .models  <- c(.models,.f2);
                eval(parse(text=sprintf("rxUse(%s, internal=internal, overwrite=overwrite, compress=compress)", .f2)),
                     envir=.env)
                .docFile  <- file.path(devtools::package_file("R"),paste0(.f2,"-doc.R"));
                if (!file.exists(.docFile)){
                    message(sprintf("Creating documentation file %s", .docFile));
                    sink(.docFile);
                    .tmp <- .env[[.f2]];
                    .mv  <- rxModelVars(.tmp);
                    cat(sprintf("##' %s RxODE model\n",.f2));
                    cat("##'\n");
                    cat(sprintf("##' @format An \\emph{RxODE} model with %s parameters, %s ODE states, and %s calc vars.\n",
                                length(.tmp$params), length(.tmp$state)+.mv$extraCmt, length(.tmp$lhs)))
                    cat("##'\n")
                    cat(sprintf("##'\\emph{Parameters (%s$params)}\n",.f2))
                    cat("##'\n")
                    cat("##' \\describe{\n")
                    .def  <- rxInits(.tmp);
                    .defs  <- paste0(" (default=",.def,")");
                    .defs[is.na(.def)]  <- "";
                    cat(paste(paste0("##'   \\item{",.tmp$params,"}{",.defs,"}\n"),collapse=""))
                    cat("##'}\n");
                    .state  <- .tmp$state
                    ##
                    if (.mv$extraCmt==2){
                        .state  <- c(.state, "depot", "central");
                    } else if (.mv$extraCmt==1){
                        .state  <- c(.state, "central");
                    }
                    if (length(.state)>0){
                        cat("##'\n")
                        cat(sprintf("##' \\emph{State %s$state}\n",.f2))
                        cat("##'\n")
                        cat("##' \\describe{\n")
                        cat(paste(paste0("##'   \\item{",.state,"}{ (=",seq_along(.state),")}\n"),collapse=""))
                        cat("##' }\n")
                    }
                    .lhs  <- .tmp$lhs
                    if (length(.lhs)>0){
                        cat("##'\n")
                        cat(sprintf("##' \\emph{Calculated Variables %s$lhs}\n",.f2))
                        cat("##'\n")
                        cat("##' \\describe{\n")
                        cat(paste(paste0("##'   \\item{",.lhs,"}{}\n"),collapse=""));
                        cat("##' }\n")
                    }
                    cat("##'\n")
                    ## cat(sprintf("##' \\emph{Model Code}\n",.f2))
                    ## cat("##'\n")
                    ## .code  <- deparse(body(eval(parse(text=paste("function(){",rxNorm(.tmp),"}")))))
                    ## .code[1]  <- "RxODE({"
                    ## .code[length(.code)]  <- "})";
                    ## cat(paste(paste0("##' ",.code,"\n"),collapse=""));
                    ## cat("##'\n")
                    cat(paste(paste0("##' @seealso \\code{\\link[RxODE]{eventTable}}, \\code{\\link[RxODE]{et}}, \\code{\\link[RxODE]{rxSolve}}, \\code{\\link[RxODE]{RxODE}}\n")))
                    cat("##' \n");
                    cat("##' @examples\n");
                    cat("##' ## Showing the model code\n");
                    cat(sprintf("##' summary(%s)\n", .f2));
                    cat("##'\n");
                    cat(sprintf('"%s"\n',.f2));
                    sink();
                }
            }
        }
        if (!dir.exists(devtools::package_file("src"))){
            dir.create(devtools::package_file("src"),recursive=TRUE);
        }
        .pkg <- basename(usethis::proj_get())
        sapply(list.files(devtools::package_file("inst/rx"),pattern="[.]c",full.names=TRUE),
               function(x){
            message(sprintf("copy %s", basename(x)))
            .f0  <- readLines(x);
            if (.pkg=="RxODE"){
                .f0[1]  <- "#include \"../inst/include/RxODE.h\"\n#include \"../inst/include/RxODE_model.h\"";
            } else {
                .f0[1]  <- "#include <RxODE.h>\n#include <RxODE_model.h>";
            }
            writeLines(text=.f0, con=file.path(devtools::package_file("src"), basename(x)));
        })
        rxClean(devtools::package_file("inst/rx"))
        .inits  <- paste0("R_init0_",.pkg,"_",.models);
        .tmp <- paste0("{\"",.pkg,"_",.models,"_model_vars\", (DL_FUNC) &",.pkg,"_",.models,"_model_vars, 0},\\");
        .tmp[length(.tmp)]  <- substr(.tmp[length(.tmp)], 0, nchar(.tmp[length(.tmp)])-1)
        .extraC  <- c("#define compiledModelCall \\",
                      .tmp,
                      paste0("SEXP ",.pkg, "_",.models,"_model_vars();"),
                      paste0("void ",.inits, "();"),
                      paste0("void R_init0_",.pkg,"_RxODE_models(){"),
                      paste0("  ",.inits, "();"),
                      "}")
        sink(file.path(devtools::package_file("src"),paste0(.pkg, "_compiled.h")))
        if (.pkg=="RxODE"){
            cat("#include <R.h>\n#include <Rinternals.h>\n#include <stdlib.h> // for NULL\n#include <R_ext/Rdynload.h>\n#include \"../inst/include/RxODE.h\"\n#include \"../inst/include/RxODE_model_shared.h\"\n")
        } else {
            cat("#include <R.h>\n#include <Rinternals.h>\n#include <stdlib.h> // for NULL\n#include <R_ext/Rdynload.h>\n#include <RxODE.h>\n#include <RxODE_model_shared.h>\n")
        }
        cat(paste(.extraC,collapse="\n"))
        cat("\n")
        sink()
        .files  <- list.files(devtools::package_file("src"));
        if (all(regexpr(paste0("^",.pkg),.files) !=-1)){
            message(sprintf("Only compiled models in this package, creating %s_init.c", .pkg));
            sink(file.path(devtools::package_file("src"),paste0(.pkg, "_init.c")))
            cat("#include <R.h>\n#include <Rinternals.h>\n#include <stdlib.h> // for NULL\n#include <R_ext/Rdynload.h>\n")
            cat('#include <RxODE.h>\n');
            cat("#include <RxODE_model.h>\n");
            cat(paste0('#include "',.pkg,'_compiled.h"\n'))
            cat(sprintf('void R_init_%s(DllInfo *info){\n', .pkg));
            cat(sprintf("  R_init0_%s_RxODE_models();\n", .pkg));
            cat("  static const R_CallMethodDef callMethods[]  = {\n  compiledModelCall\n  {NULL, NULL, 0}\n  };\n");
            cat("  R_registerRoutines(info, NULL, callMethods, NULL, NULL);\n")
            cat("  R_useDynamicSymbols(info,FALSE);\n")
            cat('}\n');
            cat(paste(paste0("void R_unload_",.pkg, "_", .models, "(DllInfo *info);\n"),collapse=""));
            cat(sprintf('void R_unload_%s(DllInfo *info){\n', .pkg));
            cat(paste(paste0("  R_unload_", .pkg,"_", .models ,"(info);\n"),collapse=""));
            cat('}\n');
            sink()
        }
        if (!file.exists(devtools::package_file("R/rxUpdated.R")) && .pkg !="RxODE"){
            sink(devtools::package_file("R/rxUpdated.R"));
            cat(".rxUpdated <- new.env(parent=emptyenv())\n");
            sink();
        }
        return(invisible(TRUE))
    } else {
        .modName <- as.character(substitute(obj));
        .pkg <- basename(usethis::proj_get())
        .env <- new.env(parent=baseenv());
        ## if (.pkg=="RxODE"){
        ##     ## Don't recompile, just update internals
        ##     obj$package <- "RxODE";
        ##     obj$modName <- paste0("RxODE_",.modName);
        ##     obj$mdir  <- devtools::package_file("inst/rx");
        ##     .updateRxModelLib(obj);
        ##     .dll <- obj$rxDll;
        ##     .mv <- rxUpdateTrans_(.dll$modVars, paste0("RxODE_",.modName,"_"),
        ##                           "RxODE");
        ##     .dll$modVars <- .mv
        ##     obj$rxDll <- .dll
        ##     rxCompile(.mv,
        ##               dir=devtools::package_file("inst/rx"),
        ##               prefix=paste0("RxODE_",.modName,"_"),
        ##               extraC = NULL,
        ##               debug = NULL,
        ##               modName = paste("RxODE_",.modName),
        ##               package="RxODE");
        ##     assign(.modName, obj, .env);
        ##     assign("internal", internal, .env)
        ##     assign("overwrite", overwrite, .env)
        ##     assign("compress", compress, .env)
        ##     eval(parse(text=sprintf("usethis::use_data(%s, internal=internal, overwrite=overwrite, compress=compress)", .modName)),
        ##          envir=.env)
        ## }
        ## else {
        obj$package <- NULL;
        assign(.modName, RxODE(rxNorm(obj), package=.pkg, modName=.modName), .env);
        assign("internal", internal, .env)
        assign("overwrite", overwrite, .env)
        assign("compress", compress, .env)
        eval(parse(text=sprintf("usethis::use_data(%s, internal=internal, overwrite=overwrite, compress=compress)", .modName)), envir=.env);
        ## }

    }
}
