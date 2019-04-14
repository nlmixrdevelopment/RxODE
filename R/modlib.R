.rxPkgInst <- function(obj){
    .wd <- getwd()
    if (regexpr(obj$package, .wd) != -1){
        .inst <- gsub(paste0("(", obj$package, ").*"), "\\1", .wd);
    } else {
        .inst <- system.file(package=obj$package)
    }
    if (regexpr("inst$", .inst) != -1) return(.inst)
    .inst2 <- file.path(.inst, "inst")
    if (file.exists(.inst2)) return(.inst2)
    .html <- file.path(.inst, "html");
    if (file.exists(.html)) return(.inst);
    return(.inst2)
}

.rxPkgDir <- function(obj){
    return(file.path(.rxPkgInst(obj), "rx"));
}

.rxPkgDll <- function(obj){
    obj$mdir <- .rxPkgDir(obj);
    file.path(obj$mdir, basename(obj$rxDll$dll))
}

##' Use model object in your package
##' @param obj model to save.
##' @export
rxUse <- function(obj, internal = FALSE, overwrite = TRUE, compress = "bzip2"){
    rxReq("usethis")
    rxReq("devtools")
    if (missing(obj)){
        .env <- new.env();
        assign("internal", internal, .env)
        assign("overwrite", overwrite, .env)
        assign("compress", compress, .env)
        for (.f in list.files(path=devtools::package_file("data"),pattern="\\.rda$",full.names=TRUE)){
            load(.f, envir=.env)
            .f2 <- basename(.f);
            .f2 <- substr(.f2, 0, nchar(.f2) - 4)
            if (is(.env[[.f2]], "RxODE")){
                message(sprintf("Recompile %s (if needed)", .f2))
                eval(parse(text=sprintf("rxUse(%s, internal=internal, overwrite=overwrite, compress=compress)", .f2)),
                     envir=.env)
            }
        }
        return(invisible(TRUE))
    } else {
        .modName <- as.character(substitute(obj));
        .pkg <- basename(usethis::proj_get())
        .env <- new.env();
        assign(.modName, RxODE(rxNorm(obj), package=.pkg, modName=.modName), .env);
        assign("internal", internal, .env)
        assign("overwrite", overwrite, .env)
        assign("compress", compress, .env)
        eval(parse(text=sprintf("usethis::use_data(%s, internal=internal, overwrite=overwrite, compress=compress)", .modName)),
             envir=.env)

    }
}


