library(devtools)
library(tidyr)
library(dplyr)
cat("Copy header to inst directory")

file.copy(devtools::package_file("src/RxODE_types.h"),
          devtools::package_file("inst/include/RxODE_types.h"),
          overwrite=TRUE);

cat("Generate header string.\n");
odec <- readLines(devtools::package_file("inst/ode.c"));
w <- which(regexpr("__ODE_SOLVER__", odec) != -1)[1];
ode <- odec[seq(1, w - 1)];
solve <- odec[seq(w, length(odec))];
hd <- sprintf("#define __HD_ODE__ \"%s\\n\"\n#define __HD_SOLVE__ \"%s\\n\"\n",
              paste(gsub("%", "%%", gsub("\"", "\\\\\"", ode)), collapse="\\n"),
              paste(gsub("%", "%%", gsub("\"", "\\\\\"", solve)), collapse="\\n"));
sink(devtools::package_file("src/ode.h"))
cat(hd);
sink();
cat("Generate dplyr and tidyr compatability functions.\n")
sink(package_file("R/rxsolve-gen.R"))
cat("## Generated code from build/refresh.R\n\n");
for (f in c("spread_", "unite_", "separate_", "gather_")){
    fn <- deparse(eval(parse(text=sprintf("args(%s)", f))));
    fn <- paste0(fn[-length(fn)], collapse="\n");
    one <- eval(parse(text=sprintf("attr(formals(%s)[1],\"names\")", f)));
    cat(sprintf("##' @importFrom tidyr %s
##' @export
%s.solveRxDll <- %s{
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$%s <- dplyr::as.tbl(%s)
  return(do.call(getFromNamespace(\"%s\",\"dplyr\"), call, envir = parent.frame(1)));
}\n\n",f, f, fn, one, one, f));
}

## merge?
for (f in c("sample_frac", "sample_n", "group_by_", "rename_", "arrange_",
            "summarise_", "transmute_", "mutate_", "distinct_", "rename_",
            "select_", "arrange_", "slice_", "filter_")){
    fn <- deparse(eval(parse(text=sprintf("args(%s)", f))));
    fn <- paste0(fn[-length(fn)], collapse="\n");
    one <- eval(parse(text=sprintf("attr(formals(%s)[1],\"names\")", f)))
    cat(sprintf("##' @importFrom dplyr %s
##' @export
%s.solveRxDll <- %s{
  call <- as.list(match.call())[-1];
  call$%s <- %s(%s);
  return(do.call(getFromNamespace(\"%s\",\"dplyr\"), call, envir = parent.frame(1)));
}\n\n",f, f, fn, one, ifelse(one == "tbl", "dplyr::as.tbl", "asTbl"), one, f));
}

for (f in c("row.names", "by", "aggregate", "anyDuplicated", "droplevels",
            "duplicated", "edit", "is.na", "Math", "rowsum",
            "split", "subset", "stack", "unstack", "unique",
            "within", "with")){
    if (f == "Math"){
        fn <- "function(x, ...)";
    } else {
        fn <- deparse(eval(parse(text=sprintf("args(%s)", f))));
        fn <- paste0(fn[-length(fn)], collapse="\n");
    }
    one <- eval(parse(text=sprintf("attr(formals(%s)[1],\"names\")", f)))
    ns <- gsub("package:", "", find(f));
    if (!any(ns == c("base"))){
        cat(sprintf("##' @importFrom %s %s\n", ns, f));
    }
    cat(sprintf("##' @export
%s.solveRxDll <- %s{
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$%s <- as.data.frame(%s)
  return(do.call(\"%s.data.frame\", call, envir = parent.frame(1)));
}\n\n",f, fn, one, one, f));
}

for (f in c("dim", "dimnames", "t")){
    fn <- deparse(eval(parse(text=sprintf("args(%s)", f))));
    fn <- paste0(fn[-length(fn)], collapse="\n");
    one <- "x"
    cat(sprintf("##' @export
%s.solveRxDll <- %s{
  call <- as.list(match.call(expand.dots=TRUE))[-1];
  call$%s <- as.matrix(%s)
  return(do.call(\"%s\", call, envir = parent.frame(1)));
}\n\n",f, fn, one, one, f));
}
sink();

cat("Update Parser c file\n");
dparser::mkdparse(devtools::package_file("inst/tran.g"),
         devtools::package_file("src/"),
         grammar_ident="RxODE");
file <- gsub("^([#]line [0-9]+ )\".*(src)/+(.*)\"","\\1\"\\2/\\3\"",
             readLines(devtools::package_file("src/tran.g.d_parser.c")))
sink(devtools::package_file("src/tran.g.d_parser.c"))
cat(paste(file,collapse="\n"));
sink();
unlink(devtools::package_file("src/tran.o"))
sink(devtools::package_file("R/version.R"))
cat("##\' Version and repository for this dparser package.
##\'
##\' @return A character vector with the version and repository.
##\' @author Matthew L. Fidler
##\' @keywords internal
##\' @export
rxVersion <- function(){return(c(version=\"");
ver <- readLines(devtools::package_file("DESCRIPTION"));
ver <- ver[regexpr("^Version:", ver) != -1]
ver <- ver[1];
cat(gsub("Version: +", "", ver))
cat("\",repo=\"");
cat("https://github.com/")
tmp <- readLines(devtools::package_file(".git/config"))
cat(gsub("\\.git$", "", gsub(".*git@github.com:", "", tmp[which(tmp == '[remote "origin"]')[1]+1])))
cat("\"))}\n");
sink();
devtools::load_all();

document();

cat("Update README\n");
owd <- getwd();
on.exit({setwd(owd)});
setwd(devtools::package_file());
knitr::knit(devtools::package_file("README.Rmd"))


create.syminv.cache <- function(n=1, xforms=c("sqrt", "log", "identity")){
    owd <- getwd();
    on.exit(setwd(owd))
    setwd(devtools::package_file("inst/inv"));
    i <- 0;
    if (n == 1){
        m <- matrix(1, 1);
        for (xform in c("sqrt", "log", "identity")){
            rxSymInvC(m, xform);
        }
    } else {
        m <- matrix(rep(1, n * n), n);
        m[lower.tri(m)] <- 0;
        m[upper.tri(m)] <- 0;
        n2 <- sum(lower.tri(m) * 1);
        grd <- eval(parse(text=sprintf("expand.grid(%s)", paste(rep("c(0,1)", n2), collapse=", "))));
        for (i in 1:length(grd[, 1])){
            v <- as.vector(m);
            vw <- which(as.vector(lower.tri(m) * 1) == 1);
            v[vw] <- grd[i, ];
            mc <- matrix(v, n);
            mc[upper.tri(mc)] <- t(mc)[upper.tri(mc)];
            for (xform in xforms){
                rxSymInvCreate(mc, xform);
                cat(".");
                i <- i + 1;
                if (i %% 5 == 0){
                    cat(i);
                }
                if (i %% 50 == 0){
                    cat("\n");
                }
            }
            ## print(mc);
            ## mc[upper.tri(mc)] <- grd[i, ]
            ## mc[lower.tri(mc)] <- t(mc)[lower.tri(mc)];
            ## print(mc);
        }
        ## print(n2)
        ## print(m);
        ## print(grd);
    }
    cat("\n");
}

create.syminv.cache();

create.syminv.cache(2);

create.syminv.cache(3);

## create.syminv.cache(4);
