.rxIndLinLine <- function(line, states, state0){
    .tmp <- symengine::expand(line) ## Expand line
    .tmp <- rxFromSE(.tmp); ## Convert SE->RxODE; Changes things like X^1 -> X
    .ret <- eval(parse(text=paste0("rxSplitPlusQ(quote(", .tmp, "))")));
    .lst <- list()
    .fullIndLin <- FALSE
    ## This adds the x to the name environment in the .lst list above.
    .addIt <- function(x, name){
        .cur <- .lst
        if (any(names(.lst) == name)){
            .cur[[name]] <- c(.cur[[name]], x)
        } else {
            .cur[[name]] <- x;
        }
        .lst <<- .cur;
    }
    ## This fixes the multiplication so that -1*x becomes -x and *1* becomes *
    .fixMult <- function(.mult){
        .mult <- .mult[.mult != "1"];
        .w <- which(.mult == "-1");
        if (length(.w) > 0){
            .mult[.w + 1] <- paste0("-", .mult[.w + 1]);
            .mult <- .mult[-.w];
        }
        return(.mult);
    }
    ## This collapses a character vector with "*" between it
    ## If we have *1/x this becomes simply /x
    .multCollapse <- function(x){
        gsub("[*]1[/] *", "/", paste(.fixMult(x), collapse="*"))
    }
    sapply(.ret, function(x){
        .mult <- eval(parse(text=paste0("rxSplitPlusQ(quote(", x, "),mult=TRUE)")))
        ## Handle -state by changing "-state" to "-1",  "state"
        if (any(.mult[1] == paste0("-", states))){
            .mult <- c("-1", substr(.mult[1], 2, nchar(.mult[1])), .mult[-1])
        }
        .num <- which(sapply(.mult, function(x){any(x == states)}))

        ## For anything with "t" in it it has to be part of the
        ## inductive part, not part of the matrix exponential.
        ## I'm not sure it is even makes sense mathematically.
        if (any(.mult == "t") || any(.mult == "-t")) .num <- -1L
        if (length(.num) == 1){
            ## This implies that this term is something * state;
            ## However the "something" may include other state items. Check that here.
            .state <- .mult[.num]
            .rest <- .mult[-.num]
            .dep <- any(sapply(.rest,
                               function(x){
                .pars <- rxModelVars(paste0("rx_expr=", x))$params
                any(sapply(.pars, function(y){any(y == states)}))
            }))
            if (.dep){
                .fullIndLin <<- TRUE
                .expr <- .multCollapse(.mult)
                .addIt(.expr, "_rxF");
            } else {
                if (length(.rest) == 0){
                    .addIt("1", .state);
                } else if (length(.rest) == 1 && .rest[1] == "-1"){
                    .addIt("-1", .state);
                } else {
                    .expr <- .multCollapse(.rest)
                    .addIt(.expr, .state);
                }
            }
        } else {
            ## Something else is here.
            .expr <- .multCollapse(.mult);
            .addIt(.expr, "_rxF")
        }
    })
    .ret <- sapply(c(states, "_rxF"), function(x){
        if (any(names(.lst) == x)) return(gsub("[+][-]", "-", paste(.lst[[x]], collapse="+")))
        return("0");
    });
    ## Check for an independet solution ie -ka*depot
    .tmp2 <- .ret[!(names(.ret) %in% c(state0, "_rxF"))];
    ## FIXME put infusions in EVERYWHERE for _rxF, that way when dropping it will get the right number.
    if (all(.tmp2 == "0")){
        ## independent solution
        .eAt <- sprintf("exp(%s*t)", .ret[state0]);
        if (.ret["_rxF"] == "0"){
            .ind <- sprintf("%s*(%s)+(%s-1)/(%s)*(rxInf[%d])", .eAt, state0, .eAt, .ret[state0], which(state0 == states) - 1);
        } else if (.ret[state0] != "0"){
            .ind <- sprintf("%s*(%s)+(%s-1)/(%s)*(%s+rxInf[%d])", .eAt, state0, .eAt, .ret[state0], .ret["_rxF"], which(state0 == states) - 1);
        } else {
            stop("Cannot calculate this by inductive linearization?");
        }
    } else {
        .ind <- ""
    }
    return(c(.ret, paste0(.fullIndLin), .ind))
}

##' This creates the inductive linearization pieces to integrate into
##' a RxODE model.
##'
##' @param model RxODE model type of object
##'
##' @param doConst Replace constants with values; By default this is
##'     \code{FALSE}.
##'
##' @return List:
##' \itemize{
##'
##' \item Matrix Exponential initial matrix A for exp(t*A)
##'
##' \item Inductive Linerization vector for F
##'
##' \item Extra RxODE code for model generation
##'
##' \item Generated C code for model variables; With ME only this will
##' be a list of size 1, otherwise it is a list of size 2.
##'
##' }
##' @author Matthew Fidler
##' @noRd
.rxIndLin <- function(model, doConst=FALSE){
    .env <- .rxLoadPrune(model, doConst=doConst)
    .states <- rxState(.env);
    .ret <- eval(parse(text=rxIndLin_(.states)))
    .fullIndLin <- any(.ret[, "indLin"] == "TRUE")
    if (.fullIndLin){
        .w <- which(.ret[, "ind"] != "")
        if (length(.w) > 0){
            .ind <- .ret[, "ind"];
            .w2 <- which(.ind == "");
            .w3 <- which(.ind != "");
            if (length(.w2) > 0){
                .ind[.w2] <- .states[.w2];
            }
            for (.i in .w3){
                .tmp <- paste0(.ret[,.states[.i]], "*", .states[.i], "+")
                .tmp[which(.ret[, .states[.i]] == "0")] <- "";
                .ret[, "_rxF"] <- gsub("[+][-]", "-", paste0(.tmp, .ret[, "_rxF"]));
            }
            .ret <- .ret[-.w3, -.w3, drop = FALSE];
            .w2 <- setNames(.w2 - 1, NULL);
            .oldStates <- .states
            .code <- paste(paste0("_rxI=", .ind, ";"), collapse="\n");
            .states <- .states[.states %in% dimnames(.ret)[[1]]];
        } else {
            .w2 <- integer(0);
            .oldStates <- .states
            .code <- NULL;
        }
    } else {
        .w2 <- integer(0);
        .oldStates <- .states
        .code <- NULL;
    }
    .ret0 <- .ret[.states, .states, drop = FALSE]
    .ret1 <- .ret[, "_rxF", drop = FALSE]
    .code <- c(.code, paste0("_rxM=", as.vector(.ret0), ";"));
    if (!all(.ret1 == "0")){
        .code <- c(.code, paste("_rxF=", as.vector(.ret1)))
        .codeSave <- c("SEXP matLst = PROTECT(allocVector(VECSXP, 3));pro++;",
                       "SEXP matNme = PROTECT(allocVector(STRSXP, 3));pro++;",
                       "SET_STRING_ELT(matNme, 0, mkChar(\"A\"));",
                       "SET_STRING_ELT(matNme, 1, mkChar(\"f\"));",
                       "SET_STRING_ELT(matNme, 2, mkChar(\"fullIndLin\"));",
                       "setAttrib(matLst, R_NamesSymbol, matNme);",
                       ## mat0 setup the nstate x nstate matrix for
                       ## inductive linearization
                       sprintf("SEXP mat0 = PROTECT(allocMatrix(STRSXP, %s, %s));pro++;", dim(.ret0)[1], dim(.ret0)[2]),
                       paste(paste0("SET_STRING_ELT(mat0, ", seq_along(as.vector(.ret0)) - 1, ",mkChar(\"", as.vector(.ret0), "\"));"), collapse="\n"),
                       sprintf("SEXP mat0nn = PROTECT(allocVector(STRSXP, %s));pro++;", dim(.ret0)[[1]]),
                       paste(paste0("SET_STRING_ELT(mat0nn, ", seq_along(dimnames(.ret0)[[1]]) - 1,
                                    ",mkChar(\"", dimnames(.ret0)[[1]], "\"));"), collapse="\n"),
                       ## Assumes state is the SEXP that we can use for
                       ## state names.
                       "SEXP mat0n = PROTECT(allocVector(VECSXP, 2));pro++;",
                       "SET_VECTOR_ELT(mat0n, 0, mat0nn);",
                       "SET_VECTOR_ELT(mat0n, 1, mat0nn);",
                       "setAttrib(mat0, R_DimNamesSymbol, mat0n);",
                       ## Now setup matF or the F matrix from inductive linearization
                       sprintf("SEXP matF = PROTECT(allocMatrix(STRSXP, %s, %s));pro++;", dim(.ret1)[1], dim(.ret1)[2]),
                       paste(paste0("SET_STRING_ELT(matF, ", seq_along(as.vector(.ret1)) - 1, ",mkChar(\"", as.vector(.ret1), "\"));"), collapse="\n"),
                       "SEXP matFd = PROTECT(allocVector(INTSXP, 2));pro++;",
                       ## Assumes state is the SEXP that we can use for state names.
                       "SEXP matFn = PROTECT(allocVector(VECSXP, 2));pro++;",
                       "SET_VECTOR_ELT(matFn, 0, state);",
                       "SEXP matrxF = PROTECT(allocVector(STRSXP,1)); pro++;",
                       "SET_STRING_ELT(matrxF, 0, mkChar(\"_rxF\"));",
                       "SET_VECTOR_ELT(matFn, 1, matrxF);",
                       "setAttrib(matF,R_DimNamesSymbol, matFn);",
                       "SET_VECTOR_ELT(matLst, 0, mat0);",
                       "SET_VECTOR_ELT(matLst, 1, matF);",
                       "SEXP fullIndLin = PROTECT(allocVector(LGLSXP, 1)); pro++;",
                       paste0("LOGICAL(fullIndLin)[0] = ", ifelse(.fullIndLin, "1", "0"), ";"),
                       "SET_VECTOR_ELT(matLst, 2, fullIndLin);",
                       "SET_VECTOR_ELT(lst,  18, matLst);"
                       );
    } else {
        ## Only ME solution, others not needed
        .codeSave <- c("SEXP matLst = PROTECT(allocVector(VECSXP, 1));pro++;",
                       "SEXP matNme = PROTECT(allocVector(STRSXP, 1));pro++;",
                       "SET_STRING_ELT(matNme, 0, mkChar(\"A\"));",
                       "setAttrib(matLst, R_NamesSymbol, matNme);",
                       ## mat0 setup the nstate x nstate matrix for
                       ## inductive linearization
                       sprintf("SEXP mat0 = PROTECT(allocMatrix(STRSXP, %s, %s));pro++;", dim(.ret0)[1], dim(.ret0)[2]),
                       paste(paste0("SET_STRING_ELT(mat0, ", seq_along(as.vector(.ret0)) - 1, ",mkChar(\"", as.vector(.ret0), "\"));"), collapse="\n"),
                       ## Assumes state is the SEXP that we can use for
                       ## state names.
                       "SEXP mat0n = PROTECT(allocVector(VECSXP, 2));pro++;",
                       "SET_VECTOR_ELT(mat0n, 0, state);",
                       "SET_VECTOR_ELT(mat0n, 1, state);",
                       "setAttrib(mat0, R_DimNamesSymbol, mat0n);",
                       "SET_VECTOR_ELT(matLst, 0, mat0);",
                       "SET_VECTOR_ELT(lst,  18, matLst);"
                       )
    }
    ## Generate extra RxODE code for generating right functions.
    .code <- paste(.code, collapse="\n");
    ## Generate C code for .ret0 and .ret1
    .ret <- list(.ret0, .ret1, .code,
                 paste(.codeSave, collapse="\n"));
    return(.ret);
}
