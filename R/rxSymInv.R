rxSymInvC <- function(mat1, diag.xform=c("sqrt", "log", "identity"), chol=FALSE){
    if (!all(as.vector(mat1) == 1 || as.vector(mat1) == 1)){
        stop("This has to be a matrix of all 1s or 0s.");
    }
    if (any(diag(mat1) == 0)){
        stop("Diagonal elements must be non-zero.");
    }
    diag.xform <- match.arg(diag.xform)
    cache.file <- file.path(ifelse(RxODE.cache.directory == ".", getwd(), RxODE.cache.directory),
                            sprintf("rx_%s.inv",
                                    digest::digest(deparse(list(mat1, diag.xform, chol)))));
    cache.file2 <- file.path(system.file("inv", package="RxODE"), cache.file);
    if (file.exists(cache.file)){
        load(cache.file);
        return(ret)
    } else if (file.exists(cache.file2)){
        load(cache.file2);
        return(ret)
    } else {
        rxCat("Diagonal form: ", diag.xform, "\n");
        rxPrint(mat1);
        num <- as.vector(mat1[upper.tri(mat1,TRUE)]);
        i <- 0;
        num <- sapply(num, function(x){
            if (x == 1){
                i <<- i + 1;
                return(i)
            } else {
                return(0)
            }
        })
        mat1[upper.tri(mat1, TRUE)] <- num -1;
        mat1[lower.tri(mat1)] <- t(mat1)[lower.tri(mat1)];
        d <- dim(mat1)[1];
        mat1 <- paste0("t", mat1);
        mat1[mat1 == "t-1"] <- "0";
        mat1 <- matrix(mat1, d);

        if (diag.xform == "sqrt"){
            diag(mat1) <- sprintf("%s**2", diag(mat1))
        } else if (diag.xform == "log"){
            diag(mat1) <- sprintf("exp(%s)", diag(mat1))
        } else {
            diag(mat1) <- sprintf("(%s)", diag(mat1));
        }
        omat <- fmat <- mat1;
        sympy.mat <- sprintf("Matrix([%s])", paste(apply(mat1, 1, function(x){
                                                 return(sprintf("[%s]", paste(x, collapse=", ")))
                                             }), collapse=", "))
        vars <- paste0("t", seq(0, i - 1))
        syms <- paste(vars, collapse=", ");
        if (length(vars) == 1){
            rxSymPyExec(sprintf("%s = Symbol('%s')",syms, syms));
        } else {
            rxSymPyExec(sprintf("%s = symbols('%s')",syms, syms));
        }
        if (chol){
            rxCat("Calculate symbolic inverse:  t(chol.mat) %*% chol.mat ...\n");
            sympy.inv <- rxSymPy(sprintf("((Matrix([%s])).transpose()).multiply(Matrix([%s]))", sympy.mat, sympy.mat));
        } else {
            rxCat("Calculate symbolic inverse...");
            sympy.inv <- rxSymPy(sprintf("(%s).inv()", sympy.mat));
            if (sympy.inv == "None"){
                stop("Inverse not calculated.");
            }
        }
        sympy.inv <- gsub("[\n\t ]*", "", sympy.inv)
        sympy.inv.txt <- sprintf("Matrix([%s])", gsub("[\n]", ", ", sympy.inv));
        sympy.inv.det.tmp1 <- sprintf("(%s).det()", sympy.inv.txt);
        mat.reg <- rex::rex(start, or(group(any_spaces, "["),
                                      group(any_spaces, "Matrix([[")), any_spaces,
                            capture(anything), any_spaces,
                            or(group("]", any_spaces, end),
                               group("]])", any_spaces, end)));
        mat.sep.reg <- rex::rex(or(group("]", any_spaces, ",", any_spaces, "["),
                                   group(any_spaces, ",", any_spaces)))
        sympy.inv <- gsub(mat.reg, "\\1", strsplit(sympy.inv, "\n")[[1]]);
        sympy.inv <- matrix(unlist(strsplit(sympy.inv, mat.sep.reg)), d, byrow=TRUE);
        rxCat("done\n");
        if (chol){
            rxCat("Calculate Omega in terms of chol(Omega^-1) parameterization...\n");
            sympy.inv.inv <- rxSymPy(sprintf("(%s).inv()", sympy.inv.txt));
            if (sympy.inv == "None"){
                stop("Inverse not calculated.");
            }
            sympy.inv.inv.txt <- sprintf("Matrix([%s])", gsub("[\n]", ", ", sympy.inv.inv));
            sympy.inv.inv <- gsub(mat.reg, "\\1", strsplit(sympy.inv.inv, "\n")[[1]]);
            sympy.inv.inv <- matrix(unlist(strsplit(sympy.inv.inv, mat.sep.reg)), d, byrow=TRUE);
            ch <- omat;
            omat <- fmat <- sympy.inv.inv;
            rxCat("done\n");
            rxCat("Calculate log(det(OMGAinv))...")
            sympy.inv.det <- paste(sprintf("log(%s)", diag(ch)), collapse=" + ");
            sympy.inv.det <- sympyC(sympy.inv.det);
            rxCat("done\n");
        }
        if (d <= 3 && !chol){
            rxCat("Calculate symbolic determinant of inverse...");
            sympy.inv.det <- rxSymPy(sympy.inv.det.tmp1)
            sympy.inv.det <- sympyC(sympy.inv.det);
            rxCat("done\n");
        } else if (!chol) {
            sympy.inv.det <- "NA_REAL";
        }
        v <- vars[1]
        rxCat("Calculate d(Omega)/d(Est) and d(Omega^-1)/d(Est)...\n");
        cnt.i <- 0;
        cnt <- function(){
            rxCat(".");
            if (cnt.i %% 5 == 0)
                rxCat(cnt.i);
            if (cnt.i %% 50 == 0)
                rxCat("\n");
            cnt.i <<- cnt.i + 1;
        }
        i <- 0;
        omega <- sprintf("  if (theta_n == 0){\n%s\n  }",
                         paste(sapply(as.vector(omat), function(x){
                             cnt();
                             ret <- sprintf("    REAL(ret)[%s] = %s;", i, sympyC(x));
                             i <<- i + 1;
                             return(ret);
                         }), collapse="\n"))
        i <- 0;
        omega1 <- sprintf("  if (theta_n == 0){\n%s\n  }",
                          paste(sapply(as.vector(sympy.inv), function(x){
                              cnt();
                              ret <- sprintf("    REAL(ret)[%s] = %s;", i, sympyC(x));
                              i <<- i + 1;
                              return(ret);
                          }), collapse="\n"))
        j <- 0;
        omegap <- paste(unlist(lapply(vars, function(v){
            i <<- 0;
            j <<- j + 1;
            return(sprintf("  if (theta_n == %s){\n%s  \n}", j,
                           paste(sapply(as.vector(omat),
                                        function(x){
                               if (x == "0"){
                                   ret <- "0";
                               } else {
                                   cnt()
                                   ret <- rxSymPy(sprintf("diff(%s,%s)", x, v));
                               }
                               ret <- sprintf("    REAL(ret)[%s] = %s;", i, sympyC(ret));
                               i <<- i + 1;
                               return(ret);
                           }),
                           collapse="\n")))})), collapse="\n");
        j <- 0;
        omega1p <- paste(unlist(lapply(vars, function(v){
            i <<- 0;
            j <<- j + 1;
            return(sprintf("  if (theta_n == %s){\n%s\n  }", j,
                           paste(sapply(as.vector(sympy.inv),
                                        function(x){
                               if (x == "0"){
                                   ret <- "0";
                               } else {
                                   cnt()
                                   ret <- rxSymPy(sprintf("diff(%s,%s)", x, v));
                               }
                               ret <- sprintf("    REAL(ret)[%s] = %s;", i, sympyC(ret));
                               i <<- i + 1;
                               return(ret);
                           }),
                           collapse="\n")))})), collapse="\n");
        rxCat("done.\n")
        ret <- paste(c("int omega_i = INTEGER(oi)[0];",
                       "int theta_n = INTEGER(tn)[0];",
                       "if (omega_i == 0){",
                       "SEXP ret = PROTECT(allocVector(REALSXP,1));",
                       "if (theta_n == 0){",
                       sprintf("REAL(ret)[0] = %d;", d),
                       "} else if (theta_n == -1){",
                       sprintf("REAL(ret)[0] = %s;", sympy.inv.det),
                       ifelse(chol, "// 0.5*log(det(Omega^-1)) Directly estimated.",
                              paste(c("if (!ISNA(REAL(ret)[0])){",
                                "if (REAL(ret)[0] > 0){",
                                "REAL(ret)[0] = 0.5*log(REAL(ret)[0]);",
                                "} else {",
                                "error(\"Omega^-1 not positive definite\");",
                                "}",
                                "}"), collapse="\n")),
                       "} else {",
                       sprintf("REAL(ret)[0] = %d;", length(vars)),
                       "}",
                       "UNPROTECT(1);",
                       "return ret;",
                       "} else {",
                       sprintf("SEXP ret = PROTECT(allocMatrix(REALSXP, %s, %s));", d, d),
                       "if (omega_i == 1){",
                       omega,
                       omegap,
                       "} else if (omega_i == -1){",
                       omega1,
                       omega1p,
                       "}",
                       "UNPROTECT(1);",
                       "return ret;}"), collapse="\n");
        for (v in vars){
            rxSymPyClear(v);
        }
        fmat <- matrix(sapply(as.vector(fmat), function(x){force(x);return(rxFromSymPy(x))}), d);
        ret <- list(ret, fmat);
        save(ret, file=cache.file);
    }
    return(ret);
}

rxSymInvCreateC <- function(c){
    inline::cfunction(signature(theta="numeric", oi="integer", tn="integer"), c)
}

rxSymInvCreateC.slow <- NULL

##' Creates an object for caluclating Omega/Omega^-1 and dervatives
##'
##' @param mat Initial Omega matrix
##' @param diag.xform transformation to diagonal elements of OMEGA. or chol(Omega^-1)
##' @param chol Boolean to state if the parameter values specify the OMEGA matrix or chol(Omega^-1)
##' @return A rxSymInv object
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymInvCreate <- function(mat,
                           diag.xform=c("sqrt", "log", "identity"),
                           chol=FALSE){
    diag.xform <- match.arg(diag.xform);
    mat2 <- mat;
    if (chol){
        mat2 <- rxInv(mat2);
        mat2 <- try({chol(mat2)});
        if (inherits(mat2, "try-error")){
            rxCat("Warning: Initial Omega matrix inverse is non-positive definite, correcting with nearPD.")
            mat2 <- as.matrix(Matrix::nearPD(mat)$mat);
            mat2 <- chol(mat2);
        }
    }
    mat3 = mat2;
    if (diag.xform == "sqrt"){
        diag(mat3) <- sqrt(diag(mat3));
    } else if (diag.xform == "log"){
        diag(mat3) <- log(diag(mat3));
    }
    elts <- as.vector(mat3)[which(as.vector(lower.tri(mat3,TRUE))*1==1)];
    th.unscaled <- c();
    for (i in 1:length(elts)){
        if (elts[i] != 0){
            th.unscaled[length(th.unscaled) + 1] <- elts[i];
        }
    }
    mat1 <-(mat != 0)*1;
    if (length(mat1) == 1){
        mat1 <- matrix(mat1, 1);
    }
    dmat <- dim(mat1)[1] -1;
    block <- list();
    last <- 1;
    if (dmat != 0){
        for (i in 1:dmat){
            if (all(mat1[rxBlockZeros(mat1,i)] == 0)){
                s <- seq(last, i);
                cur <-matrix(as.double(mat[s, s]), length(s));
                last <- i + 1;
                block[[length(block) + 1]] <- cur;
            }
        }
    }
    if (length(block) != 0){
        s <- seq(last, dmat + 1);
        cur <- matrix(as.double(mat[s, s]), length(s));
        block[[length(block) + 1]] <- cur;
    }
    if (length(block) == 0){
        ret <-rxSymInvC(mat1=mat1,
                        diag.xform=diag.xform, chol=chol);
        th <- th.unscaled;
        ret <- list(fmat=ret[[2]],
                    chol=chol,
                    th=th.unscaled,
                    fn=rxSymInvCreateC(ret[[1]]));
        class(ret) <- "rxSymInv";
        return(ret);
    } else {
        mat <- Matrix::.bdiag(block);
        matI <- lapply(block, rxSymInvCreate, diag.xform=diag.xform, chol=chol);
        ret <- list(mat=mat,
                    matI=matI,
                    th=th.unscaled,
                    chol=chol);
        class(ret) <- "rxSymInvBlock"
        return(ret);
    }
}

##' @export
print.rxSymInvBlock <- function(x, ...){
    d <- dim(x$mat)[1]
    rxCat(sprintf("Object to create Omega and Omega^-1 & derivitaves for a %sx%s matrix:\n", d, d))
    rxPrint(x$mat);
    rxCat("Use `rxSymInv' for the matrix.\n");
}

##' @export
print.rxSymInv <- function(x, ...){
    d <- dim(x$fmat)[1]
    rxCat(sprintf("Object to create Omega and Omega^-1 & derivitaves for a %sx%s matrix:\n", d, d))
    rxPrint(x$fmat);
    rxCat("Use `rxSymInv' for the matrix.\n");
}

##' Get the theta numbers for the diagonal elements.
##'
##' These start at element 1 like in R instead of element 0 like in C
##'
##' @param x object created with rxSymInvCreate
##' @return Theta numbers for diagonal elements
##' @author Matthew L. Fidler
##' @export
rxSymDiag <- function(x){
    if (class(x) == "rxSymInv"){
        return(as.numeric(gsub(rex::rex(anything, "t", capture(any_numbers),anything), "\\1", diag(x$fmat))) + 1)
    } else if (class(x) == "rxSymInvBlock"){
        mx <- 0;
        ret <- c()
        for (i in seq_along(x$matI)){
            ret <- c(ret, rxSymDiag(x$matI[[i]]) + mx);
            mx <- ret[length(ret)];
        }
        return(ret);
    }
}

##' Get Omega and Omega^-1 and derivatives
##'
##' @param invobj rxSymInv object
##' @param theta Thetas to be used for calculation
##' @param pow Power of Omega, either 1 or -1 for inverse.  If 0,
##'     return an environment with all possible deriatives and Omegas
##'     caluclated
##' @param dTheta the theta number to take the derivative.  That is
##'     d(Omega)/dTheta.  When dTheta = 0, just return the Omega or
##'     Omega^-1 matrix.  Also can be 0.5.  When 0.5 and pow is -1,
##'     return -1/2*log(det(omegaInv)).
##' @return Matrix based on parameters or environment with all the
##'     matrixes calculated in variables omega, omegaInv, dOmega,
##'     dOmegaInv.
##' @author Matthew L. Fidler
##' @export
rxSymInv <- function(invobj, theta, pow=0, dTheta=0){
    if (class(invobj) == "rxSymInv"){
        if (length(theta) == 1){
            if (theta == "n"){
                return(invobj$fn(0, 0L, 1L));
            }
        }
        if (!any(pow == c(1, 0, -1))){
            stop("The power can only be 1, 0, or -1");
        }
        if (pow == -1 && dTheta == 0.5){
            return(invobj$fn(theta, 0L, -1L))
        }
        ntheta <- invobj$fn(0, 0L, 1L);
        if (dTheta < 0 || round(dTheta) != dTheta || dTheta > ntheta) {
            stop(sprintf("dTheta can be any integer from 0 to %s", ntheta));
        }
        if (length(theta) != ntheta){
            stop(sprintf("This requires %d theta estimates.", ntheta));
        }
        if (pow == 0){
            ret <- new.env(parent=emptyenv());
            theta <- as.double(theta);
            ret$omega <- invobj$fn(theta, 1L, 0L);
            ret$omegaInv <- invobj$fn(theta, -1L, 0L);
            ret$dOmega <- list();
            ret$dOmegaInv <- list();
            for (i in seq(1, ntheta)){
                ret$dOmega[[i]] <- invobj$fn(theta, 1L, as.integer(i));
                ret$dOmegaInv[[i]] <- invobj$fn(theta, -1L, as.integer(i));
            }
            tryCatch({ret$log.det.OMGAinv.5 <- invobj$fn(theta,0L, -1L)},
                     error=function(e){
                rxCat("Warning: Omega^-1 not positive definite (correcting with nearPD).\n");
                old <- ret$omegaInv
                ret$omegaInv <- as.matrix(Matrix::nearPD(ret$omegaInv)$mat);
                RxODE_finalize_log_det_OMGAinv_5(ret);
                ret$omegaInv <- old;
            })
            if (is.na(ret$log.det.OMGAinv.5)){
                old <- ret$omegaInv
                ret$omegaInv <- as.matrix(Matrix::nearPD(ret$omegaInv)$mat);
                RxODE_finalize_log_det_OMGAinv_5(ret);
                ret$omegaInv <- old;
            }
            RxODE_finalize_focei_omega(ret);
            return(ret)
        } else {
            return(invobj$fn(as.double(theta), as.integer(pow), as.integer(dTheta)));
        }
    } else if (class(invobj) == "rxSymInvBlock"){
        if (!any(pow == c(1, 0, -1))){
            stop("The power can only be 1, 0, or -1");
        }
        nthetas <- lapply(invobj$matI, function(x){
            x$fn(0, 0L, 1L)
        });
        tot.nthetas <- sum(unlist(nthetas))
        if (length(theta) == 1){
            if (theta == "n"){
                return(tot.nthetas);
            }
        }
        if (dTheta < 0 || round(dTheta) != dTheta || dTheta > tot.nthetas) {
            stop(sprintf("dTheta can be any integer from 0 to %s.", ntheta));
        }
        if (length(theta) != tot.nthetas){
            stop(sprintf("This requires %d theta estimates.", tot.nthetas));
        }
        all.theta <- theta;
        theta.decomp <- lapply(invobj$matI, function(x) {
            s <- seq(1, x$fn(0, 0L, 1L));
            ret <- all.theta[s]
            all.theta <<- all.theta[-s];
            return(ret);
        });
        nn <- 1;
        theta.i <- lapply(invobj$matI, function(x) {
            s <- seq(nn, nn + x$fn(0, 0L, 1L) - 1);
            nn <<- max(s) + 1;
            return(s);
        });
        blockfn <- function(pow, dTheta){
            if (pow == 0L && dTheta == -1L){
                ##
                ret <- lapply(1:length(invobj$matI), function(x){
                    tmp <- invobj$matI[[x]];
                    th <- theta.decomp[[x]];
                    return(tryCatch(tmp$fn(as.double(th),0L, -1L), error=function(e){return(NA)}))
                })
                ret <- unlist(ret)
                return(sum(ret));
            } else {
                ret <- lapply(1:length(invobj$matI), function(x){
                    tmp <- invobj$matI[[x]];
                    th <- theta.decomp[[x]];
                    if (dTheta > 0){
                        thi <- theta.i[[x]];
                        if (dTheta > max(thi) || dTheta < min(thi)){
                            d <- dim(tmp$fmat)[1];
                            return(matrix(rep(0, d * d), d));
                        } else {
                            return(rxSymInv(tmp, th, pow, dTheta - min(thi) + 1))
                        }
                    }
                    return(rxSymInv(tmp, th, pow, dTheta));
                });
                return(as.matrix(Matrix::bdiag(ret)));
            }

        }
        ret <- new.env(parent=emptyenv());
        if (pow == 0){
            ret <- new.env(parent=emptyenv());
            theta <- as.double(theta);
            ret$omega <- blockfn(1L, 0L);
            ret$omegaInv <- blockfn(-1L, 0L);
            ret$dOmega <- list();
            ret$dOmegaInv <- list();
            for (i in seq(1, tot.nthetas)){
                ret$dOmega[[i]] <- blockfn(1L, as.integer(i));
                ret$dOmegaInv[[i]] <- blockfn(-1L, as.integer(i));
            }
            ret$log.det.OMGAinv.5 <- blockfn(0L, -1L);
            if (is.na(ret$log.det.OMGAinv.5)){
                rxCat("Warning: Omega^-1 may not positive definite (correcting with nearPD)\n");
                old <- ret$omegaInv
                ret$omegaInv <- as.matrix(Matrix::nearPD(ret$omegaInv)$mat);
                RxODE_finalize_log_det_OMGAinv_5(ret);
                ret$omegaInv <- old;
            }
            RxODE_finalize_focei_omega(ret);
            return(ret);
        } else {
            return(blockfn(as.integer(pow), as.integer(dTheta)));
        }
    } else {
        stop("This needs to be applied on an object created with 'rxSymInvCreate'.");
    }
}

##' Creates a logical matrix for block matrixes.
##'
##' @param mat Matrix
##' @param i Row/column where block matrix should be setup.
##'
##' @return A logical matrix returning where the elements should be
##'     zero.
##'
##' @keywords internal
##' @export
rxBlockZeros <- function(mat, i){
    return(!((row(mat) > i & col(mat) > i) | (row(mat) <= i & col(mat) <= i)))
}
rxIsBlock <- function(mat, i){
    if (missing(i)){
        for (j in 1:(dim(mat) - 1)){
            if (rxIsBlock(mat, j)){
                return(TRUE)
            } else {
                return(FALSE);
            }
        }
    } else {
        return(all(mat[rxBlockZeros(mat, i)] == 0));
    }
}
