rxSymInvC.slow <- NULL;  ## Memoize
rxSymInvC <- function(mat1, diag.xform=c("sqrt", "log", "identity")){
    if (!all(as.vector(mat1) == 1 || as.vector(mat1) == 1)){
        stop("This has to be a matrix of all 1s or 0s.");
    }
    if (any(diag(mat1) == 0)){
        stop("Diagonal elements must be non-zero.");
    }
    diag.xform <- match.arg(diag.xform)
    cache.file <- sprintf("rx_%s.inv",
                          digest::digest(deparse(list(mat1, diag.xform))));
    cache.file2 <- file.path(system.file("inv", package="RxODE"), cache.file);
    if (file.exists(cache.file)){
        load(cache.file);
        return(ret)
    } else if (file.exists(cache.file2)){
        load(cache.file2);
        return(ret)
    } else {
        cat("Diagional form: ", diag.xform, "\n");
        print(mat1);
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
        }
        omat <- fmat <- mat1;
        sympy.mat <- sprintf("Matrix([%s])", paste(apply(mat1, 1, function(x){
                                                 return(sprintf("[%s]", paste(x, collapse=", ")))
                                             }), collapse=", "))
        vars <- paste0("t", seq(0, i - 1))
        syms <- paste(vars, collapse=", ");
        if (length(vars) == 1){
            rxSympyExec(sprintf("%s = Symbol('%s')",syms, syms));
        } else {
            rxSympyExec(sprintf("%s = symbols('%s')",syms, syms));
        }
        rxCat("Calculate symbolic inverse...");
        sympy.inv <- rxSymPy(sprintf("(%s).inv()", sympy.mat));
        sympy.inv.txt <- sprintf("Matrix([%s])", gsub("[\n]", ", ", sympy.inv));
        sympy.inv.det.tmp <- sprintf("(%s).det()", sympy.inv.txt);
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
        rxCat("Calculate symbolic determinant of inverse...");
        sympy.inv.det <- "NA_REAL";
        R.utils::withTimeout({
                     sympy.inv.det <- rxSymPy(sympy.inv.det.tmp);
                     sympy.inv.det <<- sympyC(sympy.inv.det);
                 }, timeout=60 * 70, onTimeout="warning");
        rxCat("done\n");
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
                       "if (!ISNA(REAL(ret)[0])){",
                       "if (REAL(ret)[0] > 0){",
                       "REAL(ret)[0] = 0.5*log(REAL(ret)[0]);",
                       "} else {",
                       "error(\"Omega^-1 not positive definite\");",
                       "}",
                       "}",
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

##' Creates an object for caluclating Omega/Omega^-1 and dervatives
##'
##' @param mat Initial Omega matrix
##' @param diag.xform transformation to diagonal elements of OMEGA.
##' @return A rxSymInv object
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymInvCreate <- function(mat,
                           diag.xform=c("sqrt", "log", "identity")){
    diag.xform <- match.arg(diag.xform);
    mat2 <- mat;
    elts <- as.vector(mat)[which(as.vector(lower.tri(mat,TRUE))*1==1)];
    th.unscaled <- c();
    for (i in 1:length(elts)){
        if (elts[i] != 0){
            th.unscaled[length(th.unscaled) + 1] <- elts[i];
        }
    }
    ret <-rxSymInvC(mat1=(mat>0)*1,
                    diag.xform=diag.xform);
    th <- th.unscaled;
    ret <- list(fmat=ret[[2]],
                th.unscaled=th.unscaled,
                th=th,
                fn=inline::cfunction(signature(theta="numeric", oi="integer", tn="integer"), ret[[1]]));
    class(ret) <- "rxSymInv";
    return(ret);
}

##' @export
print.rxSymInv <- function(x, ...){
    d <- dim(x$fmat)[1]
    cat(sprintf("Object to create Omega and Omega^-1 & derivitaves for a %sx%s matrix:\n", d, d))
    print(x$fmat);
    cat("Use `rxSymInv' for the matrix.\n");
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
    if (class(invobj) != "rxSymInv"){
        stop("This needs to be applied on an object created with 'rxSymInvCreate'.");
    } else {
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
                cat("Warning: Omega^-1 not positive definite (correcting with nearPD)\n");
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
    }
}


rxBlockZeros <- function(mat, i){
    return(!((row(mat) > i & col(mat) > i) | (row(mat) <= i & col(mat) <= i)))
}

rxIsBlock <- function(mat, i){
    all(mat[rxBlockZeros(mat, i)] == 0);
}
