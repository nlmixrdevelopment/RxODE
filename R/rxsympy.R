regIdentifier1 <- rex::rex(one_of("_", "a":"z", "A":"Z"), any_of("_", "a":"z", "A":"Z", "0":"9", "."))
regIdentifier2 <- rex::rex(at_least(".",1), one_of("_", "a":"z", "A":"Z"), any_of("_", "a":"z", "A":"Z", "0":"9", "."));
regIdentifier <- rex::rex(or(regIdentifier1, regIdentifier2));
regSens <- rex::rex("__sens_", capture(regIdentifier), "_",  capture(regIdentifier), "__");
regToSens1 <- rex::rex( capture(regIdentifier), or("_", ".", ":"),  capture(regIdentifier));
regToSens2 <- rex::rex( "d/dt(d", capture(regIdentifier), "/d",  capture(regIdentifier), ")");
regJac <- rex::rex( "df(", capture(regIdentifier), ")/dy(",  capture(regIdentifier), ")");
regFloat1 <- rex::rex(or(group(some_of("0":"9"), ".", any_of("0":"9")),
                         group(any_of("0":"9"), ".", some_of("0":"9"))),
                      at_most(group(or("E", "e"), at_most(or("+", "-"), 1), some_of("0":"9")), 1));
regFloat2 <- rex::rex(some_of("0":"9"), or("E", "e"), at_most(or("-", "+"), 1), some_of("0":"9"));
regDecimalint <- rex::rex(or("0", group("1":"9", any_of("0":"9"))))
regNum <- sprintf("-?(?:%s|%s|%s)", regDecimalint, regFloat1, regFloat2)

## Currently unsupported sympy functions
sympy.bad <- c(
    "Chi",
    "Ci",
    "DiracDelta",
    "E1",
    "Eijk",
    "Heaviside",
    "KroneckerDelta",
    "LeviCivita",
    "Li",
    "Shi",
    "Si",
    "Ynm",
    "Ynm_c",
    "Znm",
    "airyai",
    "airyaiprime",
    "airybi",
    "airybiprime",
    "assoc_laguerre",
    "assoc_legendre",
    "bspline_basis",
    "bspline_basis_set",
    "chebyshevt",
    "chebyshevt_poly",
    "chebyshevt_root",
    "chebyshevu",
    "chebyshevu_poly",
    "chebyshevu_root",
    "dirichlet_eta",
    "ei",
    "elliptic_k",
    "erf2",
    "erf2inv",
    "erfcinv",
    "erfi",
    "erfinv",
    "expint",
    "fresnelc",
    "fresnels",
    "gegenbauer",
    "gegenbauer",
    "gegenbauer",
    "gegenbauer_poly",
    "hankel1",
    "hankel2",
    "hermite",
    "hermite_poly",
    "hyper",
    "jacobi",
    "jacobi_normalized",
    "jacobi_poly",
    "jn",
    "jn_zeros",
    "laguerre",
    "laguerre_poly",
    "legendre",
    "legendre_poly",
    "lerchphi",
    "li",
    "loggamma",
    "lowergamma",
    "mathieuc",
    "mathieucprime",
    "mathieus",
    "mathieusprime",
    "meijerg",
    "polylog",
    "uppergamma",
    "yn",
    "zeta")

regBadSymPy <- rex::rex(capture(or(sympy.bad)), "(")

##' Converts model specification into a SymPy model lines
##'
##' @param model RxODE family of objects
##' @return Lines for sympy
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxToSymPy <- function(model){
    if (class(model) == "character"){
        line <- model
        is.char <- TRUE;
    } else {
        line <- strsplit(rxModelVars(model)$model["normModel"],"\n")[[1]]
        is.char <- FALSE;
    }
    ret <- gsub(rex::rex(";", any_spaces, end), "", line);
    ret <- gsub(rex::rex("^"), "**", ret);
    id <- regIdentifier;
    ret <- gsub(rex::rex("d/dt(", any_spaces, capture(id), any_spaces, ")"), "__d_dt_\\1__", ret);
    ret <- gsub(rex::rex("df(",any_spaces, capture(id), any_spaces, ")/dy(", any_spaces, capture(id), any_spaces, ")"), "__d_df_\\1_dy_\\2__", ret);
    if (!is.char){
        ret <- ret[regexpr(rex::rex(start,any_spaces,or(names(rxInits(model))),any_spaces,"="), ret) == -1];
        ret <- ret[regexpr(rex::rex(start,any_spaces,or(names(rxInits(model))),"(0)", any_spaces,"="), ret) == -1];
    }
    return(ret);
}
##' Convert SymPy syntax to RxODE syntax
##'
##' @param text SymPy text to convert to RxODE syntax
##' @return RxODE syntax
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxFromSymPy <- function(text){
    ret <- paste0(gsub(rex::rex("**"), "^", text), ";");
    if (regexpr(regBadSymPy, ret) != -1){
        stop(sprintf("This caluclation requies the sympy function '%s' which is not currently implemented in RxODE",
                     gsub(rex::rex(anything, regBadSymPy, anything), "\\1",ret)))
    }
    return(ret);
}
rxSymPy.vars <- c();
##' Setup sympy variables
##'
##' This creates sympy variables for later evaulation in the CAS sympy.
##'
##' @param model RxODE family of objects
##' @return NULL
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPyVars <- function(model){
    if (!exists(".Jython", .GlobalEnv))
        rSymPy::sympyStart();
    if (!class(model) == "character"){
        vars <- c(rxParams(model), rxState(model));
    } else {
        vars <- model
    }
    known <- c(rxSymPy.vars, vars);
    assignInMyNamespace("rxSymPy.vars", known);
    .Jython$exec(sprintf("%s = symbols('%s')", paste(vars, collapse=", "), paste(vars, collapse=" ")));
    return(invisible());
}

##' Setup a sympy environment that sets up the specified RxODE model.
##'
##' @param model RxODE family of objects
##' @return NULL
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPySetup <- function(model){
    rxSymPyVars(model);
    for (line in rxToSymPy(model)){
        if (regexpr(rex::rex(any_spaces, capture(regIdentifier), any_spaces, "="), line) != -1){
            assignInMyNamespace("rxSymPy.vars", c(rxSymPy.vars, gsub(rex::rex(any_spaces, capture(regIdentifier), any_spaces, "=",anything), "\\1", line)));
        }
        .Jython$exec(line);
    }
    return(invisible());
}

##' Calculate df/dy derivatives
##'
##' @param model RxODE models
##' @param df is a string for the state in the df(.)/dy(.).  If
##'     missing and dy is missing, all the df(.)/dy(.) components are
##'     calulated according to the \code{vars} parameter below.
##' @param dy is a string for the state or varaible in the df(.)/dy(.).
##' @param vars is a boolean indicating if parameters will be included
##'     for the dy component in the df(.)/dy(.), instead of just state
##'     variables (required for sensitivity equations).
##' @return RxODE syntax lines for the df(.)/dy(.)
##' @author Matthew L. Fidler
##' @export
##' @keywords internal
rxSymPyDfDy <- function(model, df, dy, vars=FALSE){
    if (missing(df) && missing(dy)){
        if (class(vars) == "logical"){
            if (vars){
                jac <- expand.grid(s1=rxState(model), s2=c(rxState(model), rxParams(model)));
            } else  {
                jac <- expand.grid(s1=rxState(model), s2=rxState(model));
            }
        }
        jac <- with(jac, data.frame(jac,
                          rx=sprintf("df(%s)/dy(%s)", s1, s2),
                          sym=sprintf("__d_df_%s_dy_%s__", s1, s2)))
        rxSymPySetup(model);
        extraLines <- c();
        cat("Calculate Jacobian.");
        for (dfdy in jac$rx){
            if (!any(dfdy == rxDfdy(model))){
                extraLines[length(extraLines) + 1] <- with(jac[jac$rx == dfdy, ], rxSymPyDfDy(NULL, s1, s2));
                cat(".")
            }
        }
        cat("done!\n");
        return(extraLines);
    } else {
        if (!is.null(model)){
            rxSymPySetup(model);
        }
        line <- rSymPy::sympy(sprintf("diff(%s,%s)", rxToSymPy(sprintf("d/dt(%s)", df)), dy));
        var <- rxToSymPy(sprintf("df(%s)/dy(%s)", df, dy));
        assignInMyNamespace("rxSymPy.vars", c(rxSymPy.vars, var));
        .Jython$exec(sprintf("%s=%s", var, line));
        return(sprintf("df(%s)/dy(%s) = %s", df, dy, rxFromSymPy(line)));
    }
}
##' Calculate the full jacobain for a model
##'
##' This expand the model to caluclate the jacobian.  This requires
##' rSymPy.
##'
##' @param model RxODE family of objects
##' @return RxODE syntax for model with jacobian specified.
##' @author Matthew L. Fidler
rxSymPyJacobian <- function(model){
    extraLines <- rxSymPyDfDy(model, vars=FALSE);
    extraLines <- extraLines[regexpr(rex::rex("=", any_spaces, "0", any_spaces, ";"), extraLines) == -1];
    model <- sprintf("%s\n%s", rxModelVars(model)$model["normModel"], paste(extraLines, collapse="\n"));
    rxSymPyClean()
    return(model);
}
##' Calculate the sensitivity equations for a model
##'
##' This expands the model to calculate sensitivities.  This requires
##' rSymPy.
##'
##' @param model RxODE family of objects
##' @return Model syntax that includes the sensitivity parameters.
##' @author Matthew L. Fidler
##' @export
rxSymPySensitivity <- function(model){
    calcSens <- rxParams(model);
    state <- rxState(model)
    extraLines <- rxSymPyDfDy(model, vars=TRUE);
    all.sens <- c();
    cat("Calculate Sensitivities.");
    for (s1 in state){
        for (sns in calcSens){
            tmp <- c()
            vars <- c();
            for (s2 in state){
                vars <- c(vars, sprintf("__sens_%s_%s__", s2, sns));
                extra <- sprintf("df(%s)/dy(%s)*__sens_%s_%s__", s1, s2, s2, sns)
                tmp <- c(tmp, extra);
            }
            tmp <- c(tmp, sprintf("df(%s)/dy(%s)", s1, sns));
            rxSymPyVars(vars);
            all.sens <- c(all.sens, vars);
            line <- rSymPy::sympy(rxToSymPy(paste(tmp, collapse=" + ")));
            var.rx <- sprintf("d/dt(__sens_%s_%s__)", s1, sns)
            var <- rxToSymPy(var.rx)
            .Jython$exec(sprintf("%s=%s", var, line));
            assignInMyNamespace("rxSymPy.vars", c(rxSymPy.vars, var));
            extraLines[length(extraLines) + 1] <- sprintf("%s=%s", var.rx, rxFromSymPy(line));
            cat(".")
        }
    }
    cat("done!\n");
    cat("Expanding Jacobian for sensitivities.")
    jac2 <- expand.grid(s1=unique(all.sens), s2=unique(c(all.sens, rxState(model))))
    for (i in 1:length(jac2$s1)){
        extraLines[length(extraLines) + 1] <- rxSymPyDfDy(NULL, jac2$s1[i], jac2$s2[i]);
        cat(".");
        if (i %% 5 == 0){
            cat(i);
        }
        if (i %% 50 == 0){
            cat("\n");
        }
    }
    cat("done!\n");
    extraLines <- extraLines[regexpr(rex::rex("=", any_spaces, "0", any_spaces, ";"), extraLines) == -1];
    model <- sprintf("%s\n%s", rxModelVars(model)$model["normModel"], paste(extraLines, collapse="\n"));
    rxSymPyClean()
    return(model);
}

##' Remove variables created by RxODE from the sympy environment.
##'
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPyClean <- function(){
    for (v in unique(rxSymPy.vars)){
        try({.Jython$exec(sprintf("del %s", v))}, silent = TRUE);
    }
    assignInMyNamespace("rxSymPy.vars", c(rxSymPy.vars, var));
}

## Supported Sympy special functions
## besseli -> besseli(nu,z) -> bessel_i(z,nu,1)
## besselj -> besselj(nu,z) -> bessel_j(z,nu)
## besselk -> besselk(nu,z) -> bessel_k(z,nu,1)
## bessely -> bessely(nu,z) -> bessel_y(z,nu)
## beta -> beta
## digamma -> digamma
## erfc -> erfc
## gamma -> gammafn
## polygamma -> polygamma(n, z) returns log(gamma(z)).diff(n + 1) = pigamma(z, n)
## trigamma -> trigamma
