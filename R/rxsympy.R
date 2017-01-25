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
    return(ret);
}
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
    if (!class(model) == "character"){
        vars <- c(rxParams(model), rxState(model));
    } else {
        vars <- model
    }
    rSymPy::sympy(sprintf("%s = symbols('%s')", paste(vars, collapse=", "), paste(vars, collapse=" ")));
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
        rSymPy::sympy(line);
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
        ## cat("Calculate Jacobian.");
        for (dfdy in jac$rx){
            if (!any(dfdy == rxDfdy(model))){
                extraLines[length(extraLines) + 1] <- with(jac[jac$rx == dfdy, ], rxSymPyDfDy(NULL, s1, s2));
                ## cat(".")
            }
        }
        ## cat("done!\n");
        return(extraLines);
    } else {
        if (!is.null(model)){
            rxSymPySetup(model);
        }
        line <- rSymPy::sympy(sprintf("diff(%s,%s)", rxToSymPy(sprintf("d/dt(%s)", df)), dy));
        rSymPy::sympy(sprintf("%s=%s", rxToSymPy(sprintf("df(%s)/dy(%s)", df, dy)), line))
        return(sprintf("df(%s)/dy(%s) = %s",
                       df, dy,
                       rxFromSymPy(rSymPy::sympy(sprintf("diff(%s,%s)", rxToSymPy(sprintf("d/dt(%s)", df)), dy)))));
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
    ## cat("Calculate Sensitivities.");
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
            line <- sprintf("d/dt(__sens_%s_%s__)=%s", s1, sns, rxFromSymPy(line));
            rSymPy::sympy(rxToSymPy(line));
            extraLines[length(extraLines) + 1] <- line
            ## cat(".")
        }
    }
    ## cat("done!\n");
    ## cat("Expanding Jacobian for sensitivities.")
    jac2 <- expand.grid(s1=unique(all.sens), s2=unique(c(all.sens, rxState(model))))
    for (i in 1:length(jac2$s1)){
        extraLines[length(extraLines) + 1] <- rxSymPyDfDy(NULL, jac2$s1[i], jac2$s2[i]);
        ## cat(".");
        ## if (i %% 5 == 0){
        ##     cat(i);
        ## }
        ## if (i %% 50 == 0){
        ##     cat("\n");
        ## }
    }
    ## cat("done!\n");
    extraLines <- extraLines[regexpr(rex::rex("=", any_spaces, "0", any_spaces, ";"), extraLines) == -1];
    model <- sprintf("%s\n%s", rxModelVars(model)$model["normModel"], paste(extraLines, collapse="\n"));
    return(model);
}
