utils::globalVariables(c(".Jython"))

## Based on normalized grammar output by parser which has
## if (expr) {
## }
## else {
## }
## for the parsing of if/else statements
regIf <- rex::rex(start, any_spaces, "if", any_spaces, "(", capture(anything), ")", any_spaces, "{", any_spaces, end);
regElse <- rex::rex(start, any_spaces, "else", any_spaces, "{", any_spaces, end);
regEnd <- rex::rex(start, any_spaces, "}", any_spaces, end);
regIfOrElse <- rex::rex(or(regIf, regElse))

##' Remove INIs
##' @param x RxODE list of lines to remove
##' @return RxODE lines with inis removed.
##' @author Matthew L. Fidler
##' @keywords internal
rxRmIni <- function(x){
    x <- x[regexpr(rex::rex(start,any_spaces,or(names(rxInits(x))),any_spaces,"="), x) == -1];
    x <- x[regexpr(rex::rex(start,any_spaces,or(names(rxInits(x))),"(0)", any_spaces,"="), x) == -1];
    return(x);
}
##' Remove Jacobian
##' @param x RxODE list of lines to remove
##' @return RxODE lines with df/dy removed.
##' @author Matthew L. Fidler
##' @keywords internal
rxRmJac <- function(x){
    return(x[regexpr(rex::rex(regJac, any_spaces, "="), x) == -1])
}
##' Remove sensitivity equations
##' @param x RxODE lines to remove
##' @return Lines with d/dt(rx_sens_...._) removed.
##' @author Matthew L. Fidler
rxRmSens <- function(x){
    return(x[regexpr(rex::rex(start, any_spaces, "d/dt(", any_spaces, regSens, any_spaces, ")", any_spaces, "="), x) == -1]);
}
##' Remove print statements
##' @param x RxODE lines to remove
##' @return RxODE with print lines removed.
##' @author Matthew L. Fidler
rxRmPrint <- function(x){
    return(x[regexpr(regPrint, x) == -1]);
}

##' Expand if/else clauses into mutiple different types of lines.
##'
##'
##' @param model Model can be a character, or a RxODE model.  It needs
##'     to have normalized syntax, that is \code{if (...)\{} has to be
##'     on the same line.  The \code{else} statement must be on its
##'     own line with the closing bracket of the \code{if} statement
##'     on the previous line.  This \code{else} statment must also
##'     contain the opening bracket, like the code \code{else \{}
##' @param removeInis A boolean indicating if parameter
##'     initializations should be removed from the model.
##' @param removePrint A boolean indicating if printing statements
##'     should be removed from the model.
##' @return A named character vector. The names of the vector are the
##'     logical conditions, the values are the lines that satisfy the
##'     logical conditions.
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxExpandIfElse <- function(model, removeInis=TRUE, removePrint=TRUE){
    ## expand if/else blocks into a list with lines for conditions that are true
    x <- strsplit(rxNorm(model, FALSE), "\n")[[1]];
    if (removeInis){
        x <- rxRmIni(x);
    }
    if (removePrint){
        x <- rxRmPrint(x);
    }
    model <- x;
    w1 <- which(regexpr(regIfOrElse, model) != -1);
    w2 <- which(regexpr(regEnd, model) != -1);
    if (length(w1) > 0){
        curr.expr <- c("");
        lst <- list();
        last <- "";
        known <- list();
        for (i in 1:length(model)){
            if (any(i == w1)){
                if (regexpr(regElse, model[i]) != -1){
                    curr.expr[length(curr.expr) + 1] <- sprintf("!(%s)", last);
                } else {
                    curr.expr[length(curr.expr) + 1] <- gsub(regIf, "!(\\1)", model[i]);
                    known[[length(known) + 1]] <- c(paste(paste0("(", curr.expr[-1], ")"), collapse=" && "), curr.expr[-1]);
                    curr.expr[length(curr.expr)] <- gsub(regIf, "\\1", model[i]);
                    known[[length(known) + 1]] <- c(paste(paste0("(", curr.expr[-1], ")"), collapse=" && "), curr.expr[-1]);
                }
                lst[[i]] <- "control";
            } else if (any(i == w2)){
                last <- curr.expr[length(curr.expr)];
                curr.expr <- curr.expr[seq(1, length(curr.expr) - 1)];
                lst[[i]] <- "control";
            } else {
                lst[[i]] <- curr.expr;
            }
        }
        ret <- list();
        rm <- c();
        for (i in 1:length(known)){
            mod <- c();
            for (j in 1:length(model)){
                if (identical(lst[[j]], c(""))){
                    mod[length(mod) + 1] <- model[j];
                } else {
                    i1 <- lst[[j]][-1];
                    i2 <- known[[i]][-1];
                    i3 <- i2[1:min(length(i1), length(i2))];
                    if (!identical(i2, i3)){
                        ## Find the expression to remove
                        for (k in 1:length(known)){
                            i4 <- known[[k]][-1];
                            if (identical(i4, i3)){
                                rm <- c(rm, known[[k]][1]);
                            }
                        }
                    }
                    if (identical(i1, i3)){
                        mod[length(mod) + 1] <- model[j];
                    }
                }
            }
            ret[[known[[i]][1]]] <- paste(mod, collapse="\n");
        }
        ret <- unlist(ret);
        ret <- ret[!(names(ret) %in% rm)];
        return(ret);
    } else {
        return(paste(model, collapse="\n"));
    }
}
rxExpandIfElse.slow <- NULL



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
    if (class(model) == "character" && length(model) > 1){
        vars <- model;
    } else {
        vars <- c(rxParams(model), rxState(model),
                  "podo", "t", "time", "tlast");
    }
    vars <- sapply(vars, function(x){return(rxToSymPy(x))});
    known <- c(rxSymPy.vars, vars);
    assignInMyNamespace("rxSymPy.vars", known);
    if (length(vars) == 1){
        .Jython$exec(sprintf("%s = Symbol('%s')", vars, vars));
    } else {
        .Jython$exec(sprintf("%s = symbols('%s')", paste(vars, collapse=", "), paste(vars, collapse=" ")));
    }
    return(invisible());
}

##' Setup a sympy environment that sets up the specified RxODE model.
##'
##' @param model RxODE family of objects
##' @return a boolean indicating if the environment is cached.
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPySetup <- function(model){
    setup <- rxToSymPy(model);
    rxSymPyVars(model)
    assignInMyNamespace("rxSymPy.vars", c(rxSymPy.vars, names(setup)))
    for (line in setup){
        tmp <- line;
        names(tmp) <- NULL;
        .Jython$exec(tmp);
    }
    return(invisible());
}

##' Setup sympy envirnoment if needed
##'
##' @param model RxODE lines to setup
##' @author Matthew L. Fidler
##' @keywords internal
##' @return model lines
##' @export
rxSymPySetupIf <- function(model){
    if (class(model) != "character"){
        return(rxSymPySetup(model));
    } else {
        lastLine <- sub(rex::rex(start, any_spaces, capture(anything), any_spaces, end),
                        "\\1", strsplit(model[length(model)], "=")[[1]][1])
        if (!rxSymPyExists(rxToSymPy(lastLine))){
            model.setup <- paste(model, collapse="\n");
            rxSymPySetup(model.setup)
        }
        return(model)
    }
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
        return(rxSymPySetupIf(rxSymPyDfDyFull(model, vars, rxCondition(model))));
    } else {
        if (!is.null(model)){
            rxSymPySetup(model);
        }
        var1 <- rxToSymPy(sprintf("d/dt(%s)", df));
        var <- rxToSymPy(sprintf("df(%s)/dy(%s)", df, dy));
        line <- sprintf("diff(%s,%s)", var1, rxToSymPy(dy));
        line <- rSymPy::sympy(line);
        .Jython$exec(sprintf("%s=%s", var, line));
        ret <- sprintf("df(%s)/dy(%s) = %s", df, dy, rxFromSymPy(line));
        return(ret);
    }
}

rxSymPyDfDyFull <- function(model, vars, cond){
    if (class(vars) == "logical"){
        if (vars){
            jac <- expand.grid(s1=rxState(model), s2=c(rxState(model), rxParams(model)),
                               stringsAsFactors=FALSE);
        } else  {
            jac <- expand.grid(s1=rxState(model), s2=rxState(model),
                               stringsAsFactors=FALSE);
        }
    }
    jac <- with(jac, data.frame(jac,
                                rx=sprintf("df(%s)/dy(%s)", s1, s2),
                                sym=sprintf("__d_df_%s_dy_%s__", s1, s2),
                                stringsAsFactors=FALSE))
    rxSymPySetup(model);
    extraLines <- c();
    rxCat("Calculate Jacobian.");
    for (dfdy in jac$rx){
        if (!any(dfdy == rxDfdy(model))){
            extraLines[length(extraLines) + 1] <- with(jac[jac$rx == dfdy, ], rxSymPyDfDy(NULL, s1, s2));
            rxCat(".")
        }
    }
    rxCat("done!\n");
    return(extraLines);
}

rxSymPyDfDyFull.slow <- NULL;

##' Calculate the full jacobain for a model
##'
##' This expand the model to caluclate the jacobian.  This requires
##' rSymPy.
##'
##' @param model RxODE family of objects
##' @return RxODE syntax for model with jacobian specified.
##' @author Matthew L. Fidler
rxSymPyJacobian <- function(model){
    cnd <- rxNorm(model, TRUE); ## Get the conditional statements
    extraLines <- c();
    if (is.null(cnd)){
        extraLines <- rxSymPyDfDy(model, vars=FALSE);
    } else {
        extraLines <- c();
        for (i in cnd){
            rxCat("################################################################################\n");
            rxCat(sprintf("## Calculate for %s\n", i));
            rxCat("################################################################################\n");
            rxCondition(model, i);
            extraLines <- c(extraLines,
                            sprintf("if %s {", i),
                            rxSymPyDfDy(model, vars=FALSE),
                            "}")
        }
        rxCondition(model, FALSE);
    }
    ## extraLines <- extraLines[regexpr(rex::rex("=", any_spaces, "0", any_spaces, or(";","")), extraLines) == -1];
    model <- sprintf("%s\n%s", rxNorm(model), paste(extraLines, collapse="\n"));
    rxSymPyClean();
    return(model);
}

rxSymPyJacobian.slow <- NULL

##' Does the varaible exists in the sympy python environment
##'
##' @param var Variable to test if it exists.
##' @return boolean
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPyExists <- function(var){
    return(regexpr(rex::rex("'", var, "'"), rSymPy::sympy("dir()")) != -1)
}
##' Delete variable if exists.
##'
##' @param var Variable to delete.
##' @author Matthew L. Fidler
##' @keywords internal
rxSymPyClear <- function(var){
    if (rxSymPyExists(var)){
        .Jython$exec(sprintf("del %s", var));
    }
}

rxSymPySensitivityFull.text <- "Calculate Sensitivities."

## Note the model and cond are not used in the function BUT are used
## to memoise the correct call. Please don't remove them :)
rxSymPySensitivityFull <- function(state, calcSens, model, cond){
    rxCat(rxSymPySensitivityFull.text);
    all.sens <- extraLines <- c();
    for (s1 in state){
        for (sns in calcSens){
            tmp <- c()
            vars <- c();
            for (s2 in state){
                vars <- c(vars, sprintf("rx__sens_%s_BY_%s__", s2, rxToSymPy(sns)));
                extra <- sprintf("df(%s)/dy(%s)*rx__sens_%s_BY_%s__", s1, rxToSymPy(s2), s2, rxToSymPy(sns))
                tmp <- c(tmp, extra);
            }
            tmp <- c(tmp, sprintf("df(%s)/dy(%s)", s1, rxToSymPy(sns)));
            rxSymPyVars(vars);
            all.sens <- c(all.sens, vars);
            line <- rxToSymPy(paste(tmp, collapse=" + "))
            line <- rSymPy::sympy(line);
            var.rx <- sprintf("d/dt(rx__sens_%s_BY_%s__)", s1, rxToSymPy(sns))
            var <- rxToSymPy(var.rx)
            .Jython$exec(sprintf("%s=%s", var, line));
            assignInMyNamespace("rxSymPy.vars", c(rxSymPy.vars, var));
            extraLines[length(extraLines) + 1] <- sprintf("%s=%s", var.rx, rxFromSymPy(line));
            rxCat(".")
        }
    }
    rxCat("\ndone!\n");
    return(list(all.sens=all.sens, extraLines=extraLines))
}

rxSymPySensitivityFull.slow <- NULL;


rxSymPySensitivity2Full_ <- function(state, s1, eta, sns, all.sens){
    v1 <- rxToSymPy(sprintf("d/dt(rx__sens_%s_BY_%s__)", s1, rxToSymPy(sns)));
    tmp <- c(sprintf("diff(%s,%s)",v1, rxToSymPy(eta)));
    vars <- c();
    for (s2 in state){
        extra <- sprintf("diff(%s,%s)*rx__sens_%s_BY_%s__", v1, rxToSymPy(s2), s2, rxToSymPy(eta))
        v2 <- sprintf("rx__sens_%s_BY_%s_BY_%s__", s2, rxToSymPy(eta), rxToSymPy(sns));
        vars <- c(vars, v2);
        tmp <- rxToSymPy(sprintf("df(%s)/dy(%s)", s1, rxToSymPy(s2)))
        extra2 <- sprintf("%s*%s", tmp, v2)
        tmp <- c(tmp, extra, extra2);
    }
    rxSymPyVars(vars);
    all.sens <- c(all.sens, vars);
    line <- paste(tmp, collapse=" + ");
    line <- rSymPy::sympy(line);
    var.rx <- sprintf("d/dt(rx__sens_%s_BY_%s_BY_%s__)", s1, rxToSymPy(eta), rxToSymPy(sns))
    var <- rxToSymPy(var.rx)
    .Jython$exec(sprintf("%s=%s", var, line));
    assignInMyNamespace("rxSymPy.vars", c(rxSymPy.vars, var));
    rxCat(".");
    return(list(all.sens=all.sens, line=sprintf("%s=%s", var.rx, rxFromSymPy(line))));
}

## Note the model and cond are not used in the function BUT are used
## to memoise the correct call. Please don't remove them :)
rxSymPySensitivity2Full <- function(state, etas, thetas, model, cond){
    all.sens <- extraLines <- c();
    rxCat(rxSymPySensitivityFull.text);
    for (s1 in state){
        for (eta in etas){
            if (identical(etas, thetas)){
                tmp <- rxSymPySensitivity2Full_(state, s1, eta, eta, all.sens);
            } else {
                for (sns in thetas){
                    tmp <- rxSymPySensitivity2Full_(state, s1, eta, sns, all.sens);
                }
            }
            extraLines[length(extraLines) + 1] <- tmp$line;
            all.sens <- tmp$all.sens;
        }
    }
    rxCat("\ndone!\n");
    return(list(all.sens=all.sens, extraLines=extraLines))
}
rxSymPySensitivity2Full.slow <- NULL;

rxSymPySensitivity.single <- function(model, calcSens, calcJac){
    rxSymPySetupIf(model);
    state <- rxState(model);
    extraLines <- rxSymPyDfDy(model, vars=TRUE);
    if (class(calcSens) == "list" && all(c("eta","theta") %in% names(calcSens))){
        eta <- calcSens$eta;
        theta <- calcSens$theta;
        assignInMyNamespace("rxSymPySensitivityFull.text", "Calculate d/dt(d(state)/d(eta)) .");
        on.exit({assignInMyNamespace("rxSymPySensitivityFull.text", "Calculate Sensitivites.")})
        ## Calculate dx/dn
        tmp <- rxSymPySensitivityFull(state, eta, model, rxCondition(model));
        all.sens <- tmp$all.sens;
        extraLines <- c(extraLines, rxSymPySetupIf(tmp$extraLines));
        ## Calculate dx/dT
        assignInMyNamespace("rxSymPySensitivityFull.text", "Calculate d/dt(d(state)/d(theta)) .");
        tmp <- rxSymPySensitivityFull(state, theta, model, rxCondition(model));
        all.sens <- c(all.sens, tmp$all.sens);
        extraLines <- c(extraLines, tmp$extraLines);
        ## Calculate d^2x/dx^2
        assignInMyNamespace("rxSymPySensitivityFull.text", "Calculate d/dt(d^2(state)/d(eta)^2) .");
        tmp <- rxSymPySensitivity2Full(state, eta, eta, model, rxCondition(model));
        all.sens <- c(all.sens, tmp$all.sens);
        extraLines <- c(extraLines, tmp$extraLines);
        ## Calculate d^2x/(dn dT)
        assignInMyNamespace("rxSymPySensitivityFull.text", "Calculate d/dt(d^2(state)/d(eta)d(theta)) .");
        tmp <- rxSymPySensitivity2Full(state, eta, theta, model, rxCondition(model));
        all.sens <- c(all.sens, tmp$all.sens);
        extraLines <- c(extraLines, tmp$extraLines);
    } else {
        tmp <- rxSymPySensitivityFull(state, calcSens, model, rxCondition(model))
        extraLines <- c(extraLines, rxSymPySetupIf(tmp$extraLines));
        all.sens <- tmp$all.sens;
    }
    if (calcJac){
        rxCat("Expanding Jacobian for sensitivities.")
        jac2 <- expand.grid(s1=unique(all.sens), s2=unique(c(all.sens, rxState(model))),
                            stringsAsFactors=FALSE)
        for (i in 1:length(jac2$s1)){
            extraLines[length(extraLines) + 1] <- rxSymPyDfDy(NULL, jac2$s1[i], jac2$s2[i]);
            rxCat(".");
            if (i %% 5 == 0){
                rxCat(i);
            }
            if (i %% 50 == 0){
                rxCat("\n");
            }
        }
        rxCat("\ndone!\n");
        ## extraLines <- extraLines[regexpr(rex::rex("=", any_spaces, "0", end), extraLines) == -1];
    } else {
        extraLines <- rxRmJac(extraLines);
    }
    extraLines <- extraLines[regexpr(rex::rex(any_spaces, regJac, any_spaces, "=", any_spaces,
                                              "0", any_spaces, or(";", ""), any_spaces), extraLines) == -1];
    ## cat(paste(extraLines, collapese="\n"), "\n")
    return(extraLines);
}

##' Calculate the sensitivity equations for a model
##'
##' This expands the model to calculate sensitivities.  This requires
##' rSymPy.
##'
##' @param model RxODE family of objects
##' @param calcSens Either a logical or list of sensitivity parameters
##'     to calculate. When \code{TRUE}, calculate the sensitivities of
##'     all the known parameters.  When \code{FALSE} raise an error.
##' @param calcJac A boolean that determines if the jacobian should be
##'     calculated.
##' @param keepState State parameters to keep the sensitivites for.
##'
##' @param collapseModel A boolean to collapse the model to not
##'     include any left handed quantities.
##' @return Model syntax that includes the sensitivity parameters.
##' @author Matthew L. Fidler
##' @export
rxSymPySensitivity <- function(model, calcSens, calcJac=FALSE, keepState=NULL,
                               collapseModel=FALSE){
    if (missing(calcSens)){
        calcSens <- rxParams(model);
    }
    if (class(calcSens) == "logical"){
        if (calcSens){
            calcSens <- rxParams(model);
        } else {
            stop("It is pointless to request a sensitivity calculation when calcSens=FALSE.")
        }
    }
    cnd <- rxNorm(model, TRUE); ## Get the conditional statements
    extraLines <- c();
    if (is.null(cnd)){
        rxSymPySetupIf(model);
        if (collapseModel){
            for (v in rxState(model)){
                tmp <- rSymPy::sympy(rxToSymPy(sprintf("d/dt(%s)", v)));
                tmp <- rxFromSymPy(tmp);
                extraLines[length(extraLines) + 1] <- sprintf("d/dt(%s)=%s", v, tmp);
            }
        }
        extraLines <- c(extraLines, rxSymPySensitivity.single(model, calcSens, calcJac));
    } else {
        extraLines <- c();
        for (i in cnd){
            rxCat("################################################################################\n");
            rxCat(sprintf("## Calculate for %s\n", i));
            rxCat("################################################################################\n");
            rxCondition(model, i);
            rxSymPySetupIf(model);
            tmpl <- c();
            if (collapseModel){
                for (v in rxState(model)){
                    tmp <- rSymPy::sympy(rxToSymPy(sprintf("d/dt(%s)", v)));
                    tmp <- rxFromSymPy(tmp);
                    tmpl[length(tmpl) + 1] <- sprintf("d/dt(%s)=%s", v, tmp);
                }
            }
            extraLines <- c(extraLines,
                            sprintf("if %s {", i),
                            tmpl,
                            rxSymPySensitivity.single(model, calcSens, calcJac),
                            "}")
        }
        rxCondition(model, FALSE);
    }
    if (!is.null(keepState)){
        w1 <- which(regexpr(rex::rex("d/dt(", regSens), extraLines) != -1);
        if (length(w1) > 0){
            w2 <- which(regexpr(rex::rex("d/dt(rx__sens_", or(keepState)), extraLines[w1]) == -1);
            if (length(w2) > 0){
                extraLines <- extraLines[-w1[w2]];
            }
        }
    }
    if (collapseModel){
        ret <- paste(extraLines, collapse="\n");
    } else {
        ret <- sprintf("%s\n%s", rxNorm(model), paste(extraLines, collapse="\n"));
    }
    return(ret);
}

rxSymPySensitivity.slow <- NULL;

##' Remove variables created by RxODE from the sympy environment.
##'
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPyClean <- function(){
    for (v in unique(rxSymPy.vars)){
        rxSymPyClear(v);
    }
    assignInMyNamespace("rxSymPy.vars", c());
}

##' Add a return statment to a function.
##'
##' @param fn Function to deparse
##' @param ret boolean stating if a return statement will be added.
##' @return Function with parens removed and add a return statment.
##' @author Matthew L. Fidler
rxAddReturn <- function(fn, ret=TRUE){
    txt <- deparse(body(fn));
    if (txt[1] == "{"){
        txt <- txt[-c(1, length(txt))];
    }
    ## FIXME, naieve assumption about functions.
    if (ret){
        if (regexpr(rex::rex(or(boundary, start), "return(", anything, ")"), txt[length(txt)], perl=TRUE) == -1){
            ## Add return statement
            if (regexpr(rex::rex("}"), txt[length(txt)]) == -1){
                txt[length(txt)] <- gsub(rex::rex(start, any_spaces, capture(anything), or(";", ""), any_spaces, end), "return(\\1);", txt[length(txt)]);
            }
        }
    }
    return(paste(txt, collapse="\n"))
}

rxSymPySetupDPred <- function(newmod, calcSens, states, prd="rx_pred_"){
    extraLines <- c();
    zeroSens <- TRUE;
    for (var in calcSens){
        newLine2 <- rSymPy::sympy(sprintf("diff(%s,%s)", prd, rxToSymPy(var)));
        tmp <- c(newLine2);
        for (state in states){
            newLine <- rSymPy::sympy(sprintf("diff(%s,%s)",prd, state));
            tmp[length(tmp) + 1] <- sprintf("(%s)*rx__sens_%s_BY_%s__",
                                            newLine, state, rxToSymPy(var));
        }
        tmp <- paste(paste0("(", tmp, ")"), collapse=" + ")
        tmp <- rSymPy::sympy(sprintf("simplify(%s(%s))",
                                     ifelse(prd == "rx_pred_", "-", ""),
                                     tmp
                                     ));
        if (tmp != "0"){
            zeroSens <- FALSE;
        }
        tmp <- sprintf("rx__sens_%s_BY_%s__ = %s", prd, rxToSymPy(var), rxFromSymPy(tmp));
        extraLines[length(extraLines) + 1] <- tmp;
    }
    attr(extraLines, "zeroSens") <- zeroSens;
    return(extraLines);
}

rxIf__ <- function(x){
    x <- unique(x);
    x <- x[x != ""];
    if (length(x) == 0){
        return("");
    } else if (length(x) == 1){
        return(sprintf("if %s {", x));
    } else {
        x <- sapply(x, function(x){
            return(substring(x, 2, nchar(x) - 1));
        });
        return(sprintf("if (%s) {", paste(x, collapse=" && ")));
    }
}
##' Setup Pred function based on RxODE object.
##'
##' This is for the so-called inner problem.
##'
##' @param obj RxODE object
##' @param predfn Prediction function
##' @param pkpars Pk Pars function
##' @param errfn Error function
##' @return RxODE object expanded with predfn and with calculated sensitivities.
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPySetupPred <- function(obj, predfn, pkpars, errfn){
    rxSymPyVars(obj);
    on.exit({rxSymPyClean()});
    if (!missing(pkpars)){
        txt <- rxAddReturn(pkpars, FALSE);
        newmod <- rxGetModel(paste0(paste(txt, collapse="\n"), "\n", rxNorm(obj)));
        obj <- newmod;
        etas <- rxParams(obj)
        etas <- etas[regexpr(regEta, etas) != -1];
        if (length(etas) > 0){
            calcSens <- etas;
        } else {
            calcSens <- rxParams(obj);
        }
    } else {
        calcSens <- rxParams(obj)
        newmod <- obj;
    }
    baseState <- rxState(obj);
    obj <- rxGetModel(rxNorm(obj), calcSens=calcSens, collapseModel=TRUE);
    base <- rxNorm(obj);
    txt <- rxAddReturn(predfn)
    txt <- rxParseErr(txt,ret="rx_pred_")
    cnd.base <- rxNorm(newmod, TRUE); ## Get the conditional statements
    if (is.null(cnd.base)) cnd.base <- "";
    pred.mod <- rxGetModel(txt);
    cnd.pred <- rxNorm(pred.mod, TRUE);
    if (is.null(cnd.pred)) cnd.pred <- "";
    one.pred <- FALSE;
    some.pred <- FALSE;
    setup.prd <- function(curr.base, curr.pred){
        rxCondition(obj, ifelse(curr.base == "", FALSE, curr.base));
        rxCondition(pred.mod, ifelse(curr.pred == "", FALSE, curr.pred));
        rxSymPyClean();
        rxSymPySetup(paste0(rxNorm(obj), "\n", rxNorm(pred.mod)));
        lines <- rxIf__(c(curr.base, curr.pred));
        return(lines)
    }
    grd <- expand.grid(cnd.base, cnd.pred)
    grd <- grd[!(paste0("(!",grd$Var1,")") == grd$Var2), ]
    grd <- grd[!(paste0("(!",grd$Var2,")") == grd$Var1), ]
    pred <- paste(apply(grd, 1, function(x){
        curr.base <- x[1];
        curr.pred <- x[2];
        lines <- c();
        lgl <- TRUE
        lines <- setup.prd(curr.base, curr.pred);
        if (lines == ""){
            lgl <- FALSE;
            lines <- "";
        }
        rxCat(sprintf("## Calculate d(f)/d(eta) %s\n", lines));
        newmod <- rxGetModel(paste0(rxNorm(obj), "\n", rxNorm(pred.mod)));
        newlines <- rxSymPySetupDPred(newmod, calcSens, baseState)
        if(attr(newlines, "zeroSens")){
            some.pred <<- TRUE;
        } else {
            one.pred <<- TRUE;
        }
        if (length(newlines) == 0){
            rxSymPyClean()
            stop();
        }
        lines <- c(lines, rxNorm(pred.mod), newlines);
        if (lgl){
            lines <- c(lines, "}")
        }
        rxCat(sprintf("## %sdone\n", ifelse(lgl, "}", "")));
        return(paste(lines, collapse="\n"));
    }), collapse="\n");
    if (!one.pred){
        stop("At least some part of your prediction function needs to depend on the state variables.")
    }
    if (some.pred){
        warning("Some of your prediction function does not depend on the state varibles.");
    }
    if (!missing(errfn)){
        pars <- rxParams(newmod);
        w <- which(sapply(pars, function(x){
            return(substr(x, 0, 2) == "TH")
        }))
        pars <- pars[w];
        mtheta <- max(as.numeric(gsub(rex::rex("THETA[", capture(numbers), "]"), "\\1", pars)))
        txt <- rxAddReturn(errfn);
        txt <- rxParseErr(txt, base.theta=mtheta + 1);
        extra.pars <- attr(txt, "ini");
        err.mod <- rxGetModel(txt);
        cnd.err <- rxNorm(err.mod, TRUE);
        if (is.null(cnd.err)) cnd.err <- "";
        grd <- expand.grid(cnd.base, cnd.err, cnd.pred)
        grd <- grd[!(paste0("(!",grd$Var1,")") == grd$Var2), ]
        grd <- grd[!(paste0("(!",grd$Var1,")") == grd$Var3), ]
        grd <- grd[!(paste0("(!",grd$Var2,")") == grd$Var1), ]
        grd <- grd[!(paste0("(!",grd$Var2,")") == grd$Var3), ]
        grd <- grd[!(paste0("(!",grd$Var3,")") == grd$Var1), ]
        grd <- grd[!(paste0("(!",grd$Var3,")") == grd$Var2), ]
        setup.err <- function(curr.base, curr.err, curr.pred){
            rxCondition(obj, ifelse(curr.base == "", FALSE, curr.base));
            rxCondition(pred.mod, ifelse(curr.pred == "", FALSE, curr.pred));
            rxCondition(err.mod, ifelse(curr.err == "", FALSE, curr.err));
            rxSymPyClean();
            rxSymPySetup(paste0(rxNorm(obj), "\n", rxNorm(pred.mod), "\n", rxNorm(err.mod)));
            lines <- rxIf__(c(curr.base, curr.pred, curr.err));
            return(lines)
        }
        err <- paste(apply(grd, 1, function(x){
            curr.base <- x[1];
            curr.err <- x[2];
            curr.pred <- x[3];
            lines <- c();
            lgl <- TRUE
            lines <- setup.err(curr.base, curr.err, curr.pred);
            if (lines == ""){
                lgl <- FALSE;
                lines <- c();
            }
            rxCat(sprintf("## Calculate d(err)/d(eta) %s\n", lines));
            newmod <- rxGetModel(paste0(rxNorm(obj), "\n", rxNorm(err.mod)));
            newlines <- rxSymPySetupDPred(newmod, calcSens, baseState, prd="rx_r_")
            lines <- c(lines, rxNorm(err.mod), newlines);
            if (lgl){
                lines <- c(lines, "}")
            }
            rxCat(sprintf("## %sdone\n", ifelse(lgl, "}", "")));
            return(paste(lines, collapse="\n"));
        }), collapse="\n")
    } else {
        err <- ""
    }
    ## txt <- paste0(rxNorm(obj), "\n", txt);
    return(RxODE(paste0(base, "\n", pred, "\n", err)));
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
## polygamma -> polygamma(n, z) returns log(gamma(z)).diff(n + 1) = psigamma(z, n)
## trigamma -> trigamma
