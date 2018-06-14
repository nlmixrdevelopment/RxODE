utils::globalVariables(c(".Jython"));
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
.rxRmIni <- function(x){
    x <- x[regexpr(rex::rex(start,any_spaces,or(names(rxInits(x))),any_spaces,or("=", "~")), x) == -1];
    x <- x[regexpr(rex::rex(start,any_spaces,or(names(rxInits(x))),"(0)", any_spaces,or("=", "~")), x) == -1];
    return(x);
}
##' Remove Jacobian
##' @param x RxODE list of lines to remove
##' @return RxODE lines with df/dy removed.
##' @author Matthew L. Fidler
##' @keywords internal
.rxRmJac <- function(x){
    return(x[regexpr(rex::rex(regJac, any_spaces, or("=", "~")), x) == -1])
}
##' Remove sensitivity equations
##' @param x RxODE lines to remove
##' @return Lines with d/dt(rx_sens_...._) removed.
##' @author Matthew L. Fidler
.rxRmSens <- function(x){
    return(x[regexpr(rex::rex(start, any_spaces, "d/dt(", any_spaces, regSens, any_spaces, ")", any_spaces, or("=", "~")), x) == -1]);
}
##' Remove print statements
##' @param x RxODE lines to remove
##' @return RxODE with print lines removed.
##' @author Matthew L. Fidler
.rxRmPrint <- function(x){
    return(x[regexpr(getFromNamespace("regPrint", "RxODE"), x) == -1]);
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
rxExpandIfElse <- memoise::memoise(function(model, removeInis=TRUE, removePrint=TRUE){
    ## expand if/else blocks into a list with lines for conditions that are true
    x <- strsplit(rxNorm(model, FALSE), "\n")[[1]];
    if (removeInis){
        x <- .rxRmIni(x);
    }
    if (removePrint){
        x <- .rxRmPrint(x);
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
                    known[[length(known) + 1]] <- c(paste(paste0("(", curr.expr[-1], ")"),
                                                          collapse=" && "), curr.expr[-1]);
                    curr.expr[length(curr.expr)] <- gsub(regIf, "\\1", model[i]);
                    known[[length(known) + 1]] <- c(paste(paste0("(", curr.expr[-1], ")"),
                                                          collapse=" && "), curr.expr[-1]);
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
})

.rxSymPy <- NULL;
rxSymPy.vars <- c();

##' Start SymPy
##'
##' @author G Grothendieck, Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPyStart <- function(){
    if (is.null(.rxSymPy)){
        assignInMyNamespace(".rxSymPy", new.env(parent = emptyenv()))
        .rxSymPy$started <- NULL;
    }
    start <- function(python){
        if (is.null(.rxSymPy$started) &&
            any(RxODE.sympy.engine == c("", python))){
            if (requireNamespace(python, quietly = TRUE)){
                tmp <- try({RxODE::rxSymPyExec("import sys", .python = python, .start = FALSE)});
                if (!inherits(tmp, "try-error")){
                    try({rxSymPyExec("import gc", .python = python, .start = FALSE)});
                    tmp <- try({rxSymPyExec("from sympy import *", .python=python, .start=FALSE)})
                    if (inherits(tmp, "try-error")){
                        ## try to install
                        tmp <- try({
                            system("python -m pip install sympy");
                            rxSymPyExec("from sympy import *", .python=python, .start=FALSE);
                            .rxSymPy$started <- python;
                        })
                        if (inherits(tmp, "try-error")){
                            rxCat("Could not install sympy in the system python.\n");
                            if (requireNamespace("rSymPy", quietly = TRUE)){
                                rxSymPyExec( paste( "sys.path.append(", system.file( "Lib", package = "rSymPy" ), ")", sep = '"' ),
                                            .python=python, .start=FALSE);
                                rxSymPyExec( "from sympy import *",
                                            .python=python, .start=FALSE);
                                rxCat(sprintf("Using sympy in rSymPy by running it in %s\n", python));
                                .rxSymPy$started <- python;
                            }
                        } else {
                            rxCat(sprintf("Successfully installed sympy\nUsing sympy via %s\n", python));
                            .rxSymPy$started <- python;
                        }
                    } else {
                        rxCat(sprintf("Using sympy via %s\n", python));
                        .rxSymPy$started <- python;
                    }
                }
            }
        }
    }
    start("SnakeCharmR");
    ## start("rPython");
    ## start("PythonInR");
    if (is.null(.rxSymPy$started)){
        if (any(RxODE.sympy.engine == c("", "rSymPy"))) {
            if (requireNamespace("rSymPy", quietly = TRUE)){
                if (!exists(".Jython", .GlobalEnv)){
                    rxCat("Using sympy via rSymPy (creating jython process)\n");
                    rSymPy::sympyStart()
                    .rxSymPy$started <- "rSymPy";
                    try({.Jython$exec("import gc")});
                } else {
                    rxCat("Using sympy via rSymPy (with exisiting jython)\n");
                    rSymPy::sympyStart()
                    .rxSymPy$started <- "rSymPy";
                    try({.Jython$exec("import gc")});
                }
            }
        }
    }

    if (is.null(.rxSymPy$started)){
        rxCat("RxODE requires SymPy for this function.\n");
        rxCat("We recommend you install SymPy for Python and then interact with Python using SnakeCharmR.\n");
        rxCat("In Windows you can have help setting this up by typing: `rxWinPythonSetup()`.\n");
        rxCat("Another option is to use the package rSymPy, which depends on Java and is a bit slower (and older) version of SymPy.\n");
        stop("Could not start SymPy");
    }
}
##' Execute Python statement without getting the return value.
##'
##' @param ... Parameters sent to Jython to execute.
##' @param .python Python to use, defaults to the started python.
##' @param .start A boolean (default TRUE) starting python first if
##'     RxODE doesn't know it has been started...
##' @return Nothing
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPyExec <- function(..., .python, .start=TRUE){
    if (.start && missing(.python)){
        rxSymPyStart();
    }
    if (missing(.python)){
        .python <- .rxSymPy$started;
    }
    if (.python == "SnakeCharmR"){
        SnakeCharmR::py.exec(...);
    }
    ## if (.python == "rPython"){
    ##     rPython::python.exec(...);
    ## }
    ## if (.python == "PythonInR"){
    ##     PythonInR::pyExecp(...);
    ## }
    if (.python == "rSymPy"){
        .Jython$exec(...);
    }
}
##' Fix SymPy expressions to be R parsable expressions
##'
##' @param var sympy expression
##' @return R valid expression
##' @author Matthew L. Fidler
rxSymPyFix <- function(var){
    ret <- gsub(rex::rex(boundary, "_"), "rx_underscore_", var, perl=TRUE);
    ret <- gsub(rex::rex(",", any_spaces, ")"), ")", ret)
    return(ret);
}

rxSymPyPreFix <- function(var){
    ret <- gsub(rex::rex(boundary, "rx_underscore_"), "_", var, perl=TRUE);
    return(ret);
}

##' Execute sympy statement.
##'
##' @param ... Parameters sent to Jython to execute.
##' @return String representing the return from sympy environment
##' @author G Grothendieck, Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPy <- function(...){
    ret <- try(rxSymPy0(...), silent=TRUE);
    if (inherits(ret, "try-error")){
        args <- list(...)
        cmd <- args[[1]]
        reg <- rex::rex("diff(", capture(except_any_of(",")), ",", capture(except_any_of(")")), ")")
        if (regexpr(reg, cmd) != -1){
            var <- gsub(reg, "\\1", cmd);
            val <- try(rxSymPy0(var), silent=TRUE);
            if (!inherits(val, "try-error")){
                stop(sprintf("Error running SymPy command:\n    %s\n    %s: %s\n\n%s",
                             cmd, var, val, attr(ret, "condition")$message))
            } else {
                stop(sprintf("Error running SymPy command:\n    %s\n    %s: %s\n\n%s",
                             cmd, var, attr(val, "condition")$message, attr(ret, "condition")$message))

            }
        }
        stop(sprintf("Error running sympy command:\n    %s\n\n%s", cmd, attr(ret, "condition")$message));
    } else {
        ret <- paste(strsplit(ret, "\n+")[[1]], collapse=" ");
    }
    return(ret)
}

rxSymPy0 <- function(...){
    rxSymPyStart();
    if (.rxSymPy$started == "SnakeCharmR"){
        SnakeCharmR::py.exec(paste("__Rsympy=None"))
        SnakeCharmR::py.exec(paste("__Rsympy=", ..., sep = ""))
        SnakeCharmR::py.exec(paste("__Rsympy = str(__Rsympy)"))
        ret <- SnakeCharmR::py.get("__Rsympy");
        return(rxSymPyFix(ret));
    }
    ## if (.rxSymPy$started == "rPython"){
    ##     rPython::python.exec(paste("__Rsympy=None"))
    ##     rPython::python.exec(paste("__Rsympy=", ..., sep = ""))
    ##     rPython::python.exec(paste("__Rsympy = str(__Rsympy)"))
    ##     ret <- rPython::python.get("__Rsympy");
    ##     return(rxSymPyFix(ret));
    ## }
    ## if (.rxSymPy$started == "PythonInR"){
    ##     PythonInR::pyExec(paste("__Rsympy=None"))
    ##     PythonInR::pyExec(paste("__Rsympy=", ..., sep = ""))
    ##     PythonInR::pyExec(paste("__Rsympy = str(__Rsympy)"))
    ##     ret <- PythonInR::pyGet("__Rsympy");
    ##     return(rxSymPyFix(ret));
    ## }
    if (.rxSymPy$started == "rSymPy"){
        ret <- rxSymPyFix(rSymPy::sympy(...))
        return(ret);
    }
}
##' Return the version of SymPy that is running
##'
##' @param numeric boolean that specifies if the major and minor
##'     release should be a number.
##' @return Verson of sympy that is running.
##' @author Matthew L. Fidler
##' @export
rxSymPyVersion <- memoise::memoise(function(numeric=TRUE){
    rxSymPyExec("import sympy");
    ret <- rxSymPy("sympy.__version__");
    if (numeric)
        ret <- as.numeric(sub("^([0-9]+[.][0-9]+).*", "\\1", ret))
    return(ret)
})

##' Return a list of reserved functions and variables from sympy
##'
##' @return List of reserevd functions and variaibles from sympy.
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPyReserved <- memoise::memoise(function(){
    rxSymPyStart();
    rxSymPyExec("import sympy");
    vars <- rxSymPy("dir(sympy)")
    vars <- eval(parse(text=sprintf("c(%s)", substr(vars, 2, nchar(vars) - 1))));
    return(vars)
})

##' Add undefined variables to SymPy
##'
##' @param txt Text to add undefined variables to...
##' @return nothing
##' @author Matthew L. Fidler
##' @export
##' @keywords internal
rxSyPyAddVars <- function(txt){
    allVars <- function(x){
        defined <- character()
        f <- function(x){
            if (is.atomic(x)) {
                character()
            } else if (is.name(x)) {
                as.character(x)
            } else if (is.call(x) || is.pairlist(x)) {
                if (identical(x[[1]], quote(`~`)) ||
                    identical(x[[1]], quote(`=`)) ||
                    identical(x[[1]], quote(`<-`))){
                    if (is.call(x[[3]])){
                        ret <- unique(unlist(lapply(x[[3]][-1], f)));
                    } else {
                        ret <- unique(unlist(lapply(x[[3]], f)));
                    }
                    ret <- ret[!(ret %in% defined)]
                    defined <<- unique(c(defined, x[[2]]))
                    return(ret)
                } else {
                    children <- lapply(x[-1], f)
                    unique(unlist(children))
                }
            } else {
                stop("Don't know how to handle type ", typeof(x), ".",
                     call. = FALSE)
            }
        }
        f(x);
    }
    vars <- allVars(eval(parse(text=sprintf("quote({%s})", txt))))
    for (v in vars){
        if (!rxSymPyExists(v)){
            known <- c(rxSymPy.vars, v);
            assignInMyNamespace("rxSymPy.vars", known);
            if (length(vars) == 1){
                rxSymPyExec(sprintf("%s = Symbol('%s')", v, v));
            }
        }
    }
    return(invisible());
}

##' Setup SymPy variables
##'
##' This creates SymPy variables for later evaulation in the CAS SymPy.
##'
##' @param model RxODE family of objects
##' @return NULL
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPyVars <- function(model){
    rxSymPyStart();
    if (rxIs(model,"character") && length(model) > 1){
        vars <- model;
    } else if (rxIs(model,"character") && length(model) == 1 && regexpr(rex::rex(or("=", "<-", "~")), model) == -1){
        vars <- model;
    } else {
        vars <- c(rxParams(model),
                  rxState(model),
                  "podo", "t", "time", "tlast",
                  "rx__PTR__", "rx1c");
    }
    vars <- sapply(vars, function(x){return(rxToSymPy(x))});
    known <- c(rxSymPy.vars, vars);
    assignInMyNamespace("rxSymPy.vars", known);
    if (length(vars) == 1){
        rxSymPyExec(sprintf("%s = Symbol('%s')", vars, vars));
    } else {
        rxSymPyExec(sprintf("%s = symbols('%s')", paste(vars, collapse=", "), paste(vars, collapse=" ")));
    }
    return(invisible());
}

##' Setup SymPy functions
##'
##' This creates SymPy unspecified functions for later evaulation in the CAS SymPy.
##'
##' @param functions a list of functions
##' @return NULL
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPyFunctions <- function(functions){
    vars <- sapply(functions, function(x){return(rxToSymPy(x))});
    known <- c(rxSymPy.vars, vars);
    assignInMyNamespace("rxSymPy.vars", known);
    for (f in vars){
        rxSymPyExec(sprintf("%s = Function('%s')", vars, vars));
        rxSymPyFEnv[[f]] <- functionOp(f);
        sympyRxFEnv[[f]] <- functionOp(f);
    }
    return(invisible());
}


##' Setup a SymPy environment that sets up the specified RxODE model.
##'
##' @param model RxODE family of objects
##' @param envir Environment to evaluate the model in (defaults to
##'     parent frame)
##' @return a boolean indicating if the environment is cached.
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPySetup <- function(model, envir=parent.frame()){
    if (identical(model, "")) return(invisible());
    setup <- rxToSymPy(model, envir=envir);
    const <- rxInits(model, rxLines=TRUE);
    if (!identical(const, "")) setup <- c(rxToSymPy(const, envir=envir), setup);
    rxSymPyVars(model)
    assignInMyNamespace("rxSymPy.vars", rxSymPy.vars)
    for (line in c(setup)){
        tmp <- line;
        names(tmp) <- NULL;
        rxSymPyExec(tmp);
    }
    return(invisible());
}

##' Setup SymPy environment if needed
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
                        "\\1", strsplit(model[length(model)], "[=~]")[[1]][1])
        if (!rxSymPyExists(rxToSymPy(lastLine))){
            model.setup <- paste(model, collapse="\n");
            rxSymPySetup(model.setup)
        }
        return(model)
    }
}
##' Split line into multiple lines at + or - breaks
##'
##' @param lhs Left handed side to construct
##' @param rhs Right handed side to construct
##' @param limit the number of characters for the expression to be
##'     before it is split.  By default this is 1100
##' @return an expression where the lhs is constructed iteratevly by
##'     splitting the lhs and adding it iteratevly to the lhs.
##'
##' For example:
##'
##' lhs = 1 + 2 + 3
##'
##' Would become
##'
##' lhs = 1
##' lhs = lhs + 2
##' lhs = lhs + 3
##'
##' When calling rxSplitLines("lhs", "1+2+3", 0)
##'
##' This is to deal with the unwieldly lines that sometimes come out
##' of SymPy.
##'
##' @author Matthew L. Fidler
##' @export
##' @keywords internal
rxSplitLines <- function(lhs, rhs, limit=1100){
    if ((nchar(lhs) + nchar(rhs) + 4) < limit){
        return(sprintf("%s = %s", lhs, rhs))
    } else {
        rh <- eval(parse(text=sprintf("rxSplitPlusQ(quote(%s))", rhs)));
        rh[-1] <- gsub("[+][-]", "-", paste0(lhs, "+", rh[-1]));
        return(paste(paste0(lhs, " = ", rh), collapse="\n"));
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
##'
##'     vars can also be a list of variables.  In this case the
##'     df(.)/dy(.) will be augmented with the variables in this list.
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
        line <- rxSymPy(line);
        rxSymPyExec(sprintf("%s=%s", var, line));
        ret <- rxSplitLines(sprintf("df(%s)/dy(%s)", df, dy), rxFromSymPy(line));
        return(ret);
    }
}

rxSymPyDfDyFull <- memoise::memoise(function(model, vars, cond){
    if (rxIs(vars,"logical")){
        if (vars){
            jac <- expand.grid(s1=rxState(model), s2=c(rxState(model), rxParams(model, FALSE)),
                               stringsAsFactors=FALSE);
        } else  {
            jac <- expand.grid(s1=rxState(model), s2=rxState(model),
                               stringsAsFactors=FALSE);
        }
    } else if (rxIs(vars,"character")){
        jac <- expand.grid(s1=rxState(model), s2=c(rxState(model), vars),
                           stringsAsFactors=FALSE)

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
    rxCat("done.\n");
    return(extraLines);
})

##' Calculate the full Jacobian for a model
##'
##' This expand the model to caluclate the Jacobian.  This requires
##' rSymPy.
##'
##' @param model RxODE family of objects
##' @return RxODE syntax for model with Jacobian specified.
##' @author Matthew L. Fidler
.rxSymPyJacobian <- memoise::memoise(function(model){
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
})

##' Does the varaible exist in the SymPy Python environment?
##'
##' @param var Variable to be tested.
##' @return boolean
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPyExists <- function(var){
    return(regexpr(rex::rex("'", var, "'"), rxSymPy("dir()")) != -1)
}
##' Delete variable if exists.
##'
##' @param var Variable to delete.
##' @author Matthew L. Fidler
##' @keywords internal
rxSymPyClear <- function(var){
    if (rxSymPyExists(var)){
        try(rxSymPyExec(sprintf("del %s", var)), silent=TRUE);
    }
    tmp <- rxToSymPy(var);
    if (rxSymPyExists(tmp)){
        try(rxSymPyExec(sprintf("del %s", tmp)), silent=TRUE)
    }
}

rxSymPySensitivityFull.text <- "Calculate sensitivities."

## Note the model and cond are not used in the function BUT are used
## to memoise the correct call. Please don't remove them :)
rxSymPySensitivityFull <- memoise::memoise(function(state, calcSens, model, cond){
    rxCat(rxSymPySensitivityFull.text);
    all.sens <- extraLines <- c();
    for (s1 in state){
        for (sns in calcSens){
            tmp <- c()
            vars <- c();
            for (s2 in state){
                v <- sprintf("rx__sens_%s_BY_%s__", s2, rxToSymPy(sns))
                if (!rxSymPyExists(v)){
                    rxSymPyExec(sprintf("%s=0", v));
                }
                vars <- c(vars, v);
                extra <- sprintf("df(%s)/dy(%s)*rx__sens_%s_BY_%s__", s1, rxToSymPy(s2), s2, rxToSymPy(sns))
                tmp <- c(tmp, extra);
            }
            v <- sprintf("df(%s)/dy(%s)", s1, rxToSymPy(sns))
            v2 <- rxToSymPy(v);
            if (!rxSymPyExists(v2)){
                rxSymPyExec(sprintf("%s=0", v2));
            }
            tmp <- c(tmp, v);
            rxSymPyVars(vars);
            all.sens <- c(all.sens, vars);
            line <- rxToSymPy(paste(tmp, collapse=" + "))
            line <- rxSymPy(line);
            var.rx <- sprintf("d/dt(rx__sens_%s_BY_%s__)", s1, rxToSymPy(sns))
            var <- rxToSymPy(var.rx)
            rxSymPyExec(sprintf("%s=%s", var, line));
            assignInMyNamespace("rxSymPy.vars", c(rxSymPy.vars, var));
            extraLines[length(extraLines) + 1] <- sprintf("%s=%s", var.rx, rxFromSymPy(line));
            ini <- sprintf("%s(0)", s1);
            ini <- rxToSymPy(ini)
            if (any(rxSymPy.vars == ini)){
                line <- sprintf("diff(%s, %s)", ini, rxToSymPy(sns));
                line <- rxSymPy(line);
                ## Some diffs are paritally implemented in RxODE translation.
                line <- rxFromSymPy(line)
                line <- rxToSymPy(line)
                tmp <- sprintf("rx__sens_%s_BY_%s__", s1, rxToSymPy(sns));
                if (sprintf("%s(0)", v) != rxFromSymPy(tmp))
                    extraLines[length(extraLines) + 1] <- sprintf("%s(0)=%s", tmp, rxFromSymPy(line));
            }
            rxCat(".")
        }
    }
    rxCat("\ndone.\n");
    return(list(all.sens=all.sens, extraLines=extraLines))
})

rxSymPySensitivity2Full_ <- function(state, s1, eta, sns, all.sens){
    v1 <- rxToSymPy(sprintf("d/dt(rx__sens_%s_BY_%s__)", s1, rxToSymPy(sns)));
    v1 <- rxSymPy(v1);
    tmp <- c(sprintf("diff(%s,%s)",v1, rxToSymPy(eta)));
    ## Some diffs are paritally implemented in RxODE translation.
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
    line <- rxSymPy(line);
    var.rx <- sprintf("d/dt(rx__sens_%s_BY_%s_BY_%s__)", s1, rxToSymPy(eta), rxToSymPy(sns))
    var <- rxToSymPy(var.rx)
    rxSymPyExec(sprintf("%s=%s", var, line));
    assignInMyNamespace("rxSymPy.vars", c(rxSymPy.vars, var));
    ini <- sprintf("%s(0)", s1);
    ini <- rxToSymPy(ini)
    if (any(rxSymPy.vars == ini)){
        line <- sprintf("diff(diff(%s, %s),%s)", ini, rxToSymPy(eta), rxToSymPy(sns));
        line <- rxSymPy(line);
        tmp <- sprintf("rx__sens_%s_BY_%s_BY_%s__", s1, rxToSymPy(eta), rxToSymPy(sns));
        if (paste0(tmp, "(0)") != rxFromSymPy(line)){
            ini.line <- sprintf("%s(0)=%s", tmp, rxFromSymPy(line));
        }
    } else {
        ini.line <- NULL;
    }
    rxCat(".");
    return(list(all.sens=all.sens, line=sprintf("%s=%s", var.rx, rxFromSymPy(line)), ini.line=ini.line));
}

## Note the model and cond are not used in the function BUT are used
## to memoise the correct call. Please don't remove them :)
rxSymPySensitivity2Full <- memoise::memoise(function(state, etas, thetas, model, cond){
    all.sens <- extraLines <- c();
    rxCat(rxSymPySensitivityFull.text);
    for (s1 in state){
        if (identical(etas, thetas)){
            for (e1 in etas){
                for (e2 in etas){
                    tmp2 <- sort(c(e1, e2));
                    tmpO <- as.numeric(gsub(rex::rex(start, "ETA", or("[", "_"), capture(numbers), or("]")), "\\1", tmp2));
                    tmp2 <- tmp2[order(tmpO)];
                    if (identical(tmp2, c(e1, e2))){
                        tmp <- rxSymPySensitivity2Full_(state, s1, e1, e2, all.sens);
                        all.sens <- tmp$all.sens;
                        extraLines[length(extraLines) + 1] <- tmp$line;
                        if (!is.null(tmp$ini.line)){
                            extraLines[length(extraLines) + 1] <- tmp$ini.line;
                        }
                    }
                }
            }
        } else {
            for (sns in thetas){
                for (eta in etas){
                    tmp <- rxSymPySensitivity2Full_(state, s1, eta, sns, all.sens);
                    all.sens <- unique(c(all.sens, tmp$all.sens));
                    extraLines[length(extraLines) + 1] <- tmp$line;
                    if (!is.null(tmp$ini.line)){
                        extraLines[length(extraLines) + 1] <- tmp$ini.line;
                    }
                }
            }
        }
    }
    rxCat("\ndone.\n");
    return(list(all.sens=all.sens, extraLines=extraLines))
})

rxSymPySensitivity.single <- function(model, calcSens, calcJac){
    rxSymPySetupIf(model);
    state <- rxState(model);
    state <- state[regexpr(rex::rex(start, "rx_"), state) == -1];
    if (rxIs(calcSens,"list") && all(c("eta","theta") %in% names(calcSens))){
        extraLines <- rxSymPyDfDy(model, vars=c(calcSens$eta, calcSens$theta));
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
        extraLines <- .rxRmJac(extraLines);
    } else {
        extraLines <- rxSymPyDfDy(model, vars=TRUE);
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
        rxCat("\ndone.\n");
        ## extraLines <- extraLines[regexpr(rex::rex("=", any_spaces, "0", end), extraLines) == -1];
    } else {
        extraLines <- .rxRmJac(extraLines);
    }
    extraLines <- extraLines[regexpr(rex::rex(any_spaces, regJac, any_spaces, or("=", "~"), any_spaces,
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
##' @param calcJac A boolean that determines if the Jacobian should be
##'     calculated.
##' @param keepState State parameters to keep the sensitivites for.
##'
##' @param collapseModel A boolean to collapse the model that each
##'     expression only depends on the unspecified parameters (instead on LHS quantities).
##'
##' @return Model syntax that includes the sensitivity parameters.
##' @author Matthew L. Fidler
##' @export
rxSymPySensitivity <- memoise::memoise(function(model, calcSens, calcJac=FALSE, keepState=NULL,
                               collapseModel=FALSE){
    if (missing(calcSens)){
        calcSens <- rxParams(model, FALSE);
    }
    if (rxIs(calcSens,"logical")){
        if (calcSens){
            calcSens <- rxParams(model, FALSE);
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
                tmp <- rxSymPy(rxToSymPy(sprintf("d/dt(%s)", v)));
                tmp <- rxFromSymPy(tmp);
                extraLines[length(extraLines) + 1] <- sprintf("d/dt(%s)=%s", v, tmp);
                ini <- sprintf("%s(0)", v);
                ini <- rxToSymPy(ini)
                tmp <- try({rxSymPy(ini)}, silent=TRUE);
                if(!inherits(tmp, "try-error")){
                    known <- c(rxSymPy.vars, ini);
                    assignInMyNamespace("rxSymPy.vars", known);
                    tmp <- rxSymPy(ini);
                    if (paste0(v, "(0)") != rxFromSymPy(tmp))
                        extraLines[length(extraLines) + 1] <- sprintf("%s(0)=%s", v, rxFromSymPy(tmp));
                } else if (any(v == names(rxInits(model)))){
                    tmp <- as.vector(rxInits(model)[v]);
                    if (paste0(v, "(0)") != tmp)
                        extraLines[length(extraLines) + 1] <- sprintf("%s(0)=%s", v, tmp);
                }
            }
            for (v in rxLhs(model)){
                tmp <- rxSymPy(rxToSymPy(sprintf("%s", v)));
                tmp <- rxFromSymPy(tmp);
                extraLines[length(extraLines) + 1] <- sprintf("%s=%s", v, tmp);
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
                    tmp <- rxSymPy(rxToSymPy(sprintf("d/dt(%s)", v)));
                    tmp <- rxFromSymPy(tmp);
                    tmpl[length(tmpl) + 1] <- sprintf("d/dt(%s)=%s", v, tmp);
                }
                for (v in rxLhs(model)){
                    tmp <- rxSymPy(rxToSymPy(sprintf("%s", v)));
                    tmp <- rxFromSymPy(tmp);
                    tmpl[length(tmpl) + 1] <- sprintf("%s=%s", v, tmp);
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
})

##' Remove variables created by RxODE from the SymPy environment.
##'
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPyClean <- function(){
    tmp <- rxSymPy("dir()");
    tmp <- eval(parse(text=sprintf("c(%s)", substr(tmp,2,nchar(tmp)-1))))
    tmp <- tmp[regexpr(rex::rex(start, "rx_"), tmp) != -1]
    for (v in unique(c(rxSymPy.vars, tmp))){
        rxSymPyClear(v);
    }
    ## rxSymPy("clear_cache()");
    assignInMyNamespace("rxSymPy.vars", c());
    rxGc();
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
            if (regexpr(rex::rex(or("=", "~", "<-", "}")), txt[length(txt)]) == -1){
                txt[length(txt)] <- gsub(rex::rex(start, any_spaces, capture(anything), or(";", ""), any_spaces, end), "return(\\1);", txt[length(txt)]);
            }
        }
    }
    return(paste(txt, collapse="\n"))
}

rxSymPySetupDPred <- function(newmod, calcSens, states, prd="rx_pred_", pred.minus.dv=TRUE){
    states <- states[regexpr(rex::rex(start, "rx_"), states) == -1]
    extraLines <- c();
    zeroSens <- TRUE;
    if (rxIs(calcSens,"list")){
        tmp1 <- rxSymPySetupDPred(newmod, calcSens$eta, states, prd)
        zeroSens <- attr(tmp1, "zeroSens");
        tmp2 <- rxSymPySetupDPred(newmod, calcSens$theta, states, prd)
        if (!zeroSens){
            zeroSens <- attr(tmp1, "zeroSens");
        }
        extraLines <- c(tmp1, tmp2);
        ## These derivations are given in Equation #37.
        ## Should be dErr^2/(dEtaL dEtaK)
        rxCat("##   d(.)^2/(d(eta1)d(eta2)) ");
        for (eL in calcSens$eta){
            for (eK in calcSens$eta){
                tmp2 <- c(eK, eL)
                tmpO <- as.numeric(gsub(rex::rex(start, "ETA", or("[", "_"), capture(numbers), or("]")), "\\1", tmp2));
                tmp2 <- tmp2[order(tmpO)];
                if (identical(tmp2, c(eK, eL))){
                    ## diff(h)/d(eK,eL)
                    newLine2 <- rxSymPy(sprintf("diff(diff(%s,%s),%s)", prd, rxToSymPy(eK), rxToSymPy(eL)));
                    ## Some diffs are paritally implemented in RxODE translation.
                    newLine2 <- rxFromSymPy(newLine2)
                    newLine2 <- rxToSymPy(newLine2)
                    tmp <- c(newLine2);
                    for (state in states){
                        ## diff(h)/d(eK, state) * dstate/deta_eL
                        ## = diff(h)/d(eK, state) * sensitivity_state_eL
                        newLine <- rxSymPy(sprintf("diff(diff(%s,%s),%s)",prd, rxToSymPy(eK), state));
                        tmp[length(tmp) + 1] <- sprintf("(%s)*rx__sens_%s_BY_%s__",
                                                        newLine, state, rxToSymPy(eL));

                        ## d^2(h)/d(x,eL) * dx/d(eK)
                        ## d^2(h)/d(x,eL) * sens_state_eK
                        newLine <- rxSymPy(sprintf("diff(diff(%s,%s),%s)",prd, state, rxToSymPy(eL)));
                        tmp[length(tmp) + 1] <- sprintf("(%s)*rx__sens_%s_BY_%s__",
                                                        newLine, state, rxToSymPy(eK));

                        ## d^2h/d^2(state) * (dstate/deta_el) * (dstate/deta_eK)
                        ## d^2h/d^2(state) * sens_state_el * sens_state_eK
                        newLine <- rxSymPy(sprintf("diff(diff(%s,%s),%s)",prd, state, state));
                        tmp[length(tmp) + 1] <- sprintf("(%s)*rx__sens_%s_BY_%s__*rx__sens_%s_BY_%s__",
                                                        newLine, state, rxToSymPy(eL), state, rxToSymPy(eK));
                        ## dh/d(state) * dstate/(eK, eL)
                        ## dh/d(state) * sens_state_eK_eL
                        newLine <- rxSymPy(sprintf("diff(%s,%s)",prd, state));
                        tmp[length(tmp) + 1] <- sprintf("(%s)*rx__sens_%s_BY_%s_BY_%s__",
                                                        newLine, state, rxToSymPy(eK), rxToSymPy(eL));
                    }
                    tmp <- paste(paste0("(", tmp, ")"), collapse=" + ")
                    tmp <- sprintf("%s", tmp);
                    if (!pred.minus.dv && prd == "rx_pred_"){
                        tmp <- paste0("-(", tmp, ")");
                    }
                    tmp <- rxSymPy(tmp);
                    if (tmp != "0"){
                        zeroSens <- FALSE;
                    }
                    tmp <- rxSplitLines(sprintf("rx__sens_%s_BY_%s_BY_%s__ ", prd, rxToSymPy(eK), rxToSymPy(eL)),
                                        rxFromSymPy(tmp));
                    extraLines[length(extraLines) + 1] <- tmp;
                    rxCat(".");
                }
            }
        }
        rxCat("\n##   d(.)^2/(d(theta)d(eta))");
        for (theta in calcSens$theta){
            for (eta in calcSens$eta){
                ## dh/d(eta,theta)
                newLine2 <- rxSymPy(sprintf("diff(diff(%s,%s),%s)", prd, rxToSymPy(eta), rxToSymPy(theta)));
                tmp <- c(newLine2);
                for (state in states){
                    ## d^2h/d(eta,state)*dstate/theta
                    ## d^2h/d(eta,state)*sens_state_theta
                    newLine <- rxSymPy(sprintf("diff(diff(%s,%s),%s)",prd, rxToSymPy(eta), state));
                    tmp[length(tmp) + 1] <- sprintf("(%s)*rx__sens_%s_BY_%s__",
                                                    newLine, state, rxToSymPy(theta));
                    ## d^2h/d(state,theta)*dstate/eta
                    ## d^2h/d(state,theta)*sens_state_eta
                    newLine <- rxSymPy(sprintf("diff(diff(%s,%s),%s)",prd, state, rxToSymPy(theta)));
                    tmp[length(tmp) + 1] <- sprintf("(%s)*rx__sens_%s_BY_%s__",
                                                    newLine, state, rxToSymPy(eta));
                    ## d^2h/d^2(state)*dstate/eta
                    ## d^2h/d^2(state)*sens_state_eta
                    newLine <- rxSymPy(sprintf("diff(diff(%s,%s),%s)",prd, state, state));
                    tmp[length(tmp) + 1] <- sprintf("(%s)*rx__sens_%s_BY_%s__*rx__sens_%s_BY_%s__",
                                                    newLine, state, rxToSymPy(theta), state, rxToSymPy(eta));
                    ## dh/d(state)*d(h)/(deta,dtheta)
                    newLine <- rxSymPy(sprintf("diff(%s,%s)",prd, state));
                    tmp[length(tmp) + 1] <- sprintf("(%s)*rx__sens_%s_BY_%s_BY_%s__",
                                                    newLine, state, rxToSymPy(eta), rxToSymPy(theta));

                }
                tmp <- paste(paste0("(", tmp, ")"), collapse=" + ")
                tmp <- sprintf("%s", tmp);
                if (!pred.minus.dv && prd == "rx_pred_") {
                    tmp <- paste0("-(", tmp, ")");
                }
                tmp <- rxSymPy(tmp);
                if (tmp != "0"){
                    zeroSens <- FALSE;
                }
                tmp <- rxSplitLines(sprintf("rx__sens_%s_BY_%s_BY_%s__", prd, rxToSymPy(eta), rxToSymPy(theta)), rxFromSymPy(tmp));
                extraLines[length(extraLines) + 1] <- tmp;
                rxCat(".");
            }
        }
        rxCat("\n");
        attr(extraLines, "zeroSens") <- zeroSens;
        return(extraLines);
    }
    rxCat("## ");

    ## Equation #21
    rxSymPyVars(calcSens);
    for (var in calcSens){
        ##
        newLine2 <- rxSymPy(sprintf("diff(%s,%s)", prd, rxToSymPy(var)));
        tmp <- c(newLine2);
        for (state in states){
            newLine <- rxSymPy(sprintf("diff(%s,%s)",prd, state));
            tmp[length(tmp) + 1] <- sprintf("(%s)*rx__sens_%s_BY_%s__",
                                            newLine, state, rxToSymPy(var));
        }
        tmp <- paste(paste0("(", tmp, ")"), collapse=" + ")
        if (!pred.minus.dv && prd == "rx_pred_"){
            tmp <- sprintf("-(%s)", tmp);
        }
        tmp <- rxSymPy(tmp);
        ## print(tmp)
        if (tmp != "0"){
            zeroSens <- FALSE;
        }
        tmp <- rxSplitLines(sprintf("rx__sens_%s_BY_%s__", prd, rxToSymPy(var)), rxFromSymPy(tmp));
        extraLines[length(extraLines) + 1] <- tmp;
        rxCat(".");
    }
    rxCat("\n");
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


genCmt0 <- function(ncmt=1, oral=FALSE){
    ## The lincmt function generates:
    ## 1 cmt: rx_k
    ## 2 cmt: rx_k12, rx_k21
    ## 3 cmt: rx_k13, rx_k31
    rx0 <- "";
    rxc <- "d/dt(rx1) = -rx_k*rx1";
    rxp <- "";
    rx3 <- ""
    if (ncmt >= 2){
        rxc <- paste(rxc, "- rx_k12*rx1 + rx_k21*rx2")
        rxp <- "d/dt(rx2) = rx_k12*rx1 - rx_k21*rx2";
    }
    if (ncmt == 3){
        rxc <- paste(rxc, "- rx_k13*rx1 + rx_k31*rx3")
        rx3 <- "d/dt(rx3) = rx_k13*rx1 - rx_k31*rx3";
    }
    if (oral){
        rxc <- paste(rxc, "+ rx_ka*rx0");
        rx0 <- "d/dt(rx0) = -rx_ka*rx0";
    }
    fin <- "rx1c = rx1/rx_v";
    ret <- c(rx0, rxc, rxp, rx3, fin);
    ret <- ret[ret != ""];
    paste(ret, collapse="\n")
}

genCmtMod <- function(mod){
    if (regexpr(rex::rex("solveLinB("), rxNorm(mod)) == -1) return(mod);
    ## Generates based on what is currently on the sympy stack.
    rxSymPySetup(mod);
    on.exit(rxSymPyClean());
    ret <- NULL
    get.var <- function(v){
        if (rxSymPyExists(v)){
            tmp1 <- rxSymPy(v);
            tmp1 <- rxFromSymPy(tmp1);
            return(sprintf("%s ~ %s", v, tmp1));
        } else {
            return(NULL);
        }
    }
    oral <- rxSymPyExists("rx_ka");
    if (oral){
        oral <- (rxSymPy("rx_ka") != "0")
    }
    tlag <- rxSymPyExists("rx_tlag");
    if (tlag){
        tlag <- (rxSymPy("rx_tlag") != "0")
    }
    if (tlag){
        stop("tlag not supported yet.");
    }
    if (rxSymPyExists("rx_k13")){
        extra <- genCmt0(3, oral);
    } else if (rxSymPyExists("rx_k12")){
        extra <- genCmt0(2, oral);
    } else if (rxSymPyExists("rx_k")){
        extra <- genCmt0(1, oral);
    } else {
        return(mod);
    }
    ## Now build model
    mv.1 <- rxModelVars(mod);
    orig.state <- mv.1$state
    orig.state.ignore <- mv.1$state.ignore
    ka <- NULL;
    if (oral)
        ka <- get.var("rx_ka")
    ret <- paste(c(get.var("rx_v"),
                   ka,
                   get.var("rx_k"),
                   get.var("rx_k13"),
                   get.var("rx_k31"),
                   get.var("rx_k12"),
                   get.var("rx_k21"),
                   extra,
                   sapply(seq_along(orig.state), function(i){
                       cur.state <- orig.state[i];
                       sep <- ifelse(orig.state.ignore[i] == 1L, "~", "=");
                       v <- rxSymPy(rxToSymPy(sprintf("d/dt(%s)", cur.state)));
                       v <- rxFromSymPy(v);
                       return(sprintf("d/dt(%s) %s %s", cur.state, sep, v));
                   }),
                   sapply(mv.1$lhs, function(v){
                       v1 <- rxSymPy(v);
                       v1 <- rxFromSymPy(v1);
                       return(sprintf("%s=%s", v, v1));
                   })), collapse="\n")
    ret
}

rxSymPySetupPred.warn <- FALSE
##' Setup Pred function based on RxODE object.
##'
##' This is for the so-called inner problem.
##'
##' @param obj RxODE object
##' @param predfn Prediction function
##' @param pkpars Pk Pars function
##' @param errfn Error function
##' @param init Initialization parameters for scaling.
##' @param grad Boolaen indicated if the the equations for the
##'     gradient be calculated
##' @param sum.prod A boolean determining if RxODE should use more
##'     numerically stable sums/products.
##' @param pred.minus.dv Boolean stating if the FOCEi objective
##'     function is based on PRED-DV (like NONMEM).  Default TRUE.
##' @param theta.derivs Boolean indicating if theta derivatives are
##'     setup
##' @param only.numeric Instead of setting up the sensitivities for
##'     the inner problem, modify the RxODE to use numeric
##'     differentiation for the numeric inner problem only.
##' @param grad.internal Internal gradient flag.  This function is
##'     recursively called, and this shouldn't be set by the user.
##' @param theta.internal Internal theta flag.  This function is
##'     recursively called and shouldn't be called by the user.
##' @param run.internal Boolean to see if the function should be run
##'     internally.
##' @return RxODE object expanded with predfn and with calculated
##'     sensitivities.
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
##' @importFrom utils find
rxSymPySetupPred <- function(obj, predfn, pkpars=NULL, errfn=NULL, init=NULL, grad=FALSE, sum.prod=FALSE, pred.minus.dv=TRUE,
                             theta.derivs=FALSE,only.numeric=FALSE,
                             grad.internal=FALSE, theta.internal=FALSE,
                             run.internal=RxODE.sympy.run.internal){
    good.fns <- c(".GlobalEnv", "package:RxODE", "package:nlmixr")
    check.good <- function(x){
        tmp <- suppressWarnings({find(deparse(substitute(x)))})
        if (!identical(tmp, character())){
            if (!any(tmp == good.fns)){
                stop(sprintf("%s is from %s and can't be used in this context.", deparse(substitute(x)), tmp))
            }
        }
    }
    check.good(predfn);
    check.good(pkpars);
    check.good(errfn);
    ##
    if (!grad.internal && !theta.internal){
        cache.file <- file.path(ifelse(RxODE.cache.directory == ".", getwd(), RxODE.cache.directory),
                                sprintf("rx_%s%s.prd",
                                        digest::digest(paste(deparse(list(rxModelVars(obj)$md5["parsed_md5"],
                                                                          ifelse(is.function(predfn), paste(deparse(body(predfn)), collapse=""), ""),
                                                                          ifelse(is.function(pkpars), paste(deparse(body(pkpars)), collapse=""), ""),
                                                                          ifelse(is.function(errfn), paste(deparse(body(errfn)), collapse=""), ""),
                                                                          init, grad, sum.prod, pred.minus.dv, theta.derivs, only.numeric)), collapse="")),
                                        .Platform$dynlib.ext));
    } else {
        cache.file <- "";
    }
    if (file.exists(cache.file) && !grad.internal && !theta.internal){
        load(file=cache.file);
        if (any(names(ret) == "warn")){
            if (ret$warn){
                warning("Some of your prediction function does not depend on the state variables.");
            }
        }
        if (!is.null(ret$inner)){
            rxLoad(ret$inner);
        }
        if (!is.null(ret$outer)){
            rxLoad(ret$outer);
        }
        if (!is.null(ret$theta)){
            rxLoad(ret$theta);
        }
        return(ret);
    } else {
        if (!run.internal && !grad.internal && !theta.internal){
            rfile <- tempfile(fileext=".R")
            dta <- tempfile(fileext=".rdata")
            save(obj, predfn, pkpars, errfn, init, grad, sum.prod, pred.minus.dv, theta.derivs, only.numeric, file=dta);
            ## on.exit({unlink(dta); unlink(rfile)});
            tmp <- options();
            tmp <- tmp[regexpr(rex::rex(start, "RxODE."), names(tmp)) != -1];
            rf <- sprintf("%s;options(RxODE.delete.unnamed=FALSE);load(%s); require(RxODE);tmp <- rxSymPySetupPred(obj, predfn, pkpars, errfn, init, grad, sum.prod, pred.minus.dv, theta.derivs, only.numeric, run.internal=TRUE);",
                          sub("options\\(\\)", "", sub(rex::rex(",", any_spaces, ".Names", anything,end),"",sub(rex::rex(start,"structure(list("),"options(",paste0(deparse(tmp),collapse="")))),
                          deparse(dta));
            rf <- gsub("^;", "", rf)
            sink(rfile);
            cat(rf);
            sink();
            sh <- "system"   # windows's default shell COMSPEC does not handle UNC paths
            cmd <- sprintf("%s/bin/Rscript %s",
                           Sys.getenv("R_HOME"), deparse(rfile));
            suppressWarnings(system(cmd, ignore.stdout=!RxODE.verbose,
                                    ignore.stderr=!RxODE.verbose));
            if (file.exists(cache.file)){
                load(file=cache.file);
                rxLoad(ret$inner);
                if (!is.null(ret$outer)){
                    rxLoad(ret$outer);
                }
                if (!is.null(ret$theta)){
                    rxLoad(ret$theta);
                }
                if (any(names(ret) == "warn")){
                    if (ret$warn){
                        warning("Some of your prediction function does not depend on the state variables.");
                    }
                }
                return(ret);
            } else {
                stop("Error Creating expanded function.");
            }
        } else {
            assignInMyNamespace("rxSymPyExpThetas", c());
            assignInMyNamespace("rxSymPyExpEtas", c());
            gobj <- obj;
            oobj <- genCmtMod(obj);
            obj <- oobj;
            rxModelVars(oobj)
            rxSymPyVars(obj);
            on.exit({rxSymPyClean()});
            lhs <- c();
            if (!is.null(pkpars)){
                txt <- as.vector(unlist(strsplit(rxParsePk(pkpars, init=init), "\n")));
                re <- rex::rex(start, "init", or("_", ".", ""), or("C", "c"), "ond", or("ition", ""),or("", "s"), any_spaces,
                               or("=", "<-", "~"), any_spaces, "c(", capture(anything), ")", any_spaces, ";", any_spaces, end);
                w <- which(regexpr(re, txt) != -1)
                if (length(w) == 1){
                    inis <- txt[w];
                    inis <- strsplit(gsub(re, "\\1", inis), " *, *")[[1]];
                    if (length(rxState(oobj)) == length(inis)){
                        inis <- paste(paste0(rxState(oobj), "(0)=", inis, ";"), collapse="\n");
                        txt[w] <- inis;
                    } else {
                        stop("Specified %s initial conditions when there are only %s states.", length(inis), length(rxState(oobj)));
                    }
                } else if (length(w) > 1){
                    stop("Multiple initCondition= found.  Please choose one.");
                }
                newmod <- rxGetModel(paste0(paste(txt, collapse="\n"), "\n", rxNorm(obj)));
                obj <- newmod;
                lhs <- c(rxLhs(newmod), rxLhs(oobj));
            } else {
                calcSens <- rxParams(obj, FALSE)
                newmod <- obj;
            }
            txt <- rxParsePred(predfn, init=init)
            pred.mod <- rxGetModel(txt);
            extra.pars <- c();
            if (!is.null(errfn)){
                pars <- rxParams(rxGetModel(paste0(rxNorm(obj), "\n", rxNorm(pred.mod))), FALSE);
                w <- which(sapply(pars, function(x){
                    return(substr(x, 0, 2) == "TH")
                }))
                pars <- pars[w];
                mtheta <- max(as.numeric(gsub(rex::rex("THETA[", capture(numbers), "]"), "\\1", pars)))
                txt <- rxParseErr(errfn, base.theta=mtheta + 1, init=init);
                extra.pars <- attr(txt, "ini");
                err.mod <- rxGetModel(txt);
            }
            if (!is.null(pkpars)){
                full <- paste0(rxNorm(obj), "\n", rxNorm(pred.mod));
                if (!is.null(errfn)){
                    full <- paste0(full, "\n", rxNorm(err.mod));
                }
                full <- rxGetModel(full);
                etas <- rxParams(full, FALSE);
                thetas <- etas[regexpr(regTheta, etas) != -1];
                etas <- etas[regexpr(regEta, etas) != -1];
                if (only.numeric){
                    calcSens <- c()
                } else if (length(etas) > 0 && !theta.internal && !grad.internal){
                    rxCat("## Calculate ETA-based prediction and error derivatives:\n")
                    calcSens <- etas;
                } else {
                    calcSens <- rxParams(full, FALSE);
                }
                if (grad.internal){
                    rxCat("## Calculate THETA/ETA-based prediction and error 1st and 2nd order derivatives:\n")
                    calcSens <- list(eta=etas);
                    calcSens$theta <- thetas;
                }
                if (theta.internal){
                    rxCat("## Calculate THETA-based prediction and error derivatives:\n")
                    calcSens <- thetas;
                }
            }
            baseState <- rxState(obj);
            oLhs <- rxLhs(obj);
            if (length(baseState) > 0){
                obj <- rxGetModel(rxNorm(obj, removeJac=TRUE, removeSens=TRUE), calcSens=calcSens, collapseModel=TRUE);
            }
            base <- rxNorm(obj);
            cnd.base <- rxNorm(newmod, TRUE); ## Get the conditional statements
            if (is.null(cnd.base)) cnd.base <- "";
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
                newmod <- rxGetModel(paste0(rxNorm(obj), "\n", rxNorm(pred.mod)));
                if (only.numeric){
                    rxSymPySetup(newmod);
                    tmp <- rxSymPy("rx_pred_");
                    prd <- gsub(rex::rex("rx_pred_=", anything), sprintf("rx_pred_=%s", rxFromSymPy(tmp)), rxNorm(pred.mod));
                    lines <- c(lines, prd)
                    states <- rxState(newmod);
                    reg <- rex::rex("rx_pred_=", anything, or(states), anything);
                    if (regexpr(reg, prd) != -1){
                        one.pred <<- TRUE
                    }
                } else {
                    rxCat(sprintf("## Calculate d(f)/d(%seta) %s\n", ifelse(theta.internal, "th", ""), lines));
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
                    newmod <- rxGetModel(rxNorm(newmod));
                    tmp <- rxSymPy("rx_pred_");
                    lines <- c(lines,
                               ## Make sure rx_pred_ is not in terms of LHS components.
                               gsub(rex::rex("rx_pred_=", anything), sprintf("rx_pred_=%s", rxFromSymPy(tmp)), rxNorm(pred.mod)),
                               newlines);
                    if (lgl){
                        lines <- c(lines, "}")
                    }
                    rxCat(sprintf("## %sdone\n", ifelse(lgl, "}", "")));
                }
                return(paste(lines, collapse="\n"));
            }), collapse="\n");
            if (!one.pred){
                rxCat(pred, "\n");
                stop("At least some part of your prediction function needs to depend on the state variables.")
            }
            if (some.pred){
                rxSymPySetupPred.warn <- TRUE
            }
            if (!is.null(errfn)){
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
                    newmod <- rxGetModel(paste0(rxNorm(obj), "\n", rxNorm(err.mod)));
                    if (only.numeric){
                        rxSymPySetup(newmod);
                        tmp <- rxSymPy("rx_r_");
                        prd <- gsub(rex::rex("rx_r_=", anything),
                                    sprintf("rx_r_=%s", rxFromSymPy(tmp)),
                                    rxNorm(err.mod));
                        lines <- c(lines, prd)
                        if (lgl){
                            lines <- c(lines, "}")
                        }
                    } else {
                        rxCat(sprintf("## Calculate d(err)/d(%seta) %s\n", ifelse(theta.internal, "th", ""), lines));
                        newlines <- rxSymPySetupDPred(newmod, calcSens, baseState, prd="rx_r_")
                        tmp <- rxSymPy("rx_r_");
                        lines <- c(lines,
                                   ## Make sure rx_r_ is not in terms of LHS components.
                                   gsub(rex::rex("rx_r_=", anything), sprintf("rx_r_=%s", rxFromSymPy(tmp)), rxNorm(err.mod)),
                                   newlines);
                        if (lgl){
                            lines <- c(lines, "}")
                        }
                        rxCat(sprintf("## %sdone\n", ifelse(lgl, "}", "")));
                    }
                    return(paste(lines, collapse="\n"));
                }), collapse="\n")
            } else {
                err <- ""
            }
            if (!grad.internal && !theta.internal){
                outer <- NULL;
                theta <- NULL;
                if (grad){
                    outer <- rxSymPySetupPred(obj=gobj, predfn=predfn, pkpars=pkpars,
                                              errfn=errfn, init=init, sum.prod=sum.prod,pred.minus.dv=pred.minus.dv,
                                              run.internal=TRUE, grad.internal=TRUE, theta.internal=FALSE);
                }
                base <- gsub(rex::rex(or(group(anything, "~", anything), group(or(oLhs), "=", anything, ";"))), "", strsplit(base, "\n")[[1]])
                base <- paste(base[base != ""], collapse="\n");
                mod <- rxGetModel(paste0(base, "\n", pred, "\n", err))
                ostate <- rxState(oobj);
                if (only.numeric){
                    mod <- strsplit(strsplit(rxNorm(mod), "\n")[[1]], "=");
                    mod <- sapply(mod, function(x){
                        tmp <- x[1];
                        var <- rxToSymPy(tmp)
                        var <- rxSymPy(var)
                        var <- rxFromSymPy(var);
                        if (regexpr(rex::rex("d/dt("), x[1]) != -1){
                            return(sprintf("%s=%s;", x[1], var));
                        } else {
                            return(sprintf("%s=%s;", x[1], gsub(";+ *$", "", x[2])));
                        }
                    })
                    mod <- paste(mod, collapse="\n");
                }
                if (sum.prod){
                    mod <- rxSumProdModel(mod);
                    base <- rxSumProdModel(base)
                }
                mod <- rxNorm(mod);
                rxGc();
                if (!grad && theta.derivs && !only.numeric){
                    rxForget();
                    theta <- rxSymPySetupPred(obj=gobj, predfn=predfn, pkpars=pkpars,
                                              errfn=errfn, init=init, sum.prod=sum.prod,pred.minus.dv=pred.minus.dv,
                                              theta.derivs=TRUE, run.internal=TRUE, grad.internal=FALSE, theta.internal=TRUE);
                    theta <- RxODE(rxNorm(theta));
                }
                keep <- c(sprintf("d/dt(%s)", ostate), sprintf("%s(0)=", ostate), "rx_pred_=", "rx_r_=");
                pred.only <- strsplit(rxNorm(mod), "\n")[[1]]
                pred.only <- paste(pred.only[regexpr(rex::rex(start, or(keep)), pred.only) != -1], collapse="\n");
                keep <- c(sprintf("d/dt(%s)", ostate), sprintf("%s(0)=", ostate));
                r <- c(do.call(`c`, lapply(baseState, function(x){
                    ini <- sprintf("%s(0)", x);
                    ddt <- sprintf("d/dt(%s)", x);
                    c(rxToSymPy(ini), rxToSymPy(ddt));
                    })), sapply(lhs, function(x){return(rxToSymPy(x))}))
                ebe <- paste(do.call(`c`, lapply(r, function(x){
                                      if (rxSymPyExists(x)){
                                          v <- rxSymPy(x)
                                          return(sprintf("%s=%s;", rxFromSymPy(x), rxFromSymPy(v)))
                                      }
                                      return(NULL)
                                      })), collapse="\n")
                if (only.numeric){
                    mod <- NULL;
                } else {
                    ## Inner should hide compartments
                    mod <- gsub(rex::rex(capture("d/dt(", except_some_of(")"), ")"), "="), "\\1~", rxNorm(mod), perl=TRUE);
                    mod <- RxODE(mod);
                }
                ret <- list(obj=oobj,
                            pred.only=RxODE(pred.only),
                            ebe=RxODE(ebe),
                            inner=mod,
                            extra.pars=extra.pars,
                            outer=outer,
                            theta=theta,
                            warn=rxSymPySetupPred.warn,
                            pred.minus.dv=pred.minus.dv,
                            log.thetas=rxSymPyExpThetas,
                            log.etas=rxSymPyExpEtas,
                            cache.file=cache.file);
                rxSymPySetupPred.warn <- FALSE;
                class(ret) <- "rxFocei";

                save(ret, file=cache.file);
                if (only.numeric){
                    rxCat("The model has been setup to calculate residuals.\n");
                } else {
                    rxCat(sprintf("The model-based sensitivities have been calculated%s.\n", ifelse(grad, " (with FOCEi Global Gradient)", "")));
                }
                if (ret$warn){
                    warning("Some of your prediction function does not depend on the state variables.");
                }
                return(ret);
            } else {
                base <- gsub(rex::rex(or(group(anything, "~", anything), group(or(oLhs), "=", anything, ";"))), "", strsplit(base, "\n")[[1]])
                base <- paste(base[base != ""], collapse="\n");
                mod <- rxGetModel(paste0(base, "\n", pred, "\n", err));
                if (sum.prod){
                    mod <- rxSumProdModel(mod);
                }
                mod <- rxNorm(mod);
                return(RxODE(mod));
            }
        }
    }
}

rxGc <- function(){
    ## gc(reset=TRUE);
    tf <- tempfile();
    sink(tf);
    on.exit({sink();unlink(tf)})
    try({rxSymPyExec("gc.collect()")});
}

##' Return the expanded expression (via SymPy)
##'
##' @param x text SymPy-compatible expression
##' @param expr SymPy expression to use.  By default uses expand.
##' @return SymPy value for rxSymPy...
##' @author Matthew L. Fidler
##' @export
##' @keywords internal
rxSymPyExpand <- function(x, expr="expand"){
    rxSyPyAddVars(x)
    return(rxSymPy(sprintf("%s(%s)", expr, rxToSymPy(x))))
}


##' Takes a model and expands it to log multiplication
##'
##' This is a numerical trick to help reduce multiplcation errors when
##' numbers are of very different magnitude.
##'
##' @param model RxODE model
##' @param expand Expand the symbolic expression using SymPy?
##' @return Lines for expanded model. (but not compiled)
##' @author Matthew L. Fidler
##' @export
rxLogifyModel <- function(model, expand=TRUE){
    rxSymPyVars(model);
    ## rxCat(rxNorm(model));
    on.exit({
        rxSymPyClean();
        assignInMyNamespace("rxSymPyAbsLog", FALSE);
        assignInMyNamespace("rxSymPyLogSign", c());
    });
    n <- 1;
    rxCat("## Logify ie a*b -> exp(log(abs(a))+log(abs(b)))*sign(a*b)\n")
    rxCat("## ")
    negReg <- rex::rex(start, any_spaces, "-", any_spaces)
    cnd <- rxNorm(model, TRUE);
    lines <- strsplit(rxNorm(model), "\n")[[1]];
    for (i in seq_along(lines)){
        if (regexpr("[=~]", lines[i])){
            l0 <- strsplit(lines[i], "[=~]")[[1]];
            if (length(l0) == 2){l1 <- l0[1];
                l2 <- eval(parse(text=sprintf("rxSplitPlusQ(quote(%s))", substr(l0[2], 1, nchar(l0[2]) - 1))));
                for (j in seq_along(l2)){
                    tmp <- l2[j];
                    neg <- FALSE;
                    if (regexpr(negReg, tmp) != -1){
                        tmp <- gsub(negReg, "", tmp)
                        neg <- TRUE;
                    }
                    tmp2 <- eval(parse(text=sprintf("rxSimpleExprP(quote(%s))", tmp)));
                    if (!tmp2){
                        tmp <- rxSymPy(sprintf(paste0("expand_log(log(", ifelse(expand, "expand(%s)", "%s"), "), force=True)"), rxToSymPy(tmp)));
                        assignInMyNamespace("rxSymPyAbsLog", TRUE);
                        assignInMyNamespace("rxSymPyLogSign", c());
                        tmp <- rxFromSymPy(tmp);
                        if (length(rxSymPyLogSign) == 0){
                            tmp <- sprintf("exp(%s)", tmp);
                        } else {
                            even.v <- c();
                            zero.v <- c();
                            for (v in rxSymPyLogSign){
                                regs <- c()
                                regs[1] <- rex::rex(start, capture(any_spaces), capture(any_numbers), any_spaces,
                                                    "*", any_spaces, "abs_log(", v, ")", anything)
                                regs[2] <- rex::rex(anything, capture(or("+", "-")), any_spaces, capture(any_numbers), any_spaces,
                                                    "*", any_spaces, "abs_log(", v, ")", anything);
                                regs[3] <- rex::rex(anything, capture("-"), capture(any_spaces), "abs_log(", v, ")", anything);
                                for (reg in regs){
                                    if (regexpr(reg, tmp, perl=TRUE) != -1){
                                        num <- gsub(reg, "\\2", tmp, perl=TRUE);
                                        pm <- gsub(reg, "\\1", tmp, perl=TRUE);
                                        num <- suppressWarnings({as.numeric(num)});
                                        if (!is.na(num)){
                                            if (num %% 2 == 0){
                                                even.v <- c(even.v, v);
                                            }
                                        }
                                        if (pm == "-"){
                                            ## Should be zero protected...
                                            zero.v <- c(zero.v, v);
                                        }
                                    }
                                }
                            }
                            vs <- rxSymPyLogSign
                            if (length(even.v) > 0) {
                                vs <- vs[!(vs %in% even.v)]
                            }
                            if (length(zero.v)){
                                vs[vs %in% zero.v] <- sprintf("safe_zero(%s)", vs[vs %in% zero.v])
                            }
                            for (v in zero.v){
                                tmp <- gsub(rex::rex("abs_log(", v, ")"),
                                            sprintf("abs_log(safe_zero(%s))", v), tmp);
                            }
                            if (length(vs) > 0){
                                tmp <- sprintf("sign_exp(%s, %s)",
                                               paste(vs, collapse="*"),
                                               tmp);
                            } else {
                                tmp <- sprintf("exp(%s)", tmp);
                            }
                        }
                    }
                    if (neg){
                        tmp <- sprintf("-%s", tmp);
                    }
                    rxCat(".");
                    if (n %% 5 == 0){
                        rxCat(n);
                    }
                    if (n %% 50 == 0){
                        rxCat("\n## ");
                    }
                    n <- n + 1;
                    l2[j] <- tmp;
                }
                l0 <- rxSplitLines(l1, gsub(rex::rex("+-"), "-", paste(l2, collapse="+")));
                ## tmp <- gsub(rex::rex(any_spaces, "=", any_spaces, "+", any_spaces), "=",
                ##             gsub(rex::rex(any_spaces, "-", any_spaces, "+", any_spaces), "-",
                ##                  gsub(rex::rex(any_spaces, "+", any_spaces, "-", any_spaces), "-",
                ##                       paste(l1, "=", c("", rep(l1, length(l2) - 1)), "+", l2))))
                ## l0 <- paste(tmp, collapse="\n");
                lines[i] <- l0;
            }
        }
    }
    rxCat("\n## done.\n");
    newMod <- paste(paste(lines, collapse="\n"), "\n");
    ## rxCat(newMod)
    return(newMod);
}
## Supported SymPy special functions
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
