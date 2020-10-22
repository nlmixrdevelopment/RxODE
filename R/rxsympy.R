utils::globalVariables(c(".Jython", "add", "id", "sim.id"));
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
##' @return RxODE lines with \code{df/dy} removed.
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

##' Expand if/else clauses into multiple different types of lines.
##'
##'
##' @param model Model can be a character, or a RxODE model.  It needs
##'     to have normalized syntax, that is \code{if (...)\{} has to be
##'     on the same line.  The \code{else} statement must be on its
##'     own line with the closing bracket of the \code{if} statement
##'     on the previous line.  This \code{else} statement must also
##'     contain the opening bracket, like the code \code{else \{}
##' @param removeInis A boolean indicating if parameter
##'     initialization should be removed from the model.
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
        for (i in seq_along(model)){
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
        for (i in seq_along(known)){
            mod <- c();
            for (j in seq_along(model)){
                if (identical(lst[[j]], c(""))){
                    mod[length(mod) + 1] <- model[j];
                } else {
                    i1 <- lst[[j]][-1];
                    i2 <- known[[i]][-1];
                    i3 <- i2[seq(1,min(length(i1), length(i2)))];
                    if (!identical(i2, i3)){
                        ## Find the expression to remove
                        for (k in seq_along(known)){
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
    start("reticulate");
    if (is.null(.rxSymPy$started)){
        rxCat("RxODE requires SymPy for this function.\n");
        rxCat("We recommend you install SymPy for Python and then interact with Python using reticulate.\n");
        rxCat("In Windows you can have help setting this up by typing: `rxWinPythonSetup()`.\n");
        rxCat("Another option is to use SnakeCharmR\n");
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
    if (.python == "reticulate"){
        fns <- paste0(...);
        sapply(fns, function(x){
            reticulate::py_run_string(code=x, convert=FALSE)
        });
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
    } else if (.rxSymPy$started != "reticulate"){
        ret <- paste(strsplit(ret, "\n+")[[1]], collapse=" ");
    }
    return(ret)
}

rxSymPy0 <- function(...){
    rxSymPyStart();
    if (.rxSymPy$started == "reticulate"){
        ## reticulate::py_run_string("__Rsympy=None",convert=FALSE)
        ret <- reticulate::py_run_string(paste0("__Rsympy=", paste(...)),convert=FALSE)
        ret <- try(reticulate::py_to_r(ret$`__Rsympy`), silent=TRUE);
        if (inherits(ret, "try-error")){
            ret <- reticulate::py_run_string(paste0("__Rsympy=str(", paste(...), ")"),convert=FALSE);
            ret <- reticulate::py_to_r(ret$`__Rsympy`);
            return(rxSymPyFix(ret));
        } else {
            return(rxSymPyFix(ret));
        }
    }
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
}
##' Return the version of SymPy that is running
##'
##' @param numeric boolean that specifies if the major and minor
##'     release should be a number.
##' @return Version of sympy that is running.
##' @author Matthew L. Fidler
##' @export
rxSymPyVersion <- function(numeric=TRUE){
    rxSymPyExec("import sympy");
    ret <- rxSymPy("sympy.__version__");
    if (numeric)
        ret <- as.numeric(sub("^([0-9]+[.][0-9]+).*", "\\1", ret))
    return(ret)
}

##' Return a list of reserved functions and variables from sympy
##'
##' @return List of reserved functions and variables from sympy.
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPyReserved <- function(){
    rxSymPyStart();
    rxSymPyExec("import sympy");
    vars <- rxSymPy("dir(sympy)")
    if (.rxSymPy$started != "reticulate"){
        vars <- eval(parse(text=sprintf("c(%s)", substr(vars, 2, nchar(vars) - 1))));
    }
    return(c(vars, "lambda", "N", "E1"))
}

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
##' This creates SymPy variables for later evaluation in the CAS SymPy.
##'
##' @param model RxODE family of objects
##' @return NULL
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPyVars <- function(model){
    if (!is.null(model)){
        rxSymPyStart();
        if (rxIs(model,"character") && length(model) > 1){
            vars <- model;
        } else if (rxIs(model,"character") && length(model) == 1 && regexpr(rex::rex(or("=", "<-", "~")), model) == -1){
            vars <- model;
        } else {
            vars <- c(.rxParams(model),
                      rxState(model),
                      "podo", "t", "time", "tlast",
                      "rx__PTR__", "rx1c",
                      "rx_lambda_", "rx_yj_",
                      sprintf("rx_underscore_xi_%s", 1:100));
        }
        vars <- sapply(vars, function(x){return(rxToSymPy(x))});
        known <- c(rxSymPy.vars, vars);
        assignInMyNamespace("rxSymPy.vars", known);
        if (length(vars) == 1){
            rxSymPyExec(sprintf("%s = Symbol('%s')", vars, vars));
        } else {
            rxSymPyExec(sprintf("%s = symbols('%s')", paste(vars, collapse=", "), paste(vars, collapse=" ")));
        }
    }
    return(invisible());
}

##' Setup SymPy functions
##'
##' This creates SymPy unspecified functions for later evaluation in the CAS SymPy.
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
    rxSymPyFunctions(c("rxTBS", "rxTBSd"));
    setup <- rxToSymPy(model, envir=envir);
    const <- rxInits(model, rxLines=TRUE);
    if (!identical(const, "")) setup <- c(rxToSymPy(const, envir=envir), setup);
    rxSymPyVars(model)
    ## rxSymPyVars(c("rx_lambda_", "rx_yj_"))
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
    if (any(class(model) == c("RxODE", "rxModelVars", "rxModelText"))){
        return(rxSymPySetup(model));
    } else if (class(model) == "character") {
        lastLine <- sub(rex::rex(start, any_spaces, capture(anything), any_spaces, end),
                        "\\1", strsplit(model[length(model)], "[=~]")[[1]][1])
        if (!rxSymPyExists(rxToSymPy(lastLine))){
            model.setup <- paste(model, collapse="\n");
            rxSymPySetup(model.setup)
        }
        return(model)
    } else {
        return("");
    }
}
##' Split line into multiple lines at + or - breaks
##'
##' @param lhs Left handed side to construct
##' @param rhs Right handed side to construct
##' @param limit the number of characters for the expression to be
##'     before it is split.  By default this is 1100
##' @return an expression where the lhs is constructed iteratively by
##'     splitting the lhs and adding it iteratively to the lhs.
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
##' This is to deal with the unwieldy lines that sometimes come out
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
##'     calculated according to the \code{vars} parameter below.
##' @param dy is a string for the state or variable in the df(.)/dy(.).
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

rxSymPyDfDyFull <- function(model, vars, cond){
    if (rxIs(vars,"logical")){
        if (vars){
            .pars  <- .rxParams(model, FALSE);
            if (any(.pars=="ETA[1]")){
                .pars  <- .pars[regexpr(rex::rex(start,"ETA[",any_numbers,"]"), .pars) != -1]
            }
            jac <- expand.grid(s1=rxState(model), s2=c(rxState(model), .pars),
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
    rxCat("Calculate Jacobian\n");
    rxProgress(length(jac$rx));
    on.exit({rxProgressAbort()});
    for (dfdy in jac$rx){
        if (!any(dfdy == rxDfdy(model))){
            extraLines[length(extraLines) + 1] <- with(jac[jac$rx == dfdy, ], rxSymPyDfDy(NULL, s1, s2));
        }
        rxTick();
    }
    rxProgressStop();
    return(extraLines);
}

##' Calculate the full Jacobian for a model
##'
##' This expand the model to calculate the Jacobian.  This requires
##' rSymPy.
##'
##' @param model RxODE family of objects
##' @return RxODE syntax for model with Jacobian specified.
##' @author Matthew L. Fidler
.rxSymPyJacobian <- function(model){
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

##' Does the variable exist in the SymPy Python environment?
##'
##' @param var Variable to be tested.
##' @return boolean
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPyExists <- function(var){
    if (.rxSymPy$started == "reticulate"){
        return(any(var == rxSymPy("dir()")))
    } else {
        return(regexpr(rex::rex("'", var, "'"), rxSymPy("dir()")) != -1)
    }
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

rxSymPySensitivityFull.text <- "Calculate Sensitivities"

## Note the model and cond are not used in the function BUT are used
## to memoise the correct call. Please don't remove them :)
rxSymPySensitivityFull <- function(state, calcSens, model, cond){
    rxCat(rxSymPySensitivityFull.text, "\n");
    all.sens <- extraLines <- c();
    if (length(calcSens) == 0L){
        calcSens <- .rxParams(model);
    }
    if (length(calcSens) == 0L){
        stop("Can not calculate sensitivities with no parameters.")
    }
    rxProgress(length(state) + length(calcSens));
    on.exit({rxProgressAbort()});
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
                if (sprintf("%s(0)", v) != rxFromSymPy(tmp)){
                    ## +0 makes everything a conditional initialization
                    .tmp2 <- rxFromSymPy(line);
                    if (.tmp2 != "0"){
                        extraLines[length(extraLines) + 1] <- sprintf("%s(0)=%s+0.0", tmp, rxFromSymPy(line));
                    }
                }
            }
            rxTick();
        }
    }
    rxProgressStop()
    return(list(all.sens=all.sens, extraLines=extraLines))
}

rxSymPySensitivity2Full_ <- function(state, s1, eta, sns, all.sens){
    v1 <- rxToSymPy(sprintf("d/dt(rx__sens_%s_BY_%s__)", s1, rxToSymPy(sns)));
    v1 <- rxSymPy(v1);
    tmp <- c(sprintf("diff(%s,%s)",v1, rxToSymPy(eta)));
    ## Some diffs are paritally implemented in RxODE translation.
    vars <- c();
    for (s2 in state){
        extra <- sprintf("diff(%s,%s)*rx__sens_%s_BY_%s__", v1, rxToSymPy(s2), rxToSymPy(s2), rxToSymPy(eta))
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
            .tmp <- rxFromSymPy(line);
            if (.tmp != "0"){
                ini.line <- sprintf("%s(0)=%s+0.0", tmp, rxFromSymPy(line));
            }
        }
    } else {
        ini.line <- NULL;
    }
    rxCat(".");
    return(list(all.sens=all.sens, line=sprintf("%s=%s", var.rx, rxFromSymPy(line)), ini.line=ini.line));
}

## Note the model and cond are not used in the function BUT are used
## to memoise the correct call. Please don't remove them :)
rxSymPySensitivity2Full <- function(state, etas, thetas, model, cond){
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
}

rxSymPySensitivity.single <- function(model, calcSens, calcJac){
    rxSymPySetupIf(model);
    state <- rxState(model);
    state <- state[regexpr(rex::rex(start, "rx_"), state) == -1];
    if (rxIs(calcSens,"list") && all(c("eta","theta") %in% names(calcSens))){
        extraLines <- rxSymPyDfDy(model, vars=c(calcSens$eta, calcSens$theta));
        eta <- calcSens$eta;
        theta <- calcSens$theta;
        assignInMyNamespace("rxSymPySensitivityFull.text", "Calculate d/dt(d(state)/d(eta)) .");
        on.exit({assignInMyNamespace("rxSymPySensitivityFull.text", "Calculate Sensitivities.")})
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
        for (i in seq_along(jac2$s1)){
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
##' @param keepState State parameters to keep the sensitivities for.
##'
##' @param collapseModel A boolean to collapse the model that each
##'     expression only depends on the unspecified parameters (instead on LHS quantities).
##'
##' @return Model syntax that includes the sensitivity parameters.
##' @author Matthew L. Fidler
##' @export
rxSymPySensitivity <- function(model, calcSens, calcJac=FALSE, keepState=NULL,
                               collapseModel=FALSE){
    if (missing(calcSens)){
        calcSens <- .rxParams(model, FALSE);
    }
    if (rxIs(calcSens,"logical")){
        if (calcSens){
            calcSens <- .rxParams(model, FALSE);
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
                    if (paste0(v, "(0)") != rxFromSymPy(tmp)){
                        .tmp2 <- rxFromSymPy(tmp);
                        if (.tmp2 != "0"){
                            extraLines[length(extraLines) + 1] <- sprintf("%s(0)=%s+0.0", v, rxFromSymPy(tmp));
                        }
                    }
                } else if (any(v == names(rxInits(model)))){
                    tmp <- as.vector(rxInits(model)[v]);
                    if (paste0(v, "(0)") != tmp){
                        .tmp2 <- rxFromSymPy(tmp);
                        if (.tmp2 != "0"){
                            extraLines[length(extraLines) + 1] <- sprintf("%s(0)=%s+0.0", v, tmp);
                        }
                    }
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
}

##' Remove variables created by RxODE from the SymPy environment.
##'
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPyClean <- function(full=FALSE){
    if (.rxSymPy$started == "reticulate"){
        sapply(rxSymPy.vars, function(x){
            try(rxSymPyExec(sprintf("del %s", x)), silent=TRUE);
        })
        tmp <- rxSymPy("dir()");
        sapply(tmp, function(x){
            if (nchar(x) > 3 && substr(x, 0, 3) == "rx_"){
                try(rxSymPyExec(sprintf("del %s", x)), silent=TRUE);
            }
        })
    } else {
        tmp <- rxSymPy("dir()");
        tmp <- eval(parse(text=sprintf("c(%s)", substr(tmp,2,nchar(tmp)-1))))
        tmp <- tmp[regexpr(rex::rex(start, "rx_"), tmp) != -1]
        for (v in unique(c(rxSymPy.vars, tmp))){
            rxSymPyClear(v);
        }
    }

    ## rxSymPy("clear_cache()");
    assignInMyNamespace("rxSymPy.vars", c());
    if (full) rxGc();
}

##' Add a return statement to a function.
##'
##' @param fn Function to deparse
##' @param ret boolean stating if a return statement will be added.
##' @return Function with parens removed and add a return statement.
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
        rxCat("   d(.)^2/(d(eta1)d(eta2)) \n");
        rxProgress(length(calcSens$eta) * 2);
        on.exit({rxProgressAbort()});
        for (eL in calcSens$eta){
            for (eK in calcSens$eta){
                tmp2 <- c(eK, eL)
                tmpO <- as.numeric(gsub(rex::rex(start, "ETA", or("[", "_"), capture(numbers), or("]")), "\\1", tmp2));
                tmp2 <- tmp2[order(tmpO)];
                if (identical(tmp2, c(eK, eL))){
                    ## diff(h)/d(eK,eL)
                    newLine2 <- rxSymPy(sprintf("diff(diff(%s,%s),%s)", rxToSymPy(prd), rxToSymPy(eK), rxToSymPy(eL)));
                    ## Some diffs are paritally implemented in RxODE translation.
                    newLine2 <- rxFromSymPy(newLine2)
                    newLine2 <- rxToSymPy(newLine2)
                    tmp <- c(newLine2);
                    for (state in states){
                        ## diff(h)/d(eK, state) * dstate/deta_eL
                        ## = diff(h)/d(eK, state) * sensitivity_state_eL
                        newLine <- rxSymPy(sprintf("diff(diff(%s,%s),%s)",prd, rxToSymPy(eK), rxToSymPy(state)));
                        tmp[length(tmp) + 1] <- sprintf("(%s)*rx__sens_%s_BY_%s__",
                                                        newLine, rxToSymPy(state), rxToSymPy(eL));

                        ## d^2(h)/d(x,eL) * dx/d(eK)
                        ## d^2(h)/d(x,eL) * sens_state_eK
                        newLine <- rxSymPy(sprintf("diff(diff(%s,%s),%s)",prd, state, rxToSymPy(eL)));
                        tmp[length(tmp) + 1] <- sprintf("(%s)*rx__sens_%s_BY_%s__",
                                                        newLine, rxToSymPy(state), rxToSymPy(eK));

                        ## d^2h/d^2(state) * (dstate/deta_el) * (dstate/deta_eK)
                        ## d^2h/d^2(state) * sens_state_el * sens_state_eK
                        newLine <- rxSymPy(sprintf("diff(diff(%s,%s),%s)",rxToSymPy(prd), rxToSymPy(state), rxToSymPy(state)));
                        tmp[length(tmp) + 1] <- sprintf("(%s)*rx__sens_%s_BY_%s__*rx__sens_%s_BY_%s__",
                                                        newLine, rxToSymPy(state), rxToSymPy(eL), state, rxToSymPy(eK));
                        ## dh/d(state) * dstate/(eK, eL)
                        ## dh/d(state) * sens_state_eK_eL
                        newLine <- rxSymPy(sprintf("diff(%s,%s)",rxToSymPy(prd), rxToSymPy(state)));
                        tmp[length(tmp) + 1] <- sprintf("(%s)*rx__sens_%s_BY_%s_BY_%s__",
                                                        newLine, rxToSymPy(state), rxToSymPy(eK), rxToSymPy(eL));
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
                    rxTest();
                }
            }
        }
        rxProgressStop()
        rxCat("\n##   d(.)^2/(d(theta)d(eta))\n");
        rxProgress(length(calcSens$theta) + length(calcSens$eta));
        for (theta in calcSens$theta){
            for (eta in calcSens$eta){
                ## dh/d(eta,theta)
                newLine2 <- rxSymPy(sprintf("diff(diff(%s,%s),%s)", rxToSymPy(prd), rxToSymPy(eta), rxToSymPy(theta)));
                tmp <- c(newLine2);
                for (state in states){
                    ## d^2h/d(eta,state)*dstate/theta
                    ## d^2h/d(eta,state)*sens_state_theta
                    newLine <- rxSymPy(sprintf("diff(diff(%s,%s),%s)",rxToSymPy(prd), rxToSymPy(eta), rxToSymPy(state)));
                    tmp[length(tmp) + 1] <- sprintf("(%s)*rx__sens_%s_BY_%s__",
                                                    rxToSymPy(newLine), rxToSymPy(state), rxToSymPy(theta));
                    ## d^2h/d(state,theta)*dstate/eta
                    ## d^2h/d(state,theta)*sens_state_eta
                    newLine <- rxSymPy(sprintf("diff(diff(%s,%s),%s)",rxToSymPy(prd), rxToSymPy(state), rxToSymPy(theta)));
                    tmp[length(tmp) + 1] <- sprintf("(%s)*rx__sens_%s_BY_%s__",
                                                    newLine, rxToSymPy(state), rxToSymPy(eta));
                    ## d^2h/d^2(state)*dstate/eta
                    ## d^2h/d^2(state)*sens_state_eta
                    newLine <- rxSymPy(sprintf("diff(diff(%s,%s),%s)",prd, state, state));
                    tmp[length(tmp) + 1] <- sprintf("(%s)*rx__sens_%s_BY_%s__*rx__sens_%s_BY_%s__",
                                                    newLine, rxToSymPy(state), rxToSymPy(theta), rxToSymPy(state), rxToSymPy(eta));
                    ## dh/d(state)*d(h)/(deta,dtheta)
                    newLine <- rxSymPy(sprintf("diff(%s,%s)",prd, state));
                    tmp[length(tmp) + 1] <- sprintf("(%s)*rx__sens_%s_BY_%s_BY_%s__",
                                                    newLine, rxToSymPy(state), rxToSymPy(eta), rxToSymPy(theta));

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
                rxTick();
            }
        }
        rxProgressStop();
        attr(extraLines, "zeroSens") <- zeroSens;
        return(extraLines);
    }
    ## rxCat("## ");

    ## Equation #21
    rxSymPyVars(calcSens);
    rxProgress(length(calcSens));
    for (var in calcSens){
        ##
        newLine2 <- rxSymPy(sprintf("diff(%s,%s)", rxToSymPy(prd), rxToSymPy(var)));
        tmp <- c(newLine2);
        for (state in states){
            newLine <- rxSymPy(sprintf("diff(%s,%s)",rxToSymPy(prd), rxToSymPy(state)));
            tmp[length(tmp) + 1] <- sprintf("(%s)*rx__sens_%s_BY_%s__",
                                            newLine, state, rxToSymPy(var));
        }
        tmp <- paste(paste0("(", tmp, ")"), collapse=" + ")
        if (!pred.minus.dv && prd == "rx_pred_"){
            tmp <- sprintf("-(%s)", tmp);
        }
        tmp <- rxSymPy(tmp);
        if (tmp != "0"){
            zeroSens <- FALSE;
        }
        tmp <- rxSplitLines(sprintf("rx__sens_%s_BY_%s__", prd, rxToSymPy(var)), rxFromSymPy(tmp));
        extraLines[length(extraLines) + 1] <- tmp;
        rxTick();
    }
    rxProgressStop();
    attr(extraLines, "zeroSens") <- zeroSens;
    return(extraLines);
}

genCmt0 <- function(ncmt=1, oral=FALSE){
    ## The lincmt function generates:
    ## 1 cmt: rx_k
    ## 2 cmt: rx_k12, rx_k21
    ## 3 cmt: rx_k13, rx_k31
    rx0 <- "";
    rxc <- "d/dt(central) ~ -rx_k*central";
    rxp <- "";
    rx3 <- ""
    if (ncmt >= 2){
        rxc <- paste(rxc, "- rx_k12*central + rx_k21*peripheral")
        rxp <- "d/dt(peripheral) ~ rx_k12*central - rx_k21*peripheral";
    }
    if (ncmt == 3){
        rxc <- paste(rxc, "- rx_k13*central + rx_k31*peripheral2")
        rx3 <- "d/dt(peripheral2) ~ rx_k13*central - rx_k31*peripheral2";
    }
    if (oral){
        rxc <- paste(rxc, "+ rx_ka*depot");
        rx0 <- "d/dt(depot) ~ -rx_ka*depot";
    }
    fin <- "rx1c = central/rx_v";
    ret <- c(rx0, rxc, rxp, rx3, fin);
    ret <- ret[ret != ""];
    paste(ret, collapse=";\n")
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
            return(sprintf("%s ~ %s;", v, tmp1));
        } else {
            return(NULL);
        }
    }
    oral <- rxSymPyExists("rx_ka");
    if (oral){
        oral <- (rxSymPy("rx_ka") != "0") && (rxSymPy("rx_ka") != "rx_ka");
    }
    dur <- rxSymPyExists("rx_dur");
    if (dur){
        dur <- (rxSymPy("rx_dur") != "0") && (rxSymPy("rx_dur") != "rx_dur")
    }
    rate <- rxSymPyExists("rx_rate");
    if (rate){
        rate <- (rxSymPy("rx_rate") != "0") && (rxSymPy("rx_rate") != "rx_rate")
    }
    tlag <- rxSymPyExists("rx_tlag");
    if (tlag){
        tlag <- (rxSymPy("rx_tlag") != "0") && (rxSymPy("rx_tlag") != "rx_tlag")
    }
    tlag2 <- rxSymPyExists("rx_tlag2");
    if (tlag2){
        tlag2 <- (rxSymPy("rx_tlag2") != "0") && (rxSymPy("rx_tlag2") != "rx_tlag2")
    }
    f2 <- rxSymPyExists("rx_F2");
    if (f2){
        f2 <- suppressWarnings(as.numeric(rxSymPy("rx_F2")));
        if (is.na(f2)) f2  <- 0;
        f2 <- (f2 != 1) && (rxSymPy("rx_F2") != "rx_F2")
    }
    f <- rxSymPyExists("rx_F");
    if (f){
        f <- suppressWarnings(as.numeric(rxSymPy("rx_F")));
        if (is.na(f)) f  <- 0;
        f <- (f != 1) && (rxSymPy("rx_F") != "rx_F")
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
    if (tlag){
        tmp1 <- rxSymPy("rx_tlag");
        tmp1 <- rxFromSymPy(tmp1);
        if (oral){
            extra <- c(extra,
                       sprintf("lag(depot) = %s;", tmp1));
        } else {
            extra <- c(extra,
                       sprintf("lag(central) = %s;", tmp1));
        }
    }
    if (tlag2){
        tmp1 <- rxSymPy("rx_tlag2");
        tmp1 <- rxFromSymPy(tmp1);
        if (oral){
            extra <- c(extra,
                       sprintf("lag(central) = %s;", tmp1));
        } else {
            stop("Cannot handle this tlag combination.");
        }
    }
    if (f){
        tmp1 <- rxSymPy("rx_F");
        tmp1 <- rxFromSymPy(tmp1);
        if (oral){
            extra <- c(extra,
                       sprintf("F(depot) = %s;", tmp1));
        } else {
            extra <- c(extra,
                       sprintf("F(central) = %s;", tmp1));
        }
    }
    if (f2){
        tmp1 <- rxSymPy("rx_F2");
        tmp1 <- rxFromSymPy(tmp1);
        if (oral){
            extra <- c(extra,
                       sprintf("F(central) = %s;", tmp1));
        } else {
            stop("Cannot handle this F(.) combination.");
        }
    }
    if (dur){
        tmp1 <- rxSymPy("rx_dur");
        tmp1 <- rxFromSymPy(tmp1);
        extra <- c(extra,
                   sprintf("dur(central) = %s;", tmp1));
    }
    if (rate){
        tmp1 <- rxSymPy("rx_rate");
        tmp1 <- rxFromSymPy(tmp1);
        extra <- c(extra,
                   sprintf("rate(central) = %s;", tmp1));
    }
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
                       return(sprintf("d/dt(%s) %s %s;", cur.state, sep, v));
                   }),
                   sapply(mv.1$lhs, function(v){
                       v1  <- rxToSymPy(v)
                       v1 <- rxSymPy(v1);
                       v1 <- rxFromSymPy(v1);
                       return(sprintf("%s=%s;", v, v1));
                   })), collapse="\n")
    ret
}

##' Setup Pred function based on RxODE object.
##'
##' This is for the so-called inner problem.
##'
##' @param obj RxODE object
##' @param predfn Prediction function
##' @param pkpars Pk Pars function
##' @param errfn Error function
##' @param init Initialization parameters for scaling.
##' @param grad Boolean indicated if the the equations for the
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
##' @param optExpression Optimize the model text for computer
##'     evaluation.
##' @param run.internal Boolean to see if the function should be run
##'     internally.
##' @param interaction Boolean to determine if \code{dR^2/deta} is
##'     calculated for FOCEi (not needed for FOCE)
##' @return RxODE object expanded with predfn and with calculated
##'     sensitivities.
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
##' @importFrom utils find
rxSymPySetupPred <- function(obj, predfn, pkpars=NULL, errfn=NULL, init=NULL, grad=FALSE, sum.prod=FALSE, pred.minus.dv=TRUE,
                             theta.derivs=FALSE,only.numeric=FALSE,
                             grad.internal=FALSE, theta.internal=FALSE,
                             optExpression=TRUE, run.internal=RxODE.sympy.run.internal,
                             interaction=TRUE){
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
    rxTempDir();
    if (!grad.internal && !theta.internal){
        if (RxODE.cache.directory == "."){
            cache.file <- "";
        } else {
            cache.file <- file.path(RxODE.cache.directory,
                                    sprintf("rx_%s%s.prd",
                                            digest::digest(paste(deparse(list(rxModelVars(obj)$md5["parsed_md5"],
                                                                              ifelse(is.function(predfn), paste(deparse(body(predfn)), collapse=""), ""),
                                                                              ifelse(is.function(pkpars), paste(deparse(body(pkpars)), collapse=""), ""),
                                                                              ifelse(is.function(errfn), paste(deparse(body(errfn)), collapse=""), ""),
                                                                              init, grad, sum.prod, pred.minus.dv, theta.derivs, only.numeric, optExpression,
                                                                              interaction)), collapse="")),
                                            .Platform$dynlib.ext))
        }
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
            save(obj, predfn, pkpars, errfn, init, grad, sum.prod, pred.minus.dv, theta.derivs, only.numeric, optExpression, interaction, file=dta);
            ## on.exit({unlink(dta); unlink(rfile)});
            tmp <- options();
            tmp <- tmp[regexpr(rex::rex(start, "RxODE."), names(tmp)) != -1];
            rf <- sprintf("%s;options(RxODE.delete.unnamed=FALSE);load(%s); require(RxODE);tmp <- rxSymPySetupPred(obj, predfn, pkpars, errfn, init, grad, sum.prod, pred.minus.dv, theta.derivs, only.numeric, optExpression=optExpression, run.internal=TRUE,interaction=interaction);",
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
            assignInMyNamespace("rxErrEnv.lambda", NULL);
            assignInMyNamespace("rxErrEnv.yj", NULL);
            assignInMyNamespace("rxSymPyExpThetas", c());
            assignInMyNamespace("rxSymPyExpEtas", c());
            .mv0 <- rxModelVars(obj);
            .curDvid <- .mv0$dvid;
            .state0  <- .mv0$state;
            if (length(.state0) > 0){
                .state0 <- paste(paste0("cmt(",.mv0$state,");\n"),collapse="");
            } else {
                .state0 <- "";
            }
            .statef  <- .mv0$stateExtra;
            if (length(.statef) > 0){
                .statef <- paste0(paste(paste0("\ncmt(",.mv0$stateExtra,");"),collapse=""),"\n");
            } else {
                .statef <- "";
            }
            if (length(.curDvid) > 1){
                .dvidF <- paste0("\ndvid(",
                                 paste(.curDvid, collapse=","),");\n");
            } else {
                .dvidF <- "";
            }

            .endpoints  <- c(.mv0$state, .mv0$stateExtra)
            if (length(.endpoints) > 0){
                names(.endpoints) <- paste0("(CMT==",seq_along(.endpoints),")")
            }
            gobj <- obj;
            oobj <- genCmtMod(obj);
            obj <- oobj;
            ##rxModelVars(oobj)
            ##rxSymPyVars(obj);
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
                        inis <- paste(paste0(rxState(oobj), "(0)=", inis, "+0.0;"), collapse="\n");
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
                calcSens <- .rxParams(obj, FALSE)
                newmod <- obj;
            }
            ##
            extra.pars <- c();
            if (is.null(errfn)){
                errfn <- function(){add(0.1)};
            }
            ## Get maximum theta.
            txt <- rxParsePred(predfn, init=init);
            pred.mod <- rxGetModel(txt);
            pars <- .rxParams(rxGetModel(paste0(rxNorm(obj), "\n", rxNorm(pred.mod))), FALSE);
            w <- which(sapply(pars, function(x){
                return(substr(x, 0, 2) == "TH")
            }))
            pars <- pars[w];
            mtheta <- max(as.numeric(gsub(rex::rex("THETA[", capture(numbers), "]"), "\\1", pars)))
            .err <- rxParseErr(errfn, base.theta=mtheta + 1, init=init)
            pred.mod <- rxParsePred(predfn, init=init, .err);
            extra.pars <- attr(pred.mod, "ini");
            pred.mod <- rxGetModel(pred.mod)
            full <- paste0(rxNorm(obj), "\n", rxNorm(pred.mod));
            full <- rxGetModel(full);
            etas <- .rxParams(full, FALSE);
            thetas <- etas[regexpr(regTheta, etas) != -1];
            etas <- etas[regexpr(regEta, etas) != -1];
            if (length(etas) > 0 && !theta.internal && !grad.internal && !only.numeric){
                rxCat("Calculate ETA-based prediction and error derivatives:\n")
                calcSens <- etas;
            } else {
                calcSens <- .rxParams(full, FALSE);
            }
            if (grad.internal){
                rxCat("Calculate THETA/ETA-based prediction and error 1st and 2nd order derivatives:\n")
                calcSens <- list(eta=etas);
                calcSens$theta <- thetas;
            }
            if (theta.internal){
                rxCat("Calculate THETA-based prediction and error derivatives:\n")
                calcSens <- thetas;
            }
            .baseState <- rxState(obj);
            .inits0 <- rxInits(obj);
            .inits <- .inits0
            .oLhs <- c(names(.inits),rxLhs(obj));
            if (length(.inits) > 0){
                .inits <- rxInits(obj, rxLines=TRUE);
            } else {
                .inits <- "";
            }

            .full <- rxGetModel(paste0(rxNorm(obj), "\n", rxNorm(pred.mod)));
            ## Now everything is setup get conditional statements and move on.
            .cond <- rxExpandIfElse(.full);
            .ncond <- names(.cond);
            .zeroSens <- FALSE;
            if (only.numeric) calcSens = FALSE;
            .mods <- lapply(seq_along(.cond), function(.i){
                if ((regexpr(rex::rex("rx_pred_="), .cond[.i]) == -1) |
                    (regexpr(rex::rex("rx_r_="), .cond[.i]) == -1)){
                    return(NULL)
                } else {
                    if (!is.null(.ncond)){
                        rxCat("################################################################################\n");
                        if (any(names(.endpoints)==.ncond[.i])){
                            rxCat(sprintf("## %s %s\n", crayon::bold("Endpoint:"), .endpoints[.ncond[.i]]));
                        } else {
                            rxCat(sprintf("## %s %s\n", crayon::bold("Condition:"), .ncond[.i]));
                        }
                        rxCat("################################################################################\n");
                    }
                    .origStates <- rxState(.cond[.i])
                    if (length(.origStates) > 0){
                        .cond <- gsub(rex::rex(capture(or("rx_lambda_", "rx_yj_")), "~"), "\\1=", .cond[.i]);
                        .full <- rxGetModel(.cond, calcSens=calcSens, collapseModel=TRUE);
                        .full <- rxGetModel(paste0(.inits, gsub(rex::rex(capture(or("rx_lambda_", "rx_yj_")), "="), "\\1~", rxNorm(.full))));
                    } else {
                        .full <- rxGetModel(paste0(.inits, .cond[.i]));
                    }
                    .fullState <- rxState(.full);
                    rxCat("Load into sympy...");
                    rxSymPySetup(.full);
                    if (!only.numeric && !is.null(.ncond)){
                        .tmp <- .ncond[.i];
                        if (regexpr(rex::rex("(!"), .tmp) != -1){
                            .tmp <- sub(rex::rex(start, any_spaces, "(!", capture(anything), ")", any_spaces, end),
                                        "\\1", .tmp);
                        }
                        .p <- rxModelVars(paste0("if ", .tmp, "{rx_tmp=1.0+2.0;}"))$params;
                        for (.v in .p){
                            .tmp <- rxToSymPy(.v);
                            if (rxSymPyExists(.tmp)){
                                .tmp <- rxSymPy(.tmp);
                                .tmp <- rxFromSymPy(.tmp);
                                if (regexpr(rex::rex(except_some_of("H"), "ETA[", any_numbers, "]"), .tmp) != -1){
                                    stop("if/else depends on ETAs, currently not supported by FOCEi/CWRES");
                                }
                            }
                        }
                    }
                    if (!is.null(.ncond)){
                        ## FIXME use parsing to fix these LHS quantities.
                        for(.v in .oLhs[order(sapply(.oLhs, function(.x){ -nchar(.x)}))]){
                            .tmp <- rxToSymPy(.v);
                            if (rxSymPyExists(.tmp)){
                                .tmp <- rxSymPy(.tmp);
                                .tmp <- rxFromSymPy(.tmp);
                                .ncond[.i] <<- gsub(rex::rex(.v), .tmp, .ncond[.i]);
                            }
                        }
                    }
                    on.exit({rxCat("Freeing Python/SymPy memory...");rxSymPyClean();rxCat("done\n")});
                    rxCat("done\n");
                    if (rxSymPyExists("rx_pred_") & rxSymPyExists("rx_r_")){
                        if (!only.numeric){
                            if (.useUtf()){
                                rxCat("Calculate \u2202(f)/\u2202(\u03B7)\n");
                            } else {
                                rxCat("Calculate d(f)/d(eta)\n");
                            }
                            .newlines <- rxSymPySetupDPred(.full, calcSens, .baseState);
                            .zeroSens <<- .zeroSens | attr(.newlines, "zeroSens")
                            if (interaction){
                                if (.useUtf()){
                                    rxCat("Calculate \u2202(R\u00B2)/\u2202(\u03B7)\n");
                                } else {
                                    rxCat("Calculate d(R^2)/d(eta)\n");
                                }
                                .r0 <- rxSymPy("rx_r_");
                                .newlinesR <- rxSymPySetupDPred(.full, calcSens, .baseState, prd="rx_r_");
                            } else {
                                ## Calculate r(ETAs=0) for FOCE
                                .sub <- sprintf("rx_r_.subs([%s])",
                                        paste(sapply(calcSens, function(x){
                                            sprintf("(%s, 0)", rxToSymPy(x));
                                        }), collapse=", "))
                                .r0 <- rxSymPy(.sub);
                                .newlinesR <- "";
                            }
                        }
                        .states <- paste(sapply(.fullState, function(x){
                            .ini <- sprintf("%s(0)", x);
                            .iniS <- rxToSymPy(.ini);
                            ## Make sure to get any initial conditions
                            if (rxSymPyExists(.iniS)){
                                .iniS <- rxSymPy(.iniS);
                                .iniS <- rxFromSymPy(.iniS);
                                if (.iniS != "0"){
                                    .iniS <- sprintf("%s=%s;\n", .ini, .iniS);
                                } else {
                                    .iniS <- "";
                                }
                            } else {
                                .iniS <- ""
                            }
                            .ddt <- sprintf("d/dt(%s)", x);
                            .sympy <- rxToSymPy(.ddt);
                            .v <- rxSymPy(.sympy);
                            .v <- rxFromSymPy(.v);
                            ## Now add f/dur/rate/alag
                            .f <- sprintf("f(%s)", x);
                            .fS <- rxToSymPy(.f);
                            if (rxSymPyExists(.fS)){
                                .fS <- rxSymPy(.fS);
                                .fS <- rxFromSymPy(.fS);
                                .fS <- sprintf("%s=%s;\n", .f, .fS);
                            } else {
                                .fS <- "";
                            }
                            .dur <- sprintf("dur(%s)", x);
                            .durS <- rxToSymPy(.dur);
                            if (rxSymPyExists(.durS)){
                                .durS <- rxSymPy(.durS);
                                .durS <- rxFromSymPy(.durS);
                                .durS <- sprintf("%s=%s;\n", .dur, .durS);
                            } else {
                                .durS <- "";
                            }
                            .rate <- sprintf("rate(%s)", x);
                            .rateS <- rxToSymPy(.rate);
                            if (rxSymPyExists(.rateS)){
                                .rateS <- rxSymPy(.rateS);
                                .rateS <- rxFromSymPy(.rateS);
                                .rateS <- sprintf("%s=%s;\n", .rate, .rateS);
                            } else {
                                .rateS <- "";
                            }
                            .lag <- sprintf("alag(%s)", x);
                            .lagS <- rxToSymPy(.lag);
                            if (rxSymPyExists(.lagS)){
                                .lagS <- rxSymPy(.lagS);
                                .lagS <- rxFromSymPy(.lagS);
                                .lagS <- sprintf("%s=%s;\n", .lag, .lagS);
                            } else {
                                .lagS <- "";
                            }
                            if (any(x == .origStates)){
                                return(paste0(.iniS,.ddt,"=",.v,";",
                                              .fS,.durS, .rateS, .lagS));
                            } else {
                                return(paste0(.iniS,.ddt,"~",.v,";",
                                              .fS,.durS, .rateS, .lagS));
                            }
                        }), collapse="\n");
                        .yj <- rxToSymPy("rx_yj_")
                        .yj <- rxSymPy(.yj);
                        .lambda <- rxToSymPy("rx_lambda_")
                        .lambda <- rxSymPy(.lambda);
                        .states <- paste0(.states,
                                          "\nrx_yj_~", rxFromSymPy(.yj), ";\n",
                                          "rx_lambda_~", rxFromSymPy(.lambda), ";\n");

                        .pred <- rxToSymPy("rx_pred_");
                        .pred <- rxSymPy(.pred);
                        .sensState <- .fullState[!(.fullState %in% .origStates)]
                        if (length(.sensState) > 1){
                            .reg <- rex::rex("d/dt(", or(.sensState), ")", any_of("=~"), except_any_of("\n;"), any_of("\n;"));
                            .reg2 <- rex::rex(or(.sensState), "(0)", any_of("="), except_any_of("\n;"), any_of("\n;"));
                            .pred.only <- paste0(gsub(.reg2, "", gsub(.reg, "", .states)),
                                                 "rx_pred_=", rxFromSymPy(.pred), ";");
                        } else {
                            .pred.only <- paste0(.states, "rx_pred_=", rxFromSymPy(.pred), ";");
                        }
                        .r <- rxToSymPy("rx_r_");
                        .r <- rxSymPy(.r);
                        .mtime <- sapply(.rxMtimes, function(x){
                            .mtime <- rxToSymPy(x);
                            if (rxSymPyExists(.mtime)){
                                .mtime <- try(rxSymPy(x), silent=TRUE);
                                if (inherits(.mtime, "try-error")) return("");
                                mtime <- try(rxFromSymPy(x), silent=TRUE)
                                if (inherits(.mtime, "try-error")) return("");
                                return(sprintf("%s=%s;", x, .mtime));
                            }
                            return("");
                        })
                        .mtime <- .mtime[.mtime != ""];

                        if (!only.numeric){
                            .inner <- paste0(.states,
                                             "rx_pred_=", rxFromSymPy(.pred), ";")
                            .inner <- paste(c(.inner, .newlines, .mtime,
                                            paste0("rx_r_=", rxFromSymPy(.r0), ";"),
                                            .newlinesR), collapse="\n");
                        } else {
                            .inner <- NULL
                        }
                        ## PRED only R has etas on it.
                        .pred.only <- paste0(.pred.only, paste(.mtime, collapse="\n"),
                                             paste0("rx_r_=", rxFromSymPy(.r), ";"));
                        .lhs <- setNames(sapply(.oLhs, function(x){
                            .lhs <- rxToSymPy(x);
                            if (rxSymPyExists(.lhs)){
                                .lhs <- try(rxSymPy(x), silent=TRUE);
                                if (inherits(.lhs, "try-error"))  return("");
                                if (x == .lhs){
                                    return("");
                                } else {
                                    .lhs <- try(rxFromSymPy(.lhs), silent=TRUE);
                                    if (inherits(.lhs, "try-error"))  return("");
                                    return(sprintf("%s=%s;", x, .lhs));
                                }
                            }
                            return("");
                         }), .oLhs);
                        .lhs <- .lhs[which(.lhs != "")];
                        .lhs <- gsub(rex::rex(capture("nlmixr_",except_any_of("=")),"="),"\\1~",.lhs)
                        .pred.only  <- paste0(.pred.only, paste(.lhs, collapse="\n"));
                        .toLines <- function(x){
                            if(is.null(x)) return(NULL);
                            strsplit(rxNorm(rxGetModel(x)), "\n")[[1]]
                        }
                        return(list(pred.only=.toLines(.pred.only),
                                    lhs=.lhs,
                                    inner=.toLines(.inner)));
                    } else {
                        rxCat("Does not have predictions or errors, skipping.\n");
                        return(NULL)
                    }
                }
            })
            w <- which(sapply(seq_along(.mods), function(x){return(is.null(.mods[[x]]))}));
            if (length(w) > 0){
                .mods <- .mods[-w];
                .cond <- .cond[-w];
                .ncond <- .ncond[-w];
            }
            if (length(.mods) == 1){
                if (.zeroSens){
                    stop("Prediction doesn't depend on the ETAs")
                }
                .mods <- .mods[[1]]
                .pred.only <- paste(.mods$pred.only, collapse="\n");
                if (!only.numeric){
                    .inner <- paste(.mods$inner, collapse="\n");
                } else {
                    .inner <- NULL;
                }
            } else {
                .ord <- order(sapply(.ncond, nchar));
                .mods <- .mods[.ord];
                .cond <- .cond[.ord];
                .ncond <- .ncond[.ord];
                .collapseModel <- function(what){
                    .lineNums <- seq_along(.mods[[1]][[what]]);
                    .len <- length(.mods);
                    .linesEqual <- function(x){
                        .lines <- sapply(seq(1, .len), function(y){
                            .mods[[y]][[what]][x]
                        })
                        all(.lines[1] == .lines[-1])
                    }
                    .eq <- sapply(.lineNums, .linesEqual);
                    .lines <- sapply(seq_along(.eq), function(i){
                        if (.eq[i]) return(.mods[[1]][[what]][i]);
                        .ret <- character(.len);
                        .oi <- i;
                        while(i <= length(.eq) & !.eq[i]){
                            for(y in seq(1, .len)){
                                .ret[y] <- paste(c(.ret[y], .mods[[y]][[what]][i]), collapse="\n");
                                .mods[[y]][[what]][i] <<- "";
                                .eq[i] <<- TRUE;
                            }
                            i <- i + 1;
                        }
                        .ret <- paste(sapply(seq(1, .len), function(y){
                            if (y > 1){
                                if (.ncond[y] == sprintf("(!%s)", .ncond[y - 1])){
                                    return(sprintf("else {%s\n}", .ret[y]))
                                }
                            }
                            ## FIXME: better branching for lhs conditions
                            ## Currently none.
                            .extra <- "";
                            return(sprintf("%sif %s {%s\n}", .extra, .ncond[y], .ret[y]));
                        }), collapse="\n");
                    })
                    .lines <- .lines[.lines != ""];
                    return(paste(.lines, collapse="\n"));
                }
                .pred.only <- .collapseModel("pred.only");
                .inner <- .collapseModel("inner");
            }
            rxCat(paste0(crayon::bold("################################################################################"), "\n"));
            toRx <- function(x, what){
                if (is.null(x)) return(NULL);
                if (x == "") return(NULL);
                .miss <- rxModelVars(x)$params
                .missEtas <- .miss[regexpr(regEta, .miss) != -1];
                if (length(setdiff(.missEtas,etas)) !=0 && length(setdiff(etas, .missEtas))!=0){
                    stop("Predictions do not depend on all the etas defined")
                }
                .miss <- rxModelVars(x)$params
                .missEtas <- .miss[regexpr(regEta, .miss) != -1];
                if (length(setdiff(.missEtas,etas)) !=0 && length(setdiff(etas, .missEtas))!=0){
                    print(setdiff(.missEtas,etas));
                    print(setdiff(etas, .missEtas));
                    stop("Predictions do not depend on all the individual deviations (etas) defined")
                }
                .missTheta <- .miss[regexpr(regTheta,.miss) != -1];
                if (length(setdiff(.missTheta,thetas)) !=0 && length(setdiff(thetas, .missTheta))!=0){
                    print(setdiff(.missTheta,thetas));
                    print(setdiff(thetas, .missTheta));
                    stop("Predictions do not depend on all the population parameters (theas) defined")
                }
                if (sum.prod){
                    rxCat(sprintf("Stabilizing round off errors in products & sums in %s model...", what));
                    x <- rxSumProdModel(x);
                    rxCat("done\n");
                }
                if (optExpression){
                    rxCat(sprintf("Optimizing expressions in %s model...", what));
                    .mod <- paste0(.state0, rxOptExpr(x), .statef, .dvidF)
                    rxCat("done\n");
                } else {
                    .mod <- paste0(.state0, x, .statef, .dvidF);
                }
                rxCat(sprintf("Compiling %s model...", what));
                .ret <- RxODE(.mod);
                rxCat("done\n");
                return(.ret);
            }
           .extraProps <- NULL
            if (!is.null(.inner)){
                ## check for power expressions for appropriate scaling.
                if (.inner == ""){
                    .inner <- NULL
                } else {
                    .extraProps <- .rxFindPow(.inner)
                }
            }
            ret <- list(obj=oobj,
                        pred.only=toRx(.pred.only, "Predictions/EBE"),
                        inner=toRx(.inner, "Inner"),
                        extra.pars=extra.pars,
                        outer=NULL,
                        theta=NULL,
                        warn=.zeroSens,
                        pred.minus.dv=pred.minus.dv,
                        log.thetas=rxSymPyExpThetas,
                        log.etas=rxSymPyExpEtas,
                        extraProps=.extraProps,
                        cache.file=cache.file)
            class(ret) <- "rxFocei";
            save(ret, file=cache.file);
            if (only.numeric){
                rxCat("Standardized prediction/ebe models produced.\n");
            } else {
                rxCat(sprintf("The model-based sensitivities have been calculated%s.\n", ifelse(grad, " (with FOCEi Global Gradient)", "")));
            }
            if (ret$warn){
                warning("Some of your prediction function does not depend on the state variables.");
            }
            return(ret)
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
##' First Order Expansion of ETA
##'
##' @param expr RxODE model
##' @return Return a RxODE model with first order Taylor expansion around ETA
##' @export
##' @author Matthew L. Fidler
rxFoExpandEta <-function(expr){
    rxSymPySetup(expr);
    .vars <- .rxParams(expr)
    .etas <- .vars[regexpr(rex::rex(start, "ETA[", any_numbers, "]", end), .vars) != -1]
    .fn <- function(line, .w="~"){
        .l2 <- strsplit(line, .w)[[1]]
        if (length(.l2) == 2){
            .l2.1 <- .l2[1];
            .lhs <- rxFromSymPy(.l2.1)
            for (.e in .etas){
                .e2 <- rxToSymPy(.e)
                .tmp <- sprintf("%s = %s.subs(%s, 0) + diff(%s, %s).subs(%s,0)*%s", .lhs, .lhs, .e2, .lhs, .e2, .e2, .e2);
                rxSymPy(.tmp)
            }
            .l3 <- rxSymPy(.lhs)
            return(sprintf("%s%s%s", .l2.1, .w, rxFromSymPy(.l3)))
        } else {
            stop("Malformed RxODE block.")
        }
    }
    .lines <- sapply(strsplit(rxNorm(expr), "\n+")[[1]],
                     function(line){
        if (regexpr("~", line) != -1){
            return(.fn(line, "~"))
        } else {
            return(.fn(line, "="))
        }
    })
    return(paste(.lines, collapse="\n"));
}
##' This creates the dv/dx RxODE model for a linear solved system
##'
##' @param model RxODE pred model
##' @param ncmt Number of compartments of the solved system
##' @param parameterization Integer representing parameterization type
##' @param optExpression boolean indicating if you should optimize the
##'     expression.
##' @return RxODE model expressing \code{dvdx}
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxSymPyLincmtDvdx <- function(model, ncmt, parameterization, optExpression=TRUE){
    .mod <- rxGetModel(model)
    rxTempDir();
    if (RxODE.cache.directory == "."){
        .file <- "";
    } else {
        .file <- file.path(RxODE.cache.directory,
                           sprintf("rx_%s.dvdx",
                                   digest::digest(paste(deparse(c(rxModelVars(.mod)$md5["parsed_md5"], ncmt, parameterization, optExpression)),
                                                        collapse=""))));
    }
    if (file.exists(.file)){
        readRDS(.file);
    } else {
        .pm <- list(
            c("CL", "V", "KA", "TLAG"),
            c("CL", "V", "CLD", "VT", "KA", "TLAG"),
            c("CL", "V", "CLD", "VT", "CLD2", "VT2", "KA", "TLAG"),
            c("KE", "V", "KA", "TLAG"),
            c("KE", "V", "K12", "K21", "KA", "TLAG"),
            c("KE", "V", "K12", "K21", "K13", "K31", "KA", "TLAG")
        )
        dim(.pm)<-c(3, 2)
        .pars <- .pm[ncmt, parameterization][[1]];

        .lhs <- rxLhs(.mod)
        .etas <- .rxParams(.mod);
        .etas <- .etas[regexpr(rex::rex(start, "ETA[",any_numbers, "]"), .etas) != -1]
        .etas <- max(as.numeric(sub(rex::rex(start, "ETA[",capture(any_numbers), "]"), "\\1", .etas)))
        rxSymPySetup(model)
        .len <- .etas * length(.pars);
        .mat <- character(.len)
        .i <- 1;
        message("Calculate d(pred)/d(eta)")
        rxProgress(.etas * length(.pars))
        for (.eta in sprintf("ETA[%d]", seq(1, .etas))){
            for (.p in .pars){
                if (any(.lhs == .p)){
                    .line <- sprintf("diff(%s,%s)", rxToSymPy(.p), rxToSymPy(.eta));
                    .line <- rxSymPy(.line);
                    .line <- rxFromSymPy(.line)
                    .mat[.i] <- .line
                } else {
                    .mat[.i] <- "0"
                }
                .i <- .i + 1
                rxTick()
            }
        }
        rxProgressStop()
        .mod <- paste(sprintf("rxMat%d=%s", seq_along(.mat) - 1, .mat), collapse="\n")
        dim(.mat) <- c(length(.pars), .etas)
        dimnames(.mat) <- list(.pars, sprintf("ETA[%d]", seq(1, .etas)))
        print(.mat)
        if (optExpression){
            rxCat(sprintf("Optimizing expressions in pred..."));
            .mod <- rxOptExpr(.mod)
            rxCat("done\n");
        }
        saveRDS(.mod, .file)
        .mod
    }
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
