.rxConvertCmt <- function(cmt, state){
    ## This converts CMT to compartment#
    if (is(cmt, "numeric")){
        return(cmt);
    } else if (is(cmt, "factor")){
        .uniqueCmt <- levels(cmt);
    } else if (is(cmt, "character")){
        .uniqueCmt <- unique(cmt);
    }
    .uniqueCmt <- .uniqueCmt[.uniqueCmt != "(default)"];
    .uniqueCmt <- .uniqueCmt[.uniqueCmt != "(obs)"];
    if (length(.uniqueCmt)==0L){
        .lvlCmt <- c("(default)","(obs)");
    } else {
        .lvlCmt <- setdiff(.uniqueCmt, state);
        .lvlCmt <- c(state, .lvlCmt, "(obs)");
    }
    .cmt <- paste0(cmt);
    .cmt <- ifelse(.cmt=="(default)",.lvlCmt[1],.cmt);
    .cmt <- as.numeric(factor(.cmt,levels=.lvlCmt));
    .cmt <- ifelse(.cmt==length(.lvlCmt),0,.cmt);
    .lvlCmt <- .lvlCmt[-length(.lvlCmt)];
    attr(.cmt,"cmt") <- .lvlCmt;
    return(.cmt);
}

##' Translate Data into solving ready data.
##'
##' @param data Data to convert
##' @param model Model to use for conversion
##' @param covs Covariate suport for backward compatibility
##' @return Converted data frame
##' @author Matthew Fidler
rxEvTrans <- function(data, model, covs){
    .d <- data
    if (!is.null(covs)){
        if (!rxIs(.d, "rxEt")){
            stop("This only works with RxODE event tables; integrate covariates into a data.frame.")
        }
        if (.d$maxId != 1){
            stop("Covariate interpolation only works with single subject data; Add the covariate information to the data.frame")
        }
        .cov <- as.matrix(covs);
        .covLen <- dim(.cov)[1];
        if (.covLen !=  length(.d$time)){
            .samplingTime <- .d$get.sampling()$time;
            if (.covLen != length(.samplingTime))
                stop("Covariate length need to match the sampling times or all the times in the event table.");
            .lst <- as.matrix(do.call("cbind", lapply(seq(1L, dim(.cov)[2]), function(i){
                                                   f <- stats::approxfun(.samplingTime, .cov[, i])
                                                   return(f(.d$time))
                                               })))
            dimnames(.lst) <- list(NULL, dimnames(.cov)[[2]]);
            .cov <- .lst;
        }
        .d <- cbind(as.data.frame(.d), as.data.frame(.cov));
    } else {
        .d <- as.data.frame(data)
    }
    .colNames <- colnames(.d)
    .colNames <- tolower(.colNames)
    .anyCmt <- which(.colNames == "cmt")
    .state <- rxModelVars(model);
    if (length(.anyCmt)==1){
        .d[,.anyCmt] <- .rxConvertCmt(.d[,.anyCmt], .state);
    }
    evTrans(.d,model)
}
