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
##' @return Converted data frame
##' @author Matthew Fidler
rxEvTrans <- function(data, model){
    .d <- as.data.frame(data)
    .colNames <- colnames(.d)
    .colNames <- tolower(.colNames)
    .anyCmt <- which(.colNames == "cmt")
    .state <- rxModelVars(model);
    if (length(.anyCmt)==1){
        .d[,.anyCmt] <- .rxConvertCmt(.d[,.anyCmt], .state);
    }
    evTrans(.d,model)
}
