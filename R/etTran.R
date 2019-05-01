.convertDvid  <- function(id, maxDvid=0L){
    .udvid  <- sort(unique(id))
    if (max(.udvid) > maxDvid){
        .ndvid  <- seq_along(.udvid);
        as.integer(factor(id,levels=.udvid,.ndvid))
    } else {
        return(id)
    }
}
.convertId <- function(id){
    .pid <- paste(id);
    .lvl <- unique(.pid);
    .lab <- paste(.lvl)
    factor(.pid, levels=.lvl, labels=.lab);
}

.convertExtra <- function(dat){
    d <- as.data.frame(dat);
    .colNames0 <- colnames(d)
    col.names <- toupper(.colNames0);
    ## Handle DATE TIME; DAT1 TIME; DAT2 TIME and DAT3 TIME
    ## Note NONMEM handles dates of the format DAY-MONTH and DAY as well for the DATE class of objects.
    ## It is too complex to handle, and not very common so it will throw an error
    do.date <- FALSE;
    bad.msg <- "Dates formatted as MONTH-DAY or DAY alone are not supported in this conversion"
    dup.date <- "Dates can only be specified by one of DATE, DAT1, DAT2, DAT3 / TIME;  This data has multiple DATE columns"
    check.bad <- function(d){
        d <- paste(d);
        if (any(unlist(lapply(strsplit(d, "[^0-9]+"), length)) != 3)){
            stop(bad.msg)
        }
        return(d)
    }
    if (any(col.names == "DATE")){
        ##  Month Day Year
        dat.reg2 <- rex::rex(start, any_spaces, capture(numbers), non_numbers, capture(numbers), non_numbers, capture(number, number), any_spaces, end)
        dat.reg4 <- rex::rex(start, any_spaces, capture(numbers), non_numbers, capture(numbers), non_numbers, capture(number, number, number, number), any_spaces, end)
        dt <- check.bad(d$DATE)
        d$DATE.TIME <- as.POSIXct(NA);
        w <- which(regexpr(dat.reg2, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(dat.reg2, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format="%m-%d-%y %H:%M")
        w <- which(regexpr(dat.reg4, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(dat.reg4, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format="%m-%d-%Y %H:%M")
        d <- d[, -which(names(d) == "DATE")];
        do.date <- TRUE;
    }
    if (any(col.names == "DAT1")){
        if (do.date){
            stop(dup.date)
        }
        ## DAT1   day month year
        dat.reg2 <- rex::rex(start, any_spaces, capture(numbers), non_numbers, capture(numbers), non_numbers, capture(number, number), any_spaces, end)
        dat.reg4 <- rex::rex(start, any_spaces, capture(numbers), non_numbers, capture(numbers), non_numbers, capture(number, number, number, number), any_spaces, end)
        dt <- check.bad(d$DAT1)
        d$DATE.TIME <- as.POSIXct(NA);
        w <- which(regexpr(dat.reg2, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(dat.reg2, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format="%d-%m-%y %H:%M")
        w <- which(regexpr(dat.reg4, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(dat.reg4, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format="%d-%m-%Y %H:%M")
        d <- d[, -which(names(d) == "DAT1")];
        do.date <- TRUE;
    }
    if (any(col.names == "DAT2")){
        ## DAT2   year month day
        if (do.date){
            stop(dup.date)
        }
        dat.reg2 <- rex::rex(start, any_spaces, capture(number, number), non_numbers, capture(numbers), non_numbers, capture(numbers), any_spaces, end)
        dat.reg4 <- rex::rex(start, any_spaces, capture(number, number, number, number), non_numbers, capture(numbers), non_numbers, capture(numbers), any_spaces, end)
        dt <- check.bad(d$DAT2)
        d$DATE.TIME <- as.POSIXct(NA);
        w <- which(regexpr(dat.reg2, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(dat.reg2, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format="%y-%m-%d %H:%M")
        w <- which(regexpr(dat.reg4, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(dat.reg4, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format="%Y-%m-%d %H:%M")
        d <- d[, -which(names(d) == "DAT2")];
        do.date <- TRUE;
    }
    if (any(col.names == "DAT3")){
        ## DAT3   year day month
        if (do.date){
            stop(dup.date)
        }
        dat.reg2 <- rex::rex(start, any_spaces, capture(number, number), non_numbers, capture(numbers), non_numbers, capture(numbers), any_spaces, end)
        dat.reg4 <- rex::rex(start, any_spaces, capture(number, number, number, number), non_numbers, capture(numbers), non_numbers, capture(numbers), any_spaces, end)
        dt <- check.bad(d$DAT3)
        d$DATE.TIME <- as.POSIXct(NA);
        w <- which(regexpr(dat.reg2, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(dat.reg2, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format="%y-%d-%m %H:%M")
        w <- which(regexpr(dat.reg4, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(dat.reg4, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format="%Y-%d-%m %H:%M")
        d <- d[, -which(names(d) == "DAT3")];
        do.date <- TRUE;
    }
    if (do.date){
        if (any(is.na(d$DATE.TIME))){
            stop("The date time format was not correctly specified")
        }
    }
    if (do.date){
        ## Sort by date/time (though this should have been done already...)
        if (!any(names(d)=="ID")){
            d$ID  <- 1L;
        }
        if (!any(names(d)=="EVID")){
            d$EVID  <- 0L;
        }
        d <- d[order(d$ID, d$DATE.TIME, -d$EVID), ];
        d$TIME <- as.vector(unlist(sapply(unique(d$ID), function(id){
            d0 <- d[d$ID == id, ];
            return(as.numeric(difftime(d0$DATE.TIME, d0$DATE.TIME[1], units="hours")))
        })))
        d <- d[, -which(names(d) == "DATE.TIME")];
    }
    if (is(d$TIME, "numeric") || is(d$TIME, "integer")) return(d)
    stop("Cannot figure out a numeric time")
}

##'@export
print.rxEtTran <- function(x,...){
    print(as.data.frame(x));
    .cls <- class(x);
    .lst <- attr(.cls, ".RxODE.lst")
    cat("\nCovariates (non time-varying):\n")
    print(.lst$cov1)
    cat("\nCompartment translation:\n");
    print(data.frame("Compartment Name"=.lst$cmtInfo,"Compartment Number"=seq_along(.lst$cmtInfo),
                     check.names=FALSE))
}
