.convertDvid  <- function(id, maxDvid=0L) {
    .udvid  <- sort(unique(id))
    if (max(.udvid) > maxDvid) {
        .ndvid  <- seq_along(.udvid);
        as.integer(factor(id, levels = .udvid, .ndvid))
    } else {
        return(id)
    }
}

##' Converts an ID column into a factor
##'
##' @param id Id column
##' @return A factor of IDs
##' @author Matthew Fidler
##' @noRd
.convertId <- function(id) {
  .pid <- paste(id);
  .lvl <- unique(.pid);
  .lab <- paste(.lvl)
  factor(.pid, levels = .lvl, labels = .lab);
}

##' Get the nesting info for a single column
##'
##' @param id Id information
##' @param col Nesting
##' @return A factor made with `.convertId()`; When the nesting is
##'   above the `id` level supplied (like site variability) nothing
##'   else is added; However, if the nesting is above the `id` level
##'   (like occasion) it has an extra attribute `nu` that shows the
##'   total information/degrees of freedom for the combination of id
##'   and the lower nesting level.
##' @author Matthew Fidler
##' @noRd
.nestingInfoSingle <- function(col, id) {
  .f1 <- factor(paste(id, col))
  .l1 <- length(levels(.f1))
  .f2 <- .convertId(col)
  .lid <- length(levels(id))
  if (.l1 == .lid) {
    ## Case:
    ##  study id
    ##  1     1
    ##  1     2
    ##  1     3
    ##  2     4
    ##  2     5
    ##
    ## The factor(paste(study,id)) will have the same number of levels
    ## as factor(paste(id))
    return(.f2)
  } else if (.l1 > .lid) {
    ## Case:
    ##  id  occ
    ##  1     1
    ##  1     2
    ##  1     3
    ##  2     1
    ##  2     2
    ##
    ## The factor(paste(occ,id)) will have more levels than
    ## factor(paste(id))
    attr(.f2, "nu") <- .l1
    return(.f2)
  } else {
    stop("un-handled nesting information")
  }
}

##' This gets the nesting information
##'
##' @param id Id information
##' @param omega Omega information
##' @param data
##' @return Nesting Information
##' @author Matthew Fidler
##' @noRd
.nestingInfo <- function(id, omega, data) {
  .env <- new.env(parent=emptyenv())
  if (!inherits(id, "factor")) {
    id <- .convertId(id)
  }
  .lowNames <- tolower(names(data))
  .w <- which(.lowNames == "id")
  if (length(.w) == 1) {
      .id <- names(data)[.w]
  } else {
      stop("unsupported data when calculating nesting information")
  }
  omega <- lotri::as.lotri(omega, default=.id)
  .lvls <- names(omega)
  .lvls <- .lvls[.lvls != .id]
  env <- new.env(parent=emptyenv())
  env$data <- data
  env$above <- c()
  env$aboveVars <- list()
  env$below <- c()
  env$belowVars <- list()
  env$extraTheta <- 0
  env$extraEta <- 0
  lapply(.lvls, function(lvl) {
      .s <- .nestingInfoSingle(data[[lvl]], id)
      .l1 <- length(levels(.s))
      .dn <- dimnames(omega[[lvl]])[[1]]
      if (is.null(attr(.s, "nu"))) {
          env$above <- c(env$above, setNames(.l1, lvl))
          env$aboveVars <- c(env$aboveVars,
                             setNames(list(.dn), lvl))
          env$extraTheta <- env$extraTheta + length(.dn) * .l1
      } else {
          env$below <- c(env$below, setNames(.l1, lvl))
          env$belowVars <- c(env$belowVars,
                             setNames(list(.dn), lvl))
          env$extraEta <- env$extraEta + length(.dn) * .l1
      }
      env$data[[lvl]] <- .s;
  })
  return(list(data=env$data,
              omega=omega,
              idName=.id,
              id=id,
              above=env$above,
              below=env$below,
              aboveVars=env$aboveVars,
              belowVars=env$belowVars,
              extraTheta=env$extraTheta,
              extraEta=env$extraEta))
}

.warnIdSort0 <- TRUE
##' Turn on/off warnings for ID sorting.
##'
##' @param warnIdSort Boolean for if the sorting warning is turned on
##'     or off.
##' @return Nothing
##' @author Matthew Fidler
##' @export
.setWarnIdSort <- function(warnIdSort=TRUE) {
    assignInMyNamespace(".warnIdSort0", warnIdSort);
    invisible()
}

.DTEnv <- NULL
.getDTEnv <- function() {
  if (is.null(.DTEnv)) {
    if (requireNamespace("data.table", quietly = TRUE)) {
      .env <- loadNamespace("data.table");
      if (utils::compareVersion(as.character(
        utils::packageVersion("data.table")),
        "1.12.4") >= 0) {
        assignInMyNamespace(".DTEnv", .env)
        return(.env)
      }
    }
    .env <- new.env(parent=emptyenv())
    assignInMyNamespace(".DTEnv", .env)
    return(.env)
  } else {
    return(.DTEnv)
  }
}

##' This function sorts the parameter or iCov data based on the event
##' table data
##'
##' @param idData This is the individual parameter data or iCov data
##' @param goodLvl These are the "good" levels based on the event
##'     table
##' @param type Type can be "iCov" or "parameter"
##' @param warnIdSort When `TRUE` warnings about merging the
##'     parameter/id with RxODE event tables are issued.
##' @return A sorted parameter table that can be used directly in the
##'     C-based routines.
##' @author Matthew Fidler
##' @noRd
.sortId <- function(idData, goodLvl, type="parameter",
                    warnIdSort, skipStop=TRUE) {
    .n <- tolower(names(idData));
    .w <- which(.n == "id");
    .nid <- length(goodLvl);
    if (length(.w) == 1) {
        if (.nid == 1 && length(idData[[.w]]) == 1) {
            .idData <- as.data.frame(idData);
            goodLvl <- paste(.idData[[.w]]);
        } else {
            .idData <- as.data.frame(idData);
        }
        .oId <- .idData[[.w]];
        .idData[[.w]] <- factor(.idData[[.w]], levels = goodLvl, labels = goodLvl);
        .wrn <- ""
        if (any(is.na(.idData[[.w]]))) {
            .w2 <- which(is.na(.idData[[.w]]));
            .oId <- unique(.oId[.w2])
            .wrn <- sprintf("Some IDs are in the %s dataset that are not in the event dataset.\nParameter information for these IDs were dropped (%s)", type, paste(.oId, collapse = ", "))
            .idData <- .idData[!is.na(.idData[[.w]]), ];
        }
        .idData <- .idData[order(.idData[[.w]]), ];
        .idData <- .idData[, -.w, drop = FALSE];
        if (length(.idData[, 1]) == 0) {
            if (skipStop) {
                return(.idData)
            } else {
                stop(sprintf("there are no individuals left to solve in %s data", type))
            }
        }
        if (.wrn != "") warning(.wrn)
        return(.idData)
    } else if (length(.w) == 0L) {
        if (length(idData[, 1]) > 1) {
            if (warnIdSort && .warnIdSort0 && .nid > 1)
                warning(sprintf("'ID' missing in '%s' dataset\nindividual parameters are assumed to have the same order as the event dataset", type));
        }
        return(as.data.frame(idData));
    } else {
        if (length(idData[, 1]) > 1) {
            warning(sprintf("unable to detect 'ID' correctly in '%s' dataset\nindividual parameters are assumed to have the same order as the dataset", type));
        }
        .idData <- idData[, -.w, drop = FALSE];
        return(as.data.frame(.idData));
    }
}

.convertExtra <- function(dat) {
    d <- as.data.frame(dat);
    .colNames0 <- colnames(d)
    .colNames <- toupper(.colNames0);
    ## Handle DATE TIME; DAT1 TIME; DAT2 TIME and DAT3 TIME

    ## Note NONMEM handles dates of the format DAY-MONTH and DAY as
    ## well for the DATE class of objects.

    ## It is too complex to handle, and not very common so it will
    ## throw an error

    .doDate <- FALSE;
    .dupDate <- "dates can only be specified by one of: 'DATE', 'DAT1', 'DAT2', 'DAT3' / 'TIME'"
    .checkBad <- function(d) {
        d <- paste(d);
        if (any(unlist(lapply(strsplit(d, "[^0-9]+"), length)) != 3)) {
            stop("dates formatted as MONTH-DAY or DAY alone are not supported in this conversion")
        }
        return(d)
    }
    if (any(.colNames == "DATE")) {
        ##  Month Day Year
        .datReg2 <- rex::rex(start, any_spaces, capture(numbers), non_numbers, capture(numbers), non_numbers, capture(number, number), any_spaces, end)
        .datReg4 <- rex::rex(start, any_spaces, capture(numbers), non_numbers, capture(numbers), non_numbers, capture(number, number, number, number), any_spaces, end)
        dt <- .checkBad(d$DATE)
        d$DATE.TIME <- as.POSIXct(NA);
        w <- which(regexpr(.datReg2, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(.datReg2, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format = "%m-%d-%y %H:%M")
        w <- which(regexpr(.datReg4, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(.datReg4, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format = "%m-%d-%Y %H:%M")
        d <- d[, -which(names(d) == "DATE")];
        .doDate <- TRUE;
    }
    if (any(.colNames == "DAT1")) {
        if (.doDate) {
            stop(.dupDate)
        }
        ## DAT1   day month year
        .datReg2 <- rex::rex(start, any_spaces, capture(numbers), non_numbers, capture(numbers), non_numbers, capture(number, number), any_spaces, end)
        .datReg4 <- rex::rex(start, any_spaces, capture(numbers), non_numbers, capture(numbers), non_numbers, capture(number, number, number, number), any_spaces, end)
        dt <- .checkBad(d$DAT1)
        d$DATE.TIME <- as.POSIXct(NA);
        w <- which(regexpr(.datReg2, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(.datReg2, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format = "%d-%m-%y %H:%M")
        w <- which(regexpr(.datReg4, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(.datReg4, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format = "%d-%m-%Y %H:%M")
        d <- d[, -which(names(d) == "DAT1")];
        .doDate <- TRUE;
    }
    if (any(.colNames == "DAT2")) {
        ## DAT2   year month day
        if (.doDate) {
            stop(.dupDate)
        }
        .datReg2 <- rex::rex(start, any_spaces, capture(number, number), non_numbers, capture(numbers), non_numbers, capture(numbers), any_spaces, end)
        .datReg4 <- rex::rex(start, any_spaces, capture(number, number, number, number), non_numbers, capture(numbers), non_numbers, capture(numbers), any_spaces, end)
        dt <- .checkBad(d$DAT2)
        d$DATE.TIME <- as.POSIXct(NA);
        w <- which(regexpr(.datReg2, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(.datReg2, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format = "%y-%m-%d %H:%M")
        w <- which(regexpr(.datReg4, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(.datReg4, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format = "%Y-%m-%d %H:%M")
        d <- d[, -which(names(d) == "DAT2")];
        .doDate <- TRUE;
    }
    if (any(.colNames == "DAT3")) {
        ## DAT3   year day month
        if (.doDate) {
            stop(.dupDate)
        }
        .datReg2 <- rex::rex(start, any_spaces, capture(number, number), non_numbers, capture(numbers), non_numbers, capture(numbers), any_spaces, end)
        .datReg4 <- rex::rex(start, any_spaces, capture(number, number, number, number), non_numbers, capture(numbers), non_numbers, capture(numbers), any_spaces, end)
        dt <- .checkBad(d$DAT3)
        d$DATE.TIME <- as.POSIXct(NA);
        w <- which(regexpr(.datReg2, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(.datReg2, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format = "%y-%d-%m %H:%M")
        w <- which(regexpr(.datReg4, dt) != -1)
        if (length(w) > 0)
            d$DATE.TIME[w] <- as.POSIXct(paste(gsub(.datReg4, "\\1-\\2-\\3", dt[w]), d$TIME[w]), format = "%Y-%d-%m %H:%M")
        d <- d[, -which(names(d) == "DAT3")];
        .doDate <- TRUE;
    }
    if (.doDate) {
        if (any(is.na(d$DATE.TIME))) {
            stop("date/time format was not correctly specified")
        }
    }
    if (.doDate) {
        ## Sort by date/time (though this should have been done already...)
        if (!any(names(d) == "ID")) {
            d$ID  <- 1L;
        }
        if (!any(names(d) == "EVID")) {
            d$EVID  <- 0L;
        }
        d <- d[order(d$ID, d$DATE.TIME, -d$EVID), ];
        d$TIME <- as.vector(unlist(sapply(unique(d$ID), function(id) {
            d0 <- d[d$ID == id, ];
            return(as.numeric(difftime(d0$DATE.TIME,
                                       d0$DATE.TIME[1],
                                       units = "hours")))
        })))
        d <- d[, -which(names(d) == "DATE.TIME")];
    }
    if (is(d$TIME, "numeric") || is(d$TIME, "integer")) return(d)
    stop("cannot figure out numeric time")
}

## nocov start
##'@export
print.rxEtTran <- function(x, ...) {
    print(as.data.frame(x));
    .cls <- class(x);
    .lst <- attr(.cls, ".RxODE.lst")
    cat("\nCovariates (non time-varying):\n")
    print(.lst$cov1)
    cat("\nCompartment translation:\n");
    print(data.frame("Compartment Name" = .lst$cmtInfo,
                     "Compartment Number" = seq_along(.lst$cmtInfo),
                     check.names = FALSE))
} ## nocov end
