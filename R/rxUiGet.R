#' S3 for getting information from UI model
#'
#' @param x list of (UIenvironment, exact).  UI environment is the
#'   parsed function for RxODE.  `exact` is a boolean that says if an
#'   exact match is required.
#' @param ... Other arguments
#' @return value that was requested from the UI object
#' @author Matthew Fidler
#' @export
rxUiGet <- function(x, ...) {
  if (!inherits(x, "rxUiGet")) {
    stop("object is wrong type for `rxUiGet`")
  }
  UseMethod("rxUiGet")
}
#' @export
#' @rdname rxUiGet
rxUiGet.theta <- function(x, ...) {
  .x <- x[[1]]
  .ini <- .x$iniDf
  .w <- !is.na(.ini$ntheta)
  setNames(.ini$est[.w], .ini$name[.w])
}

#' @export
#' @rdname rxUiGet
rxUiGet.omega <- function(x, ...) {
  .x <- x[[1]]
  .lotri <- eval(as.call(list(quote(`lotri`), .x$.ini)))
  if (inherits(.lotri, "matrix")) {
    attr(.lotri, "lotriEst") <- NULL
    class(.lotri) <- NULL
  } else {
    attr(.lotri, "lotriEst") <- NULL
    class(.lotri) <- class(.lotri)[-1]
  }
  .lotri
}

#' @export
#' @rdname rxUiGet
rxUiGet.muRefTable <- function(x, ...) {
  .x <- x[[1]]
  .exact <- x[[2]]
  .muRef <- get("muRefDataFrame", .x)
  if (length(.muRef$theta) == 0) return(NULL)
  .muRefCov <- get("muRefCovariateDataFrame", .x)
  if (length(.muRefCov$theta) > 0) {
    .env <- new.env(parent=emptyenv())
    lapply(seq_along(.muRefCov$theta), function(i){
      .theta <- .muRefCov$theta[i]
      .cov <- paste0(.muRefCov$covariate[i], "*", .muRefCov$covariateParameter[i])
      if (exists(.theta, .env)) {
        assign(.theta, c(get(.theta, .env), .cov), .env)
      } else {
        assign(.theta, .cov, .env)
      }
    })
    .muRef$covariates <- vapply(.muRef$theta, function(theta) {
      if (exists(theta, .env)) {
        return(paste(get(theta, .env), collapse=" + "))
      }
      return("")
    }, character(1), USE.NAMES=FALSE)
  }
  if (requireNamespace("huxtable", quietly = TRUE)) {
    .muRef <- huxtable::hux(.muRef) %>%
      ## huxtable::add_colnames() %>%
      huxtable::set_bold(row = 1, col = huxtable::everywhere, value = TRUE) %>%
      huxtable::set_position("center") %>%
      huxtable::set_all_borders(TRUE)
  }
}

#' @rdname rxUiGet
#' @export
rxUiGet.multipleEndpoint <- function(x, ...) {
  .x <- x[[1]]
  .exact <- x[[2]]
  .info <- get("predDf", .x)
  if (length(.info$cond) == 1) return(NULL)
  if (getOption("RxODE.combine.dvid", TRUE)) {
    .info <- .info[order(.info$dvid), ]
  }
  .info <- with(.info, data.frame(
    variable = paste(var, "~", ifelse(use.utf(), "\u2026", "...")),
    cmt = paste0("cmt='", cond, "' or cmt=", cmt),
    "dvid*" = ifelse(is.na(dvid), "",
                     paste0("dvid='", cond, "' or dvid=", dvid)),
    check.names = FALSE,
    stringsAsFactors=FALSE))
  if (!getOption("RxODE.combine.dvid", TRUE)) {
    .info <- .info[, names(.info) != "dvid*"]
  }
  if (requireNamespace("huxtable", quietly = TRUE)) {
    .hux <- huxtable::hux(.info) %>%
      huxtable::add_colnames() %>%
      huxtable::set_bold(row = 1, col = huxtable::everywhere, value = TRUE) %>%
      huxtable::set_position("center") %>%
      huxtable::set_all_borders(TRUE)
    if (getOption("RxODE.combine.dvid", TRUE)) {
      .hux <- .hux %>%
        huxtable::add_footnote("* If dvids are outside this range, all dvids are re-numered sequentially, ie 1,7, 10 becomes 1,2,3 etc")
    }
  } else {
    .hux <- .info
  }
  .hux
}

#' @rdname rxUiGet
#' @export
rxUiGet.funPrint <- function(x, ...) {
  .x <- x[[1]]
  .ls <- ls(.x$meta, all=TRUE)
  .ret <- vector("list", length(.ls) + 3)
  .ret[[1]] <- quote(`{`)
  for (.i in seq_along(.ls)) {
    .ret[[.i + 1]] <- eval(parse(text=paste("quote(", .ls[.i], "<-", deparse1(.x$meta[[.ls[.i]]]), ")")))
  }
  .len <- length(.ls)
  .ret[[.len + 2]] <- .x$iniFun
  .ret[[.len + 3]] <- .x$modelFun
  .ret
}

#' @export
#' @rdname rxUiGet
rxUiGet.fun <- function(x, ...) {
  .ret <- rxUiGet.funPrint(x, ...)
  .ret2 <- function(){}
  body(.ret2) <- as.call(.ret)
  .ret2
}

#' @export
#' @rdname rxUiGet
rxUiGet.ini <- function(x, ...) {
  get("iniDf", x[[1]])
}

#'@export
#' @rdname rxUiGet
rxUiGet.iniFun <- function(x, ...) {
  .x <- x[[1]]
  .arg <- class(x)[1]
  bquote(ini(.(.x$.ini)))
}

#' @export
#' @rdname rxUiGet
rxUiGet.modelFun <- function(x, ...) {
  .x <- x[[1]]
  bquote(model(.(as.call(c(quote(`{`),.x$lstExpr)))))
}

#' @export
#' @rdname rxUiGet
rxUiGet.default <- function(x, ...) {
  .arg <- class(x)[1]
  if (!exists(.arg, envir=x[[1]])) return(NULL)
  get(.arg, x[[1]])
}

#' @export
`$.rxUi` <- function(obj, arg, exact = TRUE) {
  .lst <- list(obj, exact)
  class(.lst) <- c(arg, "rxUiGet")
  rxUiGet(.lst)
}

#' @export
.DollarNames.rxUi <- function(x, pattern) {
  .cmp <- vapply(as.character(methods("rxUiGet")), function(x){substr(x,9, nchar(x))}, character(1), USE.NAMES = FALSE)
  .cmp <- c(.cmp[.cmp != "default"], ls(x, all=TRUE))
  grep(pattern, .cmp, value = TRUE)
}

