#' rxTheme is the RxODE theme for plots
#'
#' @inheritParams ggplot2::theme_grey
#'
#' @param grid a Boolean indicating if the grid is on (`TRUE`) or off
#'   (`FALSE`). This could also be a character indicating `x` or `y`.
#'
#' @return ggplot2 theme used in RxODE
#'
#' @export
rxTheme <- function(base_size = 11, base_family = "",
                    base_line_size = base_size / 22,
                    base_rect_size = base_size / 22,
                    grid = TRUE) {
  half_line <- base_size / 2
  .greyTextAxisX <- ggplot2::element_text(
    color = "#808078",
    margin = ggplot2::margin(t = 0.8 * half_line / 2), vjust = 1
  )
  .greyTextAxisY <- ggplot2::element_text(
    color = "#808078",
    margin = ggplot2::margin(r = 0.8 * half_line / 2), hjust = 1
  )
  .greyLabTextX <- ggplot2::element_text(
    color = "#808078", face = "bold",
    margin = ggplot2::margin(t = half_line / 2),
    size = ggplot2::rel(1.1)
  )
  .greyLabTextY <- ggplot2::element_text(
    color = "#808078", face = "bold", angle = 90,
    margin = ggplot2::margin(r = half_line / 2),
    size = ggplot2::rel(1.1)
  )
  .title <- ggplot2::element_text(
    colour = "#808078", face = "bold", hjust = 0,
    size = ggplot2::rel(1.2),
    margin = ggplot2::margin(b = half_line)
  )
  .subTitle <- ggplot2::element_text(
    colour = "#808078", face = "bold", hjust = 0,
    margin = ggplot2::margin(b = half_line)
  )

  .greyTick <- ggplot2::element_line(color = "#808078")
  .panelGrid <-
    .greyX <- NULL
  .greyY <- NULL
  .blankGrid <- NULL
  if (inherits(grid, "character") | grid == TRUE) {
    .greyMajor <- ggplot2::element_line(color = "#BFBFB4")
    .panelGrid <- .greyMajor
    .greyMinor <- ggplot2::element_line(color = "#E6E6D8")
    .greyMajorX <- .greyMajor
    .greyMinorX <- .greyMinor
    .greyMajorY <- .greyMajor
    .greyMinorY <- .greyMinor
    if (inherits(grid, "character")) {
      if (regexpr("X", grid)[1] < 0) .greyMajorX <- ggplot2::element_blank()
      if (regexpr("Y", grid)[1] < 0) .greyMajorY <- ggplot2::element_blank()
      if (regexpr("x", grid)[1] < 0) .greyMinorX <- ggplot2::element_blank()
      if (regexpr("y", grid)[1] < 0) .greyMinorY <- ggplot2::element_blank()
    }
  } else {
    .panelGrid <- ggplot2::element_blank()
    .greyMajor <- .panelGrid
    .greyMinor <- .panelGrid
    .greyMajorX <- .greyMajor
    .greyMinorX <- .greyMinor
    .greyMajorY <- .greyMajor
    .greyMinorY <- .greyMinor
  }
  .theme <- ggplot2::theme_bw(
    base_size = base_size, base_family = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size
  ) %+replace%
    ggplot2::theme(
      plot.title = .title,
      plot.subtitle = .title,
      panel.border = ggplot2::element_blank(),
      ## panel.background = ggplot2::element_rect(fill = "#FFFFF7", colour = NA),
      panel.grid = .panelGrid,
      panel.grid.minor = .greyMinor,
      panel.grid.minor.x = .greyMinorX,
      panel.grid.minor.y = .greyMinorY,
      panel.grid.major.x = .greyMajorX,
      panel.grid.major.y = .greyMajorY,
      panel.grid.major = .greyMajor,
      axis.text.x = .greyTextAxisX,
      axis.text.y = .greyTextAxisY,
      axis.title.x = .greyLabTextX,
      axis.title.y = .greyLabTextY,
      axis.ticks.x = .greyTick,
      axis.ticks.y = .greyTick,
      strip.text = ggplot2::element_text(
        color = "#FFFFF7", face = "bold",
        size = ggplot2::rel(1.2),
        margin = ggplot2::margin(0.5 * half_line, 0.5 * half_line, 0.5 * half_line, 0.5 * half_line)
      ),
      strip.background = ggplot2::element_rect(fill = "#808078", color = NA)
    )
  ## If above ggplot2 3.2.1 then add plot.title.position="plot"
  if (utils::packageVersion("ggplot2") > "3.2.1") {
    .theme <- .theme %+replace% ggplot2::theme(plot.title.position = "plot")
  }
}

.dropUnits <- function(x) {
  if (inherits(x, "units")) {
    return(units::drop_units(x))
  }
  return(x)
}

.plotTime <- function(.dat, xlab) {
  .xgxr <- getOption("RxODE.xgxr", TRUE) &&
    requireNamespace("xgxr", quietly = TRUE)
  .timex <- NULL
  .unit <- NULL
  .xgxrT <- function(.unit) {
    if (!.xgxr) {
      return(NULL)
    }
    .timex <- xgxr::xgx_scale_x_time_units(.unit)
    if (inherits(.timex, "list")) {
      .w <- which(sapply(seq_along(.timex), function(x) {
        inherits(.timex[[x]], "labels")
      }))
      if (length(.w) > 0) {
        .timex <- .timex[-.w]
      }
    }
    return(.timex)
  }
  if (inherits(.dat$time, "units")) {
    .unit <- as.character(units(.dat$time))
    .dat$time <- .dropUnits(.dat$time)
    .timex <- .xgxrT(.unit)
    .xlab <- xlab(sprintf("%s [%s]", xlab, .unit))
  } else {
    .timex <- .xgxrT("h")
    .xlab <- xlab(xlab)
  }
  list(.timex, .xlab, .dat)
}

.plotLog <- function(.dat, .timex, log = "") {
  .logx <- NULL
  .logy <- NULL
  .xgxr <- getOption("RxODE.xgxr", TRUE) &&
    requireNamespace("xgxr", quietly = TRUE)
  if (is.character(log) && length(log) == 1) {
    if (log == "x") {
      .dat <- .dat[.dat$time > 0, ]
      if (.xgxr) {
        .logx <- xgxr::xgx_scale_x_log10()
      } else {
        .logx <- ggplot2::scale_x_log10()
      }
      .timex <- NULL
    } else if (log == "y") {
      if (.xgxr) {
        .logy <- xgxr::xgx_scale_y_log10()
      } else {
        .logy <- ggplot2::scale_y_log10()
      }
    } else if (log == "xy" || log == "yx") {
      if (.xgxr) {
        .logx <- xgxr::xgx_scale_x_log10()
        .logy <- xgxr::xgx_scale_y_log10()
      } else {
        .logx <- ggplot2::scale_x_log10()
        .logy <- ggplot2::scale_y_log10()
      }
      .dat <- .dat[.dat$time > 0, ]
      .timex <- NULL
    } else if (log != "") {
      stop(sprintf("'log=\"%s\"' not supported", log))
    }
  }
  return(list(.timex, .logx, .logy, .dat))
}

#' @export
plot.rxSolve <- function(x, y, ..., log = "",
                         xlab = "Time", ylab = "") {
  .data <- NULL
  .y <- as.character(substitute(y))
  .call0 <- match.call()[-(1:2)]
  .call <- as.list(.call0)
  .w <- names(.call) %in% c("x", "y", "log")
  if (length(.w) > 0) {
    .call <- .call[-.w]
  }
  .cmts <- c(
    as.character(substitute(y)),
    names(sapply(as.character(.call), `c`))
  )

  .cmts <- .cmts[!duplicated(.cmts)]
  .cmts <- intersect(.cmts, c(rxState(x), rxLhs(x)))
  if (length(.cmts) == 0) {
    .cmts <- NULL
  }
  .dat <- rxStack(x, .cmts)
  .nlvl <- 1L
  if (any(names(.dat) == "id")) {
    if (any(names(.dat) == "sim.id")) {
      if (any(names(.dat) == "resetno")) {
        .dat$id <- factor(paste0("id=", .dat$id, ", sim.id=", .dat$sim.id, ", resetno=", .dat$resetno))
      } else {
        .dat$id <- factor(paste0("id=", .dat$id, ", sim.id=", .dat$sim.id))
      }
    } else {
      if (any(names(.dat) == "resetno")) {
        .dat$id <- factor(paste0("id=", .dat$id, ", resetno=", .dat$resetno))
      } else {
        .dat$id <- factor(.dat$id)
      }
    }
    .nlvl <- length(levels(.dat$id))
    .dat2 <- .dat[rev(seq_along(.dat$id)), ]
    .dat2$label <- .dat$id
    .dat2$time <- .dropUnits(.dat2$time)
    row.names(.dat2) <- NULL
    .dat2 <- .dat2[!duplicated(paste0(.dat2$id, .dat2$trt)), ]
    .aes <- aes(.data$time, .data$value, color = .data$id)
    .aesG <- aes(.data$time, .data$value, group = .data$id)
    .aesLab <- aes(label = .data$label)
  } else if (any(names(.dat) == "sim.id")) {
    if (any(names(.dat) == "resetno")) {
      .dat$sim.id <- factor(paste0("sim.id=", .dat$sim.id, ", resetno=", .dat$resetno))
    } else {
      .dat$sim.id <- factor(.dat$sim.id)
    }
    .nlvl <- length(levels(.dat$sim.id))
    .dat2 <- .dat[rev(seq_along(.dat$sim.id)), ]
    .dat2$label <- .dat$sim.id
    .dat2$time <- .dropUnits(.dat2$time)
    row.names(.dat2) <- NULL
    .dat2 <- .dat2[!duplicated(paste0(.dat2$sim.id, .dat2$trt)), ]
    .aes <- aes(.data$time, .data$value, color = .data$sim.id)
    .aesG <- aes(.data$time, .data$value, group = .data$sim.id)
    .aesLab <- aes(label = .data$label)
  } else if (any(names(.dat) == "resetno")) {
    .dat$resetno <- factor(.dat$resetno)
    .nlvl <- length(levels(.dat$resetno))
    .dat2 <- .dat[rev(seq_along(.dat$resetno)), ]
    .dat2$label <- .dat$resetno
    .dat2$time <- .dropUnits(.dat2$time)
    row.names(.dat2) <- NULL
    .dat2 <- .dat2[!duplicated(paste0(.dat2$resetno, .dat2$trt)), ]
    .aes <- aes(.data$time, .data$value, color = .data$resetno)
    .aesG <- aes(.data$time, .data$value, group = .data$resetno)
    .aesLab <- aes(label = .data$label)
  } else {
    .aes <- aes(.data$time, .data$value)
  }
  .facet <- facet_wrap(~trt, scales = "free_y")
  if (length(.cmts) == 1) .facet <- NULL
  .ylab <- ylab(ylab)
  .theme <- rxTheme()
  if (!getOption("RxODE.theme", TRUE)) .theme <- NULL
  .repel <- NULL
  .legend <- NULL
  .line <- geom_line(size = 1.2)
  .rxSpaghetti <- getOption("RxODE.spaghetti", 7L)
  .ggrepel <- getOption("RxODE.ggrepel", TRUE) &&
    requireNamespace("ggrepel", quietly = TRUE)
  if (.nlvl > 1 && .nlvl < .rxSpaghetti && .ggrepel && is.null(.facet)) {
    .repel <- ggrepel::geom_label_repel(.aesLab,
      data = .dat2, nudge_x = 1,
      fontface = "bold", size = 5
    )
    .legend <- ggplot2::guides(color = "none")
  } else {
    if (.nlvl < .rxSpaghetti) {
      .legend <- ggplot2::theme(legend.title = ggplot2::element_blank())
    } else {
      .legend <- ggplot2::guides(color = "none")
      .line <- geom_line(size = 1.2, alpha = 0.2)
      .aes <- .aesG
    }
  }
  .lst <- .plotTime(.dat, xlab)
  .timex <- .lst[[1]]
  .xlab <- .lst[[2]]
  .dat <- .lst[[3]]
  .lst <- .plotLog(.dat, .timex, log)
  .timex <- .lst[[1]]
  .logx <- .lst[[2]]
  .logy <- .lst[[3]]
  .dat <- .lst[[4]]
  ggplot(.dat, .aes) +
    .line +
    .facet +
    .theme +
    .repel +
    .timex +
    .logx +
    .logy +
    .ylab +
    .xlab +
    .legend ->
  .gg
  .gg
}


#' @export
plot.rxSolveConfint1 <- function(x, y, ..., xlab = "Time", ylab = "", log = "") {
  .data <- NULL
  .lvl <- attr(class(x), ".rx")$lvl
  .parm <- attr(class(x), ".rx")$parm
  .aes <- aes(.data$time, .data$eff)
  .facet <- NULL
  .dat <- x
  .lst <- .plotTime(.dat, xlab)
  .timex <- .lst[[1]]
  .xlab <- .lst[[2]]
  .dat <- .lst[[3]]
  .lst <- .plotLog(.dat, .timex, log)
  .timex <- .lst[[1]]
  .logx <- .lst[[2]]
  .logy <- .lst[[3]]
  .dat <- .lst[[4]]
  if (length(.parm) > 1) {
    .facet <- facet_wrap(~trt, scales = "free_y")
  }
  .line <- geom_line(size = 1.2, show.legend = !is.null(.facet))
  .theme <- NULL
  if (getOption("RxODE.theme_bw", TRUE)) {
    .theme <- rxTheme()
  }
  .ylab <- ylab(ylab)
  .dat$p <- .dat$Percentile
  .d2 <- .dat[!duplicated(.dat$Percentile), ]
  .d2 <- .d2[.d2$Percentile != "50%", ]
  .dat0 <- .dat[.dat$Percentile == "50%", ]
  .dat0$low <- .dat$eff[.dat$Percentile == .d2$Percentile[1]]
  .dat0$up <- .dat$eff[.dat$Percentile == .d2$Percentile[2]]
  .aesR <- ggplot2::aes(ymin = .data$low, ymax = .data$up)
  .ribbon <- ggplot2::geom_ribbon(.aesR, alpha = 0.5, fill = "gray50")
  ## .title <- ggplot2::ggtitle(paste0("50% [", .d2$Percentile[1], ", ",
  ##                                 .d2$Percentile[2], "]"))
  ggplot2::ggplot(.dat0, .aes) +
    .ribbon +
    .line +
    .facet +
    .timex +
    .logx +
    .logy +
    .xlab +
    .ylab +
    .theme ->
  .ret
  return(.ret)
}


#' @export
plot.rxSolveConfint2 <- function(x, y, ..., xlab = "Time", ylab = "", log = "") {
  .data <- NULL
  .lvl <- attr(class(x), ".rx")$lvl
  .parm <- attr(class(x), ".rx")$parm
  .aes <- aes(.data$time, .data$p50,
    color = .data$Percentile,
    fill = .data$Percentile,
    group = .data$p
  )
  .aesR <- ggplot2::aes_string(ymin = .lvl[1], ymax = .lvl[3])
  .facet <- NULL
  .dat <- x
  .lst <- .plotTime(.dat, xlab)
  .timex <- .lst[[1]]
  .xlab <- .lst[[2]]
  .dat <- .lst[[3]]
  .lst <- .plotLog(.dat, .timex, log)
  .timex <- .lst[[1]]
  .logx <- .lst[[2]]
  .logy <- .lst[[3]]
  .dat <- .lst[[4]]
  if (length(.parm) > 1) {
    .facet <- facet_wrap(~trt, scales = "free_y")
  }
  .line <- geom_line(size = 1.1, show.legend = FALSE)
  .theme <- NULL
  if (getOption("RxODE.theme_bw", TRUE)) {
    .theme <- rxTheme()
  }
  .ylab <- ylab(ylab)
  .dat$p <- .dat$Percentile
  .d2 <- .dat[!duplicated(.dat$Percentile), ]
  .d2 <- .d2[.d2$Percentile != "50%", ]
  .d2 <- sub("[%]", "", paste0(.d2$Percentile, collapse = " / "))
  .dat$Percentile <- paste(.dat$Percentile)
  .dat$Percentile[.dat$Percentile != "50%"] <- .d2
  .dat$Percentile <- factor(.dat$Percentile, c("50%", .d2))
  .ribbon <- ggplot2::geom_ribbon(.aesR, alpha = 0.5, col = NA, show.legend = FALSE)
  .leg1 <- ggplot2::scale_color_manual(values = c("black", "gray"))
  .leg2 <- ggplot2::scale_fill_manual(values = c("black", "gray"))
  ggplot2::ggplot(.dat, .aes) +
    .ribbon +
    .line +
    .facet +
    .timex +
    .logx +
    .logy +
    .xlab +
    .ylab +
    .leg1 +
    .leg2 +
    .theme ->
  .ret
  ## p1 <- time <- eff <-Percentile <-sim.id <-id <-p2 <-p50 <-p05 <- p95 <- . <- NULL
  ## .lvl <- attr(class(x), ".rx")$lvl
  ## .parm <- attr(class(x), ".rx")$parm
  ## .ret <- ggplot2::ggplot(x,ggplot2::aes(time,p50,col=Percentile,fill=Percentile)) +
  ##     ggplot2::geom_ribbon(ggplot2::aes_string(ymin=.lvl[1],ymax=.lvl[3]),alpha=0.5)+
  ##     ggplot2::geom_line(size=1.2);
  ## if (length(.parm) > 1){
  ##     .ret <- .ret + facet_wrap( ~ trt, scales="free_y")
  ## }
  ## if (getOption("RxODE.theme_bw", TRUE)){
  ##     .ret <- .ret + ggplot2::theme_bw()
  ## }
  return(.ret)
}
