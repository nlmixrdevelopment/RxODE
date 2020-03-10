#' @importFrom ggplot2 %+replace%
`%+replace%`
rx_theme <- function(){
    .greyText <- ggplot2::element_text(color="#808078")
    .greyLabTextX <- ggplot2::element_text(color="#808078", face="bold")
    .greyLabTextY <- ggplot2::element_text(color="#808078", face="bold", angle=90)
    .title <- ggplot2::element_text(colour = "#808078", face="bold", hjust=0)
    .subTitle <- ggplot2::element_text(colour = "#808078", face="bold", hjust=0)
    .greyTick <- ggplot2::element_line(color="#808078")
    .greyMajor <- ggplot2::element_line(color="#BFBFB4")
    .greyMinor <- ggplot2::element_line(color="#E6E6D8")
    .theme <- ggplot2::theme_bw() %+replace%
        ggplot2::theme(plot.title = .title,
                       plot.subtitle = .title,
                       panel.border = ggplot2::element_blank(),
                       ## panel.background = ggplot2::element_rect(fill = "#FFFFF7", colour = NA),
                       panel.grid.minor=.greyMinor,
                       panel.grid.major=.greyMajor,
                       axis.text.x=.greyText,
                       axis.text.y=.greyText,
                       axis.title.x=.greyLabTextX,
                       axis.title.y=.greyLabTextY,
                       axis.ticks.x=.greyTick,
                       axis.ticks.y=.greyTick,
                       strip.text=ggplot2::element_text(color="#FFFFF7", size=14, face="bold"),
                       strip.background =ggplot2::element_rect(fill="#808078", color=NA)
                       )
    ## If above ggplot2 3.2.1 then add plot.title.position="plot"
    if (utils::packageVersion("ggplot2") > "3.2.1"){
        .theme <- .theme %+replace% ggplot2::theme(plot.title.position="plot")
    }
}

##'@export
plot.rxSolve <- function(x,y,..., log="") {
    .data <- NULL
    .call <- as.list(match.call()[-(1:3)])
    .call <- .call[!names(.call) %in% c("x", "y", "log")]
    .cmts <- c(as.character(substitute(y)),
               names(sapply(as.character(.call),`c`)))
    if (length(.cmts)==1 &&.cmts[1]=="") {
        .cmts <- NULL
    }
    .dat <- rxStack(x,.cmts);
    .nlvl <- 1L
    if (any(names(.dat) == "id")) {
        if (any(names(.dat) == "sim.id")){
            .dat$id <- factor(paste0("id=", .dat$id, ", sim.id=", .dat$sim.id))
        } else {
            .dat$id <- factor(.dat$id)
        }
        .nlvl <- length(levels(.dat$id))
        .dat2 <- .dat[rev(seq_along(.dat$id)), ];
        .dat2$label <- .dat$id
        .dat2$time <- units::drop_units(.dat2$time)
        row.names(.dat2) <- NULL
        .dat2 <- .dat2[!duplicated(paste0(.dat2$id, .dat2$trt)), ];
        .aes <- aes(.data$time, .data$value, color=.data$id)
        .aesG <- aes(.data$time, .data$value, group=.data$id)
        .aesLab <- aes(label=.data$label);
    } else if (any(names(.dat) == "sim.id")){
        .dat$sim.id <- factor(.dat$sim.id)
        .nlvl <- length(levels(.dat$sim.id))
        .dat2 <- .dat[rev(seq_along(.dat$sim.id)), ];
        .dat2$label <- .dat$id
        .dat2$time <- units::drop_units(.dat2$time)
        row.names(.dat2) <- NULL
        .dat2 <- .dat2[!duplicated(paste0(.dat2$sim.id, .dat2$trt)), ];
        .aes <- aes(.data$time, .data$value, color=.data$sim.id)
        .aesG <- aes(.data$time, .data$value, group=.data$sim.id)
        .aesLab <- aes(label=.data$label);
    } else {
        .aes <- aes(.data$time, .data$value)
    }
    .facet <- facet_wrap( ~ trt, scales="free_y")
    if (length(.cmts) == 1) .facet <- NULL
    .ylab <- ylab("")
    .theme <- rx_theme()
    if (!getOption("RxODE.theme_bw", TRUE)) .theme <- NULL
    .repel <- NULL
    .legend <- NULL
    .line <- geom_line(size=1.2)
    .rxSpaghetti <- getOption("RxODE.spaghetti", 7L)
    .ggrepel <- getOption("RxODE.ggrepel", TRUE) &&
        requireNamespace("ggrepel", quietly = TRUE)
    if (.nlvl > 1 && .nlvl < .rxSpaghetti && .ggrepel && is.null(.facet)) {
        .repel <- ggrepel::geom_label_repel(.aesLab, data=.dat2, nudge_x=1,
                                            fontface="bold", size=5)
        .legend <- ggplot2::guides(color = FALSE)
    } else {
        if (.nlvl < .rxSpaghetti) {
            .legend <- ggplot2::theme(legend.title=ggplot2::element_blank())
        } else {
            .legend <-ggplot2::guides(color = FALSE)
            .line <- geom_line(size=1.2, alpha=0.2)
            .aes <- .aesG
        }
    }

    .logx <- NULL
    .logy <- NULL
    .xgxr <- getOption("RxODE.xgxr", TRUE) &&
        requireNamespace("xgxr", quietly = TRUE)
    .timex <- NULL
    .unit <- NULL
    .xgxrT <- function(.unit){
        if (!.xgxr) return(NULL)
        .timex <- xgxr::xgx_scale_x_time_units(.unit)
        if (inherits(.timex, "list")){
            .w <- which(sapply(seq_along(.timex), function(x){
                       inherits(.timex[[x]], "labels")
                   }))
            if (length(.w) > 0){
                .timex <- .timex[-.w]
            }
        }
        return(.timex)
    }
    if (inherits(.dat$time, "units")) {
        .unit <- as.character(units(.dat$time))
        .dat$time <- units::drop_units(.dat$time)
        .timex <- .xgxrT(.unit)
        .xlab <- xlab(sprintf("Time [%s]", .unit))
    } else {
        .timex <- .xgxrT("h")
        .xlab <- xlab("Time")
    }
    if (is.character(log) && length(log) == 1){
        if (length(.cmts) == 2) .facet <- NULL
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
        } else if (log != ""){
            stop(sprintf("'log=\"%s\"' not supported", log))
        }
    }
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
