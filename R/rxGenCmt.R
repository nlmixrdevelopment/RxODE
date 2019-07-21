.genCmt0 <- function(ncmt=1, oral=FALSE){
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

.genCmtMod <- function(mod, doConst=TRUE){
    if (!mod$..solveLinB) return(mod);
    ## Generates based on what is currently on the sympy stack.
    ret <- NULL
    get.var <- function(v){
        if (exists(v, envir=mod)){
            tmp1 <- get(v, envir=mod);
            return(sprintf("%s ~ %s;", v, rxFromSE(tmp1)));
        } else {
            return(NULL);
        }
    }
    oral <- exists("rx_ka", envir=mod);
    if (oral){
        oral <- get("rx_ka", envir=mod);
        oral <- rxFromSE(oral);
        oral <- (oral != "0") && (oral != "rx_ka");
    }
    dur <- exists("rx_dur", envir=mod);
    if (dur){
        dur <- get("rx_dur", envir=mod);
        dur <- rxFromSE(dur);
        dur <- (dur != "0") && (dur != "rx_dur");
    }
    rate <- exists("rx_rate", envir=mod);
    if (rate){
        rate <- get("rx_rate", envir=mod);
        rate <- rxFromSE(rate);
        rate <- (rate != "0") && (rate != "rx_rate");
    }
    tlag <- exists("rx_tlag", envir=mod);
    if (tlag){
        tlag <- get("rx_tlag", envir=mod);
        tlag <- rxFromSE(tlag);
        tlag <- (tlag != "0") && (tlag != "rx_tlag");
    }
    tlag2 <- exists("rx_tlag2", envir=mod);
    if (tlag2){
        tlag2 <- get("rx_tlag2", envir=mod);
        tlag2 <- rxFromSE(tlag2);
        tlag2 <- (tlag2 != "0") && (tlag2 != "rx_tlag2");
    }
    f2 <- exists("rx_F2", envir=mod);
    if (f2){
        f2. <- get("rx_F2", envir=mod)
        f2 <- suppressWarnings(as.numeric(f2.));
        if (is.na(f2)) f2  <- 0;
        f2 <- (f2 != 1) && (f2. != "rx_F2")
    }
    f <- exists("rx_F", envir=mod);
    if (f){
        f. <- get("rx_F", envir=mod)
        f <- suppressWarnings(as.numeric(f.));
        if (is.na(f)) f  <- 0;
        f <- (f != 1) && (f. != "rx_F")
    }
    if (exists("rx_k13", envir=mod)){
        extra <- .genCmt0(3, oral);
    } else if (exists("rx_k12", envir=mod)){
        extra <- .genCmt0(2, oral);
    } else if (exists("rx_k", envir=mod)){
        extra <- .genCmt0(1, oral);
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
        tmp1 <- mod1$rx_tlag;
        tmp1 <- rxFromSE(tmp1);
        if (oral){
            extra <- c(extra,
                       sprintf("lag(depot) = %s;", tmp1));
        } else {
            extra <- c(extra,
                       sprintf("lag(central) = %s;", tmp1));
        }
    }
    if (tlag2){
        tmp1 <- mod$rx_tlag2;
        tmp1 <- rxFromSE(tmp1);
        if (oral){
            extra <- c(extra,
                       sprintf("lag(central) = %s;", tmp1));
        } else {
            stop("Cannot handle this tlag combination.");
        }
    }
    if (f){
        tmp1 <- mod$rx_F
        tmp1 <- rxFromSE(tmp1);
        if (oral){
            extra <- c(extra,
                       sprintf("F(depot) = %s;", tmp1));
        } else {
            extra <- c(extra,
                       sprintf("F(central) = %s;", tmp1));
        }
    }
    if (f2){
        tmp1 <- mod$rx_F2;
        tmp1 <- rxFromSE(tmp1);
        if (oral){
            extra <- c(extra,
                       sprintf("F(central) = %s;", tmp1));
        } else {
            stop("Cannot handle this F(.) combination.");
        }
    }
    if (dur){
        tmp1 <- mod$rx_dur;
        tmp1 <- rxFromSE(tmp1);
        extra <- c(extra,
                   sprintf("dur(central) = %s;", tmp1));
    }
    if (rate){
        tmp1 <- mod$rx_rate
        tmp1 <- rxFromSE(tmp1);
        extra <- c(extra,
                   sprintf("rate(central) = %s;", tmp1));
    }
    assignInMyNamespace(".rxAssignLinB", TRUE);
    on.exit(assignInMyNamespace(".rxAssignLinB", FALSE))
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
                       v <- rxToSE(sprintf("d/dt(%s)", cur.state));
                       v <- get(v, envir=mod)
                       v <- rxFromSE(v);
                       return(sprintf("d/dt(%s) %s %s;", cur.state, sep, v));
                   }),
                   sapply(mv.1$lhs, function(v){
                       v1  <- rxToSE(v)
                       v1 <- get(v, envir=mod);
                       v1 <- rxFromSE(v1);
                       return(sprintf("%s=%s;", v, v1));
                   })
                   ), collapse="\n")
    return(rxS(ret, doConst))
}
