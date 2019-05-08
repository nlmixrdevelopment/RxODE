##' Find the assignments in R expression
##'
##' @param x R expression
##' @return list of assigned parameters
##' @author Hadley Wickham and Matthew L. Fidler
##' @keywords internal
##' @export
findLhs <- function(x) {
    ## Modified from http://adv-r.had.co.nz/Expressions.html find_assign4
    if (is.atomic(x) || is.name(x)) {
        character()
    } else if (is.call(x)) {
        if ((identical(x[[1]], quote(`<-`)) ||
             identical(x[[1]], quote(`=`)) ||
             identical(x[[1]], quote(`~`))) &&
             is.name(x[[2]])) {
            .lhs <- as.character(x[[2]])
        } else {
            .lhs <- character()
        }
        unique(c(.lhs, unlist(lapply(x, RxODE::findLhs))))
    } else if (is.pairlist(x)) {
        unique(unlist(lapply(x, RxODE::findLhs)))
    } else {
        stop("Don't know how to handle type ", typeof(x),
             call. = FALSE)
    }
}

##' This translates the parameters specified by the model in the correct type of solving.
##'
##' @param modText model text
##' @return Translated solve for RxODE linear compartments
##' @author Matthew L. Fidler
##' @export
##' @keywords internal
rxLinCmtTrans <- function(modText){
    .vars <- c();
    .old  <- getOption("RxODE.syntax.require.ode.first", TRUE)
    if (.old){
        options(RxODE.syntax.require.ode.first=FALSE);
        on.exit(options(RxODE.syntax.require.ode.first=TRUE));
    }
    .tmpTxt  <- rex::rex(any_spaces, or("f", "F"),  any_spaces, "(",
                               any_spaces, "depot", any_spaces,
                               ")",any_spaces,or("=","<-"),any_spaces,capture(anything))
    .txt <- try({strsplit(rxNorm(modText), "\n")[[1]]}, silent = TRUE);
    if (inherits(.txt, "try-error")){
        ## nocov start
        return(modText);
        ## nocov end
    }
    .re <- rex::rex(capture(anything), boundary, "linCmt(", capture(except_any_of(")")), ")", capture(anything));
    .w <- which(regexpr(.re, .txt, perl=TRUE) != -1)
    .isDirect <- FALSE
    if (length(.w) == 0){
        ## nocov start
        return(modText);
        ## nocov end
    } else if (length(.w) == 1){
        .regIniCenter  <- rex::rex(start,any_spaces,"central",any_spaces,
                                   "(",any_spaces,"0",any_spaces, ")", any_spaces,
                                   or("=", "<-"))
        .regIniDepot  <- rex::rex(start,any_spaces,"depot",any_spaces,
                                  "(",any_spaces,"0",any_spaces, ")", any_spaces,
                                  or("=", "<-"))
        .regFdepot <- rex::rex(start, any_spaces, or("f", "F"),  any_spaces, "(",
                               any_spaces, "depot", any_spaces,
                               ")",any_spaces,or("=","<-"),any_spaces,capture(anything));
        .regFcenter <- rex::rex(start, any_spaces, or("f", "F"),  any_spaces, "(",
                                any_spaces,"central", any_spaces,
                                ")",any_spaces,or("=","<-"),any_spaces,capture(anything));

        .regLagDepot <- rex::rex(start, any_spaces, or("lag", "alag"),  any_spaces, "(",
                               any_spaces, "depot", any_spaces,
                               ")",any_spaces,or("=","<-"),any_spaces,capture(anything));
        .regLagCenter <- rex::rex(start, any_spaces, or("lag", "alag"),  any_spaces, "(",
                                   any_spaces,"central", any_spaces,
                                  ")",any_spaces,or("=","<-"),any_spaces,capture(anything));

        .regRateDepot <- rex::rex(start, any_spaces, or("rate", "r"),  any_spaces, "(",
                               any_spaces, "depot", any_spaces,
                               ")",any_spaces,or("=","<-"),any_spaces,capture(anything));
        .regRateCenter <- rex::rex(start, any_spaces, or("rate", "r"),  any_spaces, "(",
                                   any_spaces,"central", any_spaces,
                                   ")",any_spaces,or("=","<-"),any_spaces,capture(anything));

        .regDurDepot <- rex::rex(start, any_spaces, or("dur", "d"),  any_spaces, "(",
                               any_spaces, "depot", any_spaces,
                               ")",any_spaces,or("=","<-"),any_spaces,capture(anything));
        .regDurCenter <- rex::rex(start, any_spaces, or("dur", "d"),  any_spaces, "(",
                                   any_spaces,"central", any_spaces,
                                   ")",any_spaces,or("=","<-"),any_spaces,capture(anything));

        .linCmt <- gsub(.re, "\\2", .txt[.w]);
        if (.linCmt == ""){
            .tmp <- rxState(modText);
            .linCmt  <- length(.tmp);
        } else {
            .linCmt <- gsub(rex::rex(any_spaces, capture(anything), any_spaces), "\\1",
                           strsplit(.linCmt, ",")[[1]])
            .vars <- .linCmt;
            .linCmt <- suppressWarnings({as.numeric(.linCmt)});
            .vars <- .vars[is.na(.linCmt)];
            .linCmt <- .linCmt[!is.na(.linCmt)];
            if (length(.linCmt) > 1){
                stop("Can't figure out what compartment the solved system will be put into.")
            } else if (length(.linCmt) == 0){
                .tmp <- rxState(modText);
                .w <- which(.tmp=="depot");
                if (length(.w) > 0) .tmp <- .tmp[-.w];
                .w <- which(.tmp=="central");
                if (length(.w) > 0) .tmp <- .tmp[-.w];
                .linCmt <- length(.tmp);
            } else {
                .linCmt <- .linCmt - 1;
            }
        }
        .vars <- unique(c(.vars, rxLhs(modText), rxParam(modText), findLhs(eval(parse(text=sprintf("quote({%s})", rxNorm(modText)))))))
        .varsUp <- toupper(.vars);
        .reg <- rex::rex(start, "V", capture(number), end);
        .w <- which(regexpr(.reg, .varsUp) != -1);
        if (any(.varsUp == "V")){
            if (length(.w) > 0){
                .nums <- as.numeric(gsub(.reg, "\\1", .varsUp[.w]));
                .nums  <- c(.nums, max(.nums)+1:4);
                .vs <- c("V", paste0("V", .nums))
            } else {
                .vs <- c("V", paste0("V", 1:4));
            }
        } else{
            if (length(.w) > 0){
                .nums <- as.numeric(gsub(.reg, "\\1", .varsUp[.w]));
                .nums  <- c(.nums, max(.nums)+1:4);
                .vs <- paste0("V", .nums)
            } else {
                .vs <- paste0("V", 1:4);
            }
        }
        .reg <- rex::rex(start, "Q", capture(number), end);
        .w <- which(regexpr(.reg, .varsUp) != -1);
        if (any(.varsUp == "Q")){
            if (length(.w) > 0){
                .nums <- as.numeric(gsub(.reg, "\\1", .varsUp[.w]));
                .nums  <- c(.nums, max(.nums)+1:4);
                .qs <- c("Q", paste0("Q", c(.nums, max(.nums)+1:4)))
            } else {
                .qs <- c("Q", paste0("Q", 1:2));
            }
        } else{
            if (length(.w) > 0){
                .nums <- as.numeric(gsub(.reg, "\\1", .varsUp[.w]));
                .nums  <- c(.nums, max(.nums)+1:4);
                .qs <- paste0("Q", c(.nums, max(.nums)+1:4))
            } else {
                .qs <- paste0("Q", 1:3);
            }
        }
        .oral <- any(.varsUp == "KA");
        .getVar <- function(var){
            if (length(var) == 1){
                .w <- which(.varsUp == var)
                if (length(.w) == 1){
                    return(.vars[.w])
                } else {
                    stop(sprintf("Requires parameter '%s'.", var))
                }
            } else {
                i <- 1;
                while (i <= length(var)){
                    .w <- which(.varsUp == var[i]);
                    if (length(.w) == 1){
                        return(.vars[.w]);
                    }
                    i <- i + 1;
                }
                stop(sprintf("Requires one of the following parameters: '%s'.", paste(var, collapse="', '")));
            }
        }
        .lines <- c();
        if (.oral){
            ka <- .getVar("KA");
            .lines[length(.lines) + 1] <- sprintf("rx_ka ~ %s", ka);
        } else {
            .lines[length(.lines) + 1] <- sprintf("rx_ka ~ 0");
        }
        .rateDepot  <- which(regexpr(.regRateDepot,.txt)!=-1)
        if (length(.rateDepot)>=1L){
            stop("rate(depot) does not work with a solved linear system");
        }
        .rateCenter  <- which(regexpr(.regRateCenter, .txt) !=-1)
        if (length(.rateCenter)==1L){
            .tmp <- .txt[.rateCenter];
            .txt <- .txt[-.rateCenter];
            .lines[length(.lines)+1]  <- sub(.regRateCenter,"rx_rate ~ \\1", .tmp);
        } else if (length(.rateCenter)>1L) {
            stop("Can only specify rate(central) once");
        } else {
            .lines[length(.lines) + 1]  <- sprintf("rx_rate ~ 0")
        }
        .durDepot  <- which(regexpr(.regDurDepot,.txt)!=-1)
        if (length(.durDepot)>=1L){
            stop("dur(depot) does not work with a solved linear system");
        }
        .durCenter  <- which(regexpr(.regDurCenter, .txt) !=-1)
        if (length(.durCenter)==1L){
            .tmp <- .txt[.durCenter];
            .txt <- .txt[-.durCenter];
            .lines[length(.lines)+1]  <- sub(.regDurCenter,"rx_dur ~ \\1", .tmp);
        } else if (length(.durCenter)>1L) {
            stop("Can only specify dur(central) once");
        } else {
            .lines[length(.lines) + 1]  <- sprintf("rx_dur ~ 0")
        }
        .iniDepot  <- which(regexpr(.regIniDepot, .txt)!=-1);
        if (length(.iniDepot)!=0L) stop("depot(0) is not supported in the solved system, use an ODE.");
        .iniCenter  <- which(regexpr(.regIniCenter, .txt)!=-1);
        if (length(.iniCenter)!=0L) stop("central(0) is not supported in the solved system use an ODE.");
        .lagDepot  <- which(regexpr(.regLagDepot,.txt)!=-1)
        if (length(.lagDepot)==1L){
            .tmp <- .txt[.lagDepot];
            .txt <- .txt[-.lagDepot];
            if (.oral){
                .lines[length(.lines)+1]  <- sub(.regLagDepot,"rx_tlag ~ \\1", .tmp);
            } else {
                stop("lag(depot) does not exist without a depot compartment, specify a 'ka' parameter");
            }
        } else if (length(.lagDepot)>1L){
            stop("lag(depot) cannot be duplicated in a model");
        } else {
            if (.oral){
                .lines[length(.lines) + 1]  <- sprintf("rx_tlag ~ 0")
            }
        }
        .lagCenter  <- which(regexpr(.regLagCenter, .txt) !=-1)
        if (length(.lagCenter)==1L){
            .tmp <- .txt[.lagCenter];
            .txt <- .txt[-.lagCenter];
            if (.oral){
                .lines[length(.lines)+1]  <- sub(.regLagCenter,"rx_tlag2 ~ \\1", .tmp);
            } else {
                .lines[length(.lines)+1]  <- sub(.regLagCenter,"rx_tlag ~ \\1", .tmp);
                .lines[length(.lines) + 1]  <- sprintf("rx_tlag2 ~ 0")
            }
        } else if (length(.lagCenter)>1L) {
            stop("Can only specify lag(central) once.");
        } else {
            if (.oral){
                .lines[length(.lines) + 1]  <- sprintf("rx_tlag2 ~ 0")
            } else {
                .lines[length(.lines) + 1]  <- sprintf("rx_tlag ~ 0")
                .lines[length(.lines) + 1]  <- sprintf("rx_tlag2 ~ 0")
            }
        }
        .fDepot  <- which(regexpr(.regFdepot,.txt)!=-1)
        if (length(.fDepot)==1L){
            .tmp <- .txt[.fDepot];
            .txt <- .txt[-.fDepot];
            if (.oral){
                .lines[length(.lines)+1]  <- sub(.regFdepot,"rx_F ~ \\1", .tmp);
            } else {
                stop("f(depot) does not exist without a depot compartment, specify a 'ka' parameter");
            }
        } else if (length(.fDepot)>1L){
            stop("f(depot) cannot be duplicated in a model");
        } else {
            if (.oral){
                .lines[length(.lines) + 1]  <- sprintf("rx_F ~ 1")
            }
        }
        .fCenter  <- which(regexpr(.regFcenter, .txt) !=-1)
        if (length(.fCenter)==1L){
            .tmp <- .txt[.fCenter];
            .txt <- .txt[-.fCenter];
            if (.oral){
                .lines[length(.lines)+1]  <- sub(.regFcenter,"rx_F2 ~ \\1", .tmp);
            } else {
                .lines[length(.lines)+1]  <- sub(.regFcenter,"rx_F ~ \\1", .tmp);
                .lines[length(.lines) + 1]  <- sprintf("rx_F2 ~ 1")
            }
        } else if (length(.fCenter)>1L) {
            stop("Can only specify f(central) once.");
        } else {
            if (.oral){
                .lines[length(.lines) + 1]  <- sprintf("rx_F2 ~ 1")
            } else {
                .lines[length(.lines) + 1]  <- sprintf("rx_F ~ 1")
                .lines[length(.lines) + 1]  <- sprintf("rx_F2 ~ 1")
            }
        }
        .ncmt <- 1;
        if (any(.varsUp == "CL")){
            .hasVp  <- any(regexpr("^(VP[1-9]*)$", .varsUp)!=-1);
            .hasVt  <- any(regexpr("^(VT[1-9]*)$", .varsUp)!=-1);
            .hasVn  <- any(regexpr("^V[1-9]+$", .varsUp)!=-1);
            if (.hasVp && .hasVn){
                stop("Cannot mix Vp and V#");
            }
            if (.hasVt && .hasVn){
                stop("Cannot mix Vt and V#");
            }
            if (.hasVp && .hasVt){
                stop("Cannot mix Vp and Vt");
            }
            if (any(.varsUp=="VP2") && !any(.varsUp=="VP")){
                stop("Defined VP2 without VP");
            }
            if (any(.varsUp=="VT2") && !any(.varsUp=="VT")){
                stop("Defined VT2 without VT");
            }
            if (any(.varsUp=="CLD2") && !any(.varsUp=="CLD")){
                stop("Defined CLD2 without CLD");
            }
            .hasQ  <- any(regexpr("^(Q[1-9]*)$", .varsUp)!=-1);
            .hasCld  <- any(regexpr("^(CLD[1-9]*)$", .varsUp)!=-1);
            if (.hasQ && .hasCld){
                stop("Cannot mix Q and Cld")
            }
            .cl <- .getVar("CL");
            .v <- .getVar(c("V", "VC", .vs[1]));
            if (toupper(.v)=="V" && any(.varsUp=="VC")){
                stop(sprintf("Ambiguous %s/%s specification", .v, .getVar("VC")))
            }
            if (toupper(.v) !=.vs[1]){
                .vs <- c(.v, .vs);
            }
            .lines[length(.lines) + 1] <- sprintf("rx_v ~ %s", .v);
            .lines[length(.lines) + 1] <- sprintf("rx_k ~ %s/%s", .cl, .v);
            .hasVss  <- FALSE
            if (any(.varsUp == "VSS")){
                .ncmt <- 2;
                .Q <- .getVar(c("CLD", .qs[1]));
                .vss <- .getVar("VSS");
                .hasVss <- TRUE
                if (.hasVp){
                    stop("Cannot define Vp with a VSS defined");
                }
                if (.hasVt){
                    stop("Cannot define Vt with a VSS defined");
                }
                .tmp <- .varsUp[.varsUp != toupper(.v)];
                if (any(regexpr("^V[0-9]*$", .tmp)!=-1)){
                    stop("Cannot have Volumes other than Vss and central volume defined.");
                }
                .lines[length(.lines) + 1] <- sprintf("rx_k12 ~ %s/%s", .Q, .v);
                .lines[length(.lines) + 1] <- sprintf("rx_k21 ~ %s/(%s-%s)", .Q, .vss, .v);
            } else if ((any(.varsUp == .vs[2]) || any(.varsUp==.vs[3]) ||
                        any(.varsUp == "VP"))){
                .ncmt <- 2;
                .Q <- .getVar(c("CLD",.qs[1]));
                .v2 <- .getVar(c("VP", .vs[2], .vs[3]));
                if (toupper(.v2)==.vs[3]){
                    .vs  <- .vs[-2];
                }
                if (toupper(.v2) !=.vs[2]){
                    .vs  <- c(.vs[1], .v2, .vs[-1]);
                }
                if (toupper(.Q) !=.qs[1]){
                    .qs  <- c(.Q,.qs);
                }
                .lines[length(.lines) + 1] <- sprintf("rx_k12 ~ %s/%s", .Q, .v);
                .lines[length(.lines) + 1] <- sprintf("rx_k21 ~ %s/%s", .Q, .v2);
            } else if (any(.varsUp == "VT")){
                .ncmt <- 2;
                .Q <- .getVar(c("CLD", .qs[1]));
                if (toupper(.Q)=="CLD"){
                    .qs  <- c("CLD",.qs);
                }
                .v2 <- .getVar("VT");
                if (toupper(.v2) !=.vs[2]){
                    .vs  <- c(.vs[1], .v2, .vs[-1]);
                }
                .lines[length(.lines) + 1] <- sprintf("rx_k12 ~ %s/%s", .Q, .v);
                .lines[length(.lines) + 1] <- sprintf("rx_k21 ~ %s/%s", .Q, .v2);
            } else if (any(.varsUp==.qs[1])){
                stop(sprintf("Defined '%s' without corresponding volume", .qs[1]));
            } else if (any(.varsUp=="CLD")){
                stop("Defined 'CLD' without corresponding volume");
            }
            if (any(.varsUp == toupper(.vs[3]))){
                if (.hasVss){
                    stop("Vss only supported with 2 compartment models.")
                }
                .ncmt <- 3;
                .v3 <- .getVar(.vs[3]);
                .q2 <- .getVar(c("CLD2", .qs[2]));
                .lines[length(.lines) + 1] <- sprintf("rx_k13 ~ %s/%s", .q2, .v);
                .lines[length(.lines) + 1] <- sprintf("rx_k31 ~ %s/%s", .q2, .v3);
            } else if (any(.varsUp=="VP2")){
                if (.hasVss){
                    stop("Vss only supported with 2 compartment models")
                }
                .ncmt <- 3;
                .v3 <- .getVar("VP2");
                .q2 <- .getVar(c(.qs[2],"CLD2"));
                .lines[length(.lines) + 1] <- sprintf("rx_k13 ~ %s/%s", .q2, .v);
                .lines[length(.lines) + 1] <- sprintf("rx_k31 ~ %s/%s", .q2, .v3);
            } else if (any(.varsUp == "VT2")) {
                if (.hasVss){
                    stop("Vss only supported with 2 compartment models")
                }
                .ncmt <- 3;
                .v3 <- .getVar("VT2");
                .q2 <- .getVar(c(.qs[2], "CLD2"));
                .lines[length(.lines) + 1] <- sprintf("rx_k13 ~ %s/%s", .q2, .v);
                .lines[length(.lines) + 1] <- sprintf("rx_k31 ~ %s/%s", .q2, .v3);
            } else if (any(.varsUp==.qs[2])){
                if (.hasVss){
                    stop("Vss only supported with 2 compartment models")
                }
                stop(sprintf("Defined '%s' without corresponding volume", .qs[2]));
            } else if (any(.varsUp=="CLD2")){
                stop("Defined 'CLD2' without corresponding volume");
            }
        } else if (any(.varsUp == "K") || any(.varsUp == "KE") || any(.varsUp == "KEL")) {
            .k <- .getVar(c("K", "KE", "KEL"))
            .v <- .getVar(c("V", "VC", .vs[1]));
            if (toupper(.v)=="V" && any(.varsUp=="VC")){
                stop(sprintf("Ambiguous %s/%s specification", .v, .getVar("VC")))
            }
            if (toupper(.v)=="VC" && any(.varsUp==.vs[1])){
                stop(sprintf("Ambiguous %s/%s specification", .v, .getVar(.vs[1])))
            }
            if (toupper(.v) !=.vs[1]){
                .vs <- c(.v, .vs);
            }
            .lines[length(.lines) + 1] <- sprintf("rx_v ~ %s", .v);
            .lines[length(.lines) + 1] <- sprintf("rx_k ~ %s", .k);
            if (any(.varsUp == "K12") || any(.varsUp == "K21")){
                .ncmt <- 2;
                .k12 <- .getVar("K12")
                .k21 <- .getVar("K21")
                .lines[length(.lines) + 1] <- sprintf("rx_k12 ~ %s", .k12);
                .lines[length(.lines) + 1] <- sprintf("rx_k21 ~ %s", .k21);
            }
            if (any(.varsUp == "K13") || any(.varsUp=="K31")){
                if (.ncmt !=2){
                    stop("K12 and K21 need to be defined for a 3 compartment model.");
                }
                .ncmt <- 3;
                .k13 <- .getVar("K13")
                .k31 <- .getVar("K31")
                .lines[length(.lines) + 1] <- sprintf("rx_k13 ~ %s", .k13);
                .lines[length(.lines) + 1] <- sprintf("rx_k31 ~ %s", .k31);
            }
        } else if (any(.varsUp == "AOB")){
            .ncmt <- 2;
            .v <- .getVar(c("V", "VC", .vs[1]));
            if (toupper(.v)=="V" && any(.varsUp=="VC")){
                stop(sprintf("Ambiguous %s/%s specification", .v, .getVar("VC")))
            }
            if (toupper(.v)=="VC" && any(.varsUp==.vs[1])){
                stop(sprintf("Ambiguous %s/%s specification", .v, .getVar(.vs[1])))
            }
            if (toupper(.v) !=.vs[1]){
                .vs <- c(.v, .vs);
            }
            .aob <- .getVar("AOB");
            .alpha <- .getVar("ALPHA");
            .beta <- .getVar("BETA");
            .lines[length(.lines) + 1] <- sprintf("rx_v ~ %s", .v);
            .lines[length(.lines) + 1] <- sprintf("rx_k21 ~ (%s*%s+%s)/(%s+1)", .aob, .beta, .alpha, .aob);
            .lines[length(.lines) + 1] <- sprintf("rx_k ~ (%s*%s)/rx_k21", .alpha, .beta)
            .lines[length(.lines) + 1] <- sprintf("rx_k12 ~ %s+%s - rx_k21 - rx_k", .alpha, .beta);
        } else if (any(.varsUp == "ALPHA") && any(.varsUp == "BETA") && any(.varsUp == "K21")){
            .ncmt <- 2;
            .k21 <- .getVar("K21");
            .alpha <- .getVar("ALPHA");
            .beta <- .getVar("BETA");
            .v <- .getVar(c("V", "VC", .vs[1]));
            if (toupper(.v)=="V" && any(.varsUp=="VC")){
                stop(sprintf("Ambiguous %s/%s specification", .v, .getVar("VC")))
            }
            if (toupper(.v)=="VC" && any(.varsUp==.vs[1])){
                stop(sprintf("Ambiguous %s/%s specification", .v, .getVar(.vs[1])))
            }
            if (toupper(.v) !=.vs[1]){
                .vs <- c(.v, .vs);
            }
            .lines[length(.lines) + 1] <- sprintf("rx_v ~ %s", .v);
            .lines[length(.lines) + 1] <- sprintf("rx_k21 ~ %s", .k21);
            .lines[length(.lines) + 1] <- sprintf("rx_k ~ (%s*%s)/rx_k21", .alpha, .beta);
            .lines[length(.lines) + 1] <- sprintf("rx_k12 ~ %s+%s - rx_k21 - rx_k", .alpha, .beta);
        } else if (any(.varsUp == "ALPHA") && any(.varsUp == "A")) {
            .isDirect  <- TRUE
        } else if (any(.varsUp=="ALPHA")){
            .ncmt <- 1;
            .v <- .getVar(c("V", "VC", .vs[1]));
            .alpha <- .getVar("ALPHA");
            if (toupper(.v)=="V" && any(.varsUp=="VC")){
                stop(sprintf("Ambiguous %s/%s specification", .v, .getVar("VC")))
            }
            if (toupper(.v)=="VC" && any(.varsUp==.vs[1])){
                stop(sprintf("Ambiguous %s/%s specification", .v, .getVar(.vs[1])))
            }
            if (any(.varsUp=="GAMMA")){
                stop("A gamma parameter requires A/B/C and alpha/beta");
            }
            if (any(.varsUp=="K21")){
                stop("K21 requires a beta parameter.")
            }
            if (any(.varsUp=="AOB")){
                stop("AOB requires a beta parameter.")
            }
            if (any(.varsUp=="BETA")){
                stop("Beta requires AOB or K21 parameter.");
            }
            .lines[length(.lines) + 1] <- sprintf("rx_v ~ %s", .v);
            .lines[length(.lines) + 1] <- sprintf("rx_k ~ %s", .alpha);

        } else {
            stop("Could not figure out the linCmt() from the defined parameters.")
        }
        if (.isDirect){
            .alpha <- .getVar("ALPHA");
            .a  <- .getVar("A");
            if (any(.varsUp=="BETA") || any(.varsUp=="B")){
                .beta <- .getVar("BETA");
                .b <- .getVar("B");
                if (any(.varsUp=="GAMMA") || any(.varsUp=="C")){
                    ## 3 cmt
                    .gamma <- .getVar("GAMMA");
                    .c <- .getVar("C");
                    .lines[length(.lines) + 1] <- sprintf("rx_alpha ~ %s;", .alpha);
                    .lines[length(.lines) + 1] <- sprintf("rx_beta ~ %s;", .beta);
                    .lines[length(.lines) + 1] <- sprintf("rx_gamma ~ %s;", .gamma);
                    .lines[length(.lines) + 1] <- sprintf("rx_A ~ %s;",.a);
                    .lines[length(.lines) + 1] <- sprintf("rx_B ~ %s;",.b);
                    .lines[length(.lines) + 1] <- sprintf("rx_C ~ %s;",.c);
                    if (.oral){
                        .lines[length(.lines) + 1] <- "rx_A2 ~ rx_A";
                        .lines[length(.lines) + 1] <- "rx_B2 ~ rx_B";
                        .lines[length(.lines) + 1] <- "rx_C2 ~ rx_C";
                        .lines[length(.lines) + 1] <- "rx_A ~ rx_ka / (rx_ka - rx_alpha) * rx_A";
                        .lines[length(.lines) + 1] <- "rx_B ~ rx_ka / (rx_ka - rx_beta) * rx_B";
                        .lines[length(.lines) + 1] <- "rx_C ~ rx_ka / (rx_ka - rx_gamma) * rx_C";
                    } else {
                        .lines[length(.lines) + 1] <- "rx_A2 ~ 0";
                        .lines[length(.lines) + 1] <- "rx_B2 ~ 0";
                        .lines[length(.lines) + 1] <- "rx_C2 ~ 0";
                    }
                } else {
                    ## 2 cmt
                    .lines[length(.lines) + 1] <- sprintf("rx_alpha ~ %s", .alpha);
                    .lines[length(.lines) + 1] <- sprintf("rx_beta ~ %s", .beta);
                    if (.oral){
                        .lines[length(.lines) + 1] <- sprintf("rx_A ~ rx_ka / (rx_ka - rx_alpha) * %s", .a);
                        .lines[length(.lines) + 1] <- sprintf("rx_B ~ rx_ka / (rx_ka - rx_beta) * %s;", .b);
                        .lines[length(.lines) + 1] <- sprintf("rx_A2 ~ %s", .a)
                        .lines[length(.lines) + 1] <- sprintf("rx_B2 ~ %s", .b);
                    } else {
                        .lines[length(.lines) + 1] <-  sprintf("rx_A ~ %s;", .a)
                        .lines[length(.lines) + 1] <- sprintf("rx_B ~ %s;", .b);
                        .lines[length(.lines) + 1] <-  "rx_A2 ~ 0"
                        .lines[length(.lines) + 1] <- "rx_B2 ~ 0";
                    }
                    .lines[length(.lines) + 1] <- "rx_gamma ~ 0";
                    .lines[length(.lines) + 1] <- "rx_C ~ 0";
                    .lines[length(.lines) + 1]  <- "rx_C2 ~ 0";

                }
            } else {
                if (any(.varsUp=="GAMMA") || any(.varsUp=="C")){
                    stop("A three compartment model requires BETA/B");
                }
                ## 1 cmt
                .lines[length(.lines) + 1] <- sprintf("rx_alpha ~ %s", .alpha);
                if (.oral){
                    .lines[length(.lines) + 1] <- sprintf("rx_A ~ rx_ka / (rx_ka - rx_alpha) * %s", .a);
                    .lines[length(.lines) + 1] <- sprintf("rx_A2 ~ %s", .a);
                } else {
                    .lines[length(.lines) + 1] <- sprintf("rx_A ~ %s", .a);
                    .lines[length(.lines) + 1] <- "rx_A2 ~ 0.0";
                }

                .lines[length(.lines) + 1] <- "rx_beta ~ 0";
                .lines[length(.lines) + 1] <- "rx_B ~ 0";
                .lines[length(.lines) + 1] <- "rx_B2 ~ 0";
                .lines[length(.lines) + 1] <- "rx_gamma ~ 0";
                .lines[length(.lines) + 1] <- "rx_C ~ 0";
                .lines[length(.lines) + 1] <- "rx_C2 ~ 0";
            }
        } else {
            if (.ncmt == 1){
                .lines[length(.lines) + 1] <- "rx_alpha ~ rx_k";
                if (.oral){
                    .lines[length(.lines) + 1] <- "rx_A ~ rx_ka / (rx_ka - rx_alpha) / rx_v";
                    .lines[length(.lines) + 1] <- "rx_A2 ~ 1.0 / rx_v";
                } else {
                    .lines[length(.lines) + 1] <- "rx_A ~ 1.0 / rx_v";
                    .lines[length(.lines) + 1] <- "rx_A2 ~ 0.0";
                }
                .lines[length(.lines) + 1] <- "rx_beta ~ 0";
                .lines[length(.lines) + 1] <- "rx_B ~ 0";
                .lines[length(.lines) + 1] <- "rx_B2 ~ 0";
                .lines[length(.lines) + 1] <- "rx_gamma ~ 0";
                .lines[length(.lines) + 1] <- "rx_C ~ 0";
                .lines[length(.lines) + 1] <- "rx_C2 ~ 0";
            } else if (.ncmt == 2){
                .lines[length(.lines) + 1] <- "rx_beta  ~ 0.5 * (rx_k12 + rx_k21 + rx_k - sqrt((rx_k12 + rx_k21 + rx_k) * (rx_k12 + rx_k21 + rx_k) - 4.0 * rx_k21 * rx_k))"
                .lines[length(.lines) + 1] <-  "rx_alpha ~ rx_k21 * rx_k / rx_beta"
                if (.oral){
                    .lines[length(.lines) + 1] <-  "rx_A ~ rx_ka / (rx_ka - rx_alpha) * (rx_alpha - rx_k21) / (rx_alpha - rx_beta) / rx_v"
                    .lines[length(.lines) + 1] <- "rx_B ~ rx_ka / (rx_ka - rx_beta) * (rx_beta - rx_k21) / (rx_beta - rx_alpha) / rx_v;"
                    .lines[length(.lines) + 1] <-  "rx_A2 ~ (rx_alpha - rx_k21) / (rx_alpha - rx_beta) / rx_v"
                    .lines[length(.lines) + 1] <- "rx_B2 ~ (rx_beta - rx_k21) / (rx_beta - rx_alpha) / rx_v;"
                } else {
                    .lines[length(.lines) + 1] <-  "rx_A ~ (rx_alpha - rx_k21) / (rx_alpha - rx_beta) / rx_v"
                    .lines[length(.lines) + 1] <- "rx_B ~ (rx_beta - rx_k21) / (rx_beta - rx_alpha) / rx_v;"
                    .lines[length(.lines) + 1] <-  "rx_A2 ~ 0"
                    .lines[length(.lines) + 1] <- "rx_B2 ~ 0";
                }
                .lines[length(.lines) + 1] <- "rx_gamma ~ 0";
                .lines[length(.lines) + 1] <- "rx_C ~ 0";
                .lines[length(.lines) + 1]  <- "rx_C2 ~ 0";
            } else if (.ncmt == 3){
                .lines[length(.lines) + 1] <- "rx_a0 ~ rx_k * rx_k21 * rx_k31";
                .lines[length(.lines) + 1] <- "rx_a1 ~ rx_k * rx_k31 + rx_k21 * rx_k31 + rx_k21 * rx_k13 + rx_k * rx_k21 + rx_k31 * rx_k12";
                .lines[length(.lines) + 1] <- "rx_a2 ~ rx_k + rx_k12 + rx_k13 + rx_k21 + rx_k31"
                .lines[length(.lines) + 1] <- "rx_p ~ rx_a1 - rx_a2 * rx_a2 / 3.0";
                .lines[length(.lines) + 1] <- "rx_q ~ 2.0 * rx_a2 * rx_a2 * rx_a2 / 27.0 - rx_a1 * rx_a2 /3.0 + rx_a0"
                .lines[length(.lines) + 1] <- "rx_r1 ~ sqrt(-rx_p * rx_p * rx_p / 27.0)";
                .lines[length(.lines) + 1] <- "rx_r2 ~ 2 * rx_r1^(1.0/3.0)";
                .lines[length(.lines) + 1] <- "rx_theta ~ acos(-rx_q / (2.0 * rx_r1)) / 3.0"
                .lines[length(.lines) + 1] <- "rx_alpha ~ -(cos(rx_theta) * rx_r2 - rx_a2 / 3.0)"
                .lines[length(.lines) + 1] <- "rx_beta ~ -(cos(rx_theta + 2.0 / 3.0 * pi) * rx_r2 - rx_a2 / 3.0)"
                .lines[length(.lines) + 1] <- "rx_gamma ~ -(cos(rx_theta + 4.0 / 3.0 * pi) * rx_r2 - rx_a2 / 3.0)";
                .lines[length(.lines) + 1] <- "rx_A ~ (rx_k21 - rx_alpha) * (rx_k31 - rx_alpha) / (rx_alpha - rx_beta) / (rx_alpha - rx_gamma) / rx_v;"
                .lines[length(.lines) + 1] <- "rx_B ~ (rx_k21 - rx_beta) * (rx_k31 - rx_beta) / (rx_beta - rx_alpha) / (rx_beta - rx_gamma) / rx_v;"
                .lines[length(.lines) + 1] <- "rx_C ~ (rx_k21 - rx_gamma) * (rx_k31 - rx_gamma) / (rx_gamma - rx_alpha) / (rx_gamma - rx_beta) / rx_v;"
                if (.oral){
                    .lines[length(.lines) + 1] <- "rx_A2 ~ rx_A";
                    .lines[length(.lines) + 1] <- "rx_B2 ~ rx_B";
                    .lines[length(.lines) + 1] <- "rx_C2 ~ rx_C";
                    .lines[length(.lines) + 1] <- "rx_A ~ rx_ka / (rx_ka - rx_alpha) * rx_A";
                    .lines[length(.lines) + 1] <- "rx_B ~ rx_ka / (rx_ka - rx_beta) * rx_B";
                    .lines[length(.lines) + 1] <- "rx_C ~ rx_ka / (rx_ka - rx_gamma) * rx_C";
                } else {
                    .lines[length(.lines) + 1] <- "rx_A2 ~ 0";
                    .lines[length(.lines) + 1] <- "rx_B2 ~ 0";
                    .lines[length(.lines) + 1] <- "rx_C2 ~ 0";
                }
            }
        }
        .solve <- sprintf("solveLinB(rx__PTR__, t, %s, rx_A, rx_A2, rx_alpha, rx_B, rx_B2, rx_beta, rx_C, rx_C2, rx_gamma, rx_ka, rx_tlag, rx_tlag2, rx_F, rx_F2, rx_rate, rx_dur)", .linCmt);
        .lines <- paste(.lines, collapse="\n");
        .txt <- paste(sub(.re, sprintf("%s\n\\1%s\\3", .lines, .solve), .txt), collapse="\n");
        ## Put in extra compartment information
        return(rxGetModel(.txt))
    } else {
        stop("Can only have one linCmt() function in the model.  Assign it to a variable if you need the concentrations more than once.");
    }
}
