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
        stop(sprintf("do not know how to handle type '%s'", typeof(x)),
             call. = FALSE)
    }
}

##' This translates the parameters specified by the model in the correct type of solving.
##'
##' @param modText model text
##' @return Translated solve for RxODE linear compartments
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
rxLinCmtTrans <- function(modText, linCmtSens=FALSE){
    .derived <- is.na(linCmtSens)
    if (.derived) linCmtSens <- FALSE
    .vars <- c()
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
                stop("cannot figure out the solved system compartment")
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
                    stop(sprintf("requires parameter '%s'", var))
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
                stop(sprintf("requires one of the following parameters: '%s'", paste(var, collapse="', '")));
            }
        }
        .lines <- c();
        if (.oral){
            ka <- .getVar("KA");
            .lines[length(.lines) + 1] <- sprintf("rx_ka %s %s", ifelse(.derived, "=", "~"), ka);
        } else {
            .lines[length(.lines) + 1] <- sprintf("rx_ka %s 0", ifelse(.derived, "=", "~"));
        }

        .rateDepot  <- which(regexpr(.regRateDepot,.txt)!=-1)
        if (length(.rateDepot)==1L){
          .tmp <- .txt[.rateDepot];
          .txt <- .txt[-.rateDepot];
          if (.oral){
            .lines[length(.lines)+1]  <- sub(.regRateDepot,"rx_rate %s \\1", ifelse(.derived, "=", "~"), .tmp);
          } else {
            stop("'rate(depot)' does not exist without a 'depot' compartment, specify a 'ka' parameter");
          }
        } else if (length(.rateDepot)>1L){
          stop("'rate(depot)' cannot be duplicated in a model");
        } else {
          if (.oral){
            .lines[length(.lines) + 1]  <- sprintf("rx_rate %s 0", ifelse(.derived, "=", "~"))
          }
        }
        .rateCenter  <- which(regexpr(.regRateCenter, .txt) !=-1)
        if (length(.rateCenter)==1L){
          .tmp <- .txt[.rateCenter];
          .txt <- .txt[-.rateCenter];
          if (.oral){
            .lines[length(.lines)+1]  <- sub(.regRateCenter,"rx_rate2 %s \\1",
                                             ifelse(.derived, "=", "~"), .tmp);
          } else {
            .lines[length(.lines)+1]  <- sub(.regRateCenter,"rx_rate %s \\1",
                                             ifelse(.derived, "=", "~"), .tmp);
            .lines[length(.lines) + 1]  <- sprintf("rx_rate2 %s 0", ifelse(.derived, "=", "~"))
          }
        } else if (length(.rateCenter)>1L) {
          stop("can only specify 'rate(central)' once")
        } else {
          if (.oral){
            .lines[length(.lines) + 1]  <- sprintf("rx_rate2 %s 0", ifelse(.derived, "=", "~"))
          } else {
            .lines[length(.lines) + 1]  <- sprintf("rx_rate %s 0", ifelse(.derived, "=", "~"))
            .lines[length(.lines) + 1]  <- sprintf("rx_rate2 %s 0", ifelse(.derived, "=", "~"))
          }
        }
        ## dur Center
        .durDepot  <- which(regexpr(.regDurDepot,.txt)!=-1)
        if (length(.durDepot)==1L){
          .tmp <- .txt[.durDepot];
          .txt <- .txt[-.durDepot];
          if (.oral){
            .lines[length(.lines)+1]  <- sub(.regDurDepot,"rx_dur %s \\1", ifelse(.derived, "=", "~"), .tmp);
          } else {
            stop("'f(depot)' does not exist without a 'depot' compartment, specify a 'ka' parameter");
          }
        } else if (length(.durDepot)>1L){
          stop("'dur(depot)' cannot be duplicated in a model")
        } else {
          if (.oral){
            .lines[length(.lines) + 1]  <- sprintf("rx_dur %s 0", ifelse(.derived, "=", "~"))
          }
        }
        .durCenter  <- which(regexpr(.regDurCenter, .txt) !=-1)
        if (length(.durCenter)==1L){
          .tmp <- .txt[.durCenter];
          .txt <- .txt[-.durCenter];
          if (.oral){
            .lines[length(.lines)+1]  <- sub(.regDurCenter,"rx_dur2 %s \\1", ifelse(.derived, "=", "~"), .tmp);
          } else {
            .lines[length(.lines)+1]  <- sub(.regDurCenter,"rx_dur %s \\1", ifelse(.derived, "=", "~"), .tmp);
            .lines[length(.lines) + 1]  <- sprintf("rx_dur2 %s 0", ifelse(.derived, "=", "~"))
          }
        } else if (length(.durCenter)>1L) {
          stop("can only specify 'dur(central)' once");
        } else {
          if (.oral){
            .lines[length(.lines) + 1]  <- sprintf("rx_dur2 %s 0", ifelse(.derived, "=", "~"))
          } else {
            .lines[length(.lines) + 1]  <- sprintf("rx_dur %s 0", ifelse(.derived, "=", "~"))
            .lines[length(.lines) + 1]  <- sprintf("rx_dur2 %s 0", ifelse(.derived, "=", "~"))
          }
        }
        .iniDepot  <- which(regexpr(.regIniDepot, .txt)!=-1);
        if (length(.iniDepot)!=0L) stop("'depot(0)' is not supported in the solved system, use an ODE");
        .iniCenter  <- which(regexpr(.regIniCenter, .txt)!=-1);
        if (length(.iniCenter)!=0L) stop("'central(0)' is not supported in the solved system, use an ODE");
        .lagDepot  <- which(regexpr(.regLagDepot,.txt)!=-1)
        if (length(.lagDepot)==1L){
            .tmp <- .txt[.lagDepot];
            .txt <- .txt[-.lagDepot];
            if (.oral){
                .lines[length(.lines)+1]  <- sub(.regLagDepot,"rx_tlag %s \\1", ifelse(.derived, "=", "~"), .tmp);
            } else {
                stop("'alag(depot)' does not exist without a 'depot' compartment, specify a 'ka' parameter");
            }
        } else if (length(.lagDepot)>1L){
            stop("'alag(depot)' cannot be duplicated in a model");
        } else {
            if (.oral){
                .lines[length(.lines) + 1]  <- sprintf("rx_tlag %s 0", ifelse(.derived, "=", "~"))
            }
        }
        .lagCenter  <- which(regexpr(.regLagCenter, .txt) !=-1)
        if (length(.lagCenter)==1L){
            .tmp <- .txt[.lagCenter];
            .txt <- .txt[-.lagCenter];
            if (.oral){
                .lines[length(.lines)+1]  <- sub(.regLagCenter,"rx_tlag2 %s \\1", ifelse(.derived, "=", "~"), .tmp);
            } else {
                .lines[length(.lines)+1]  <- sub(.regLagCenter,"rx_tlag %s \\1", ifelse(.derived, "=", "~"), .tmp);
                .lines[length(.lines) + 1]  <- sprintf("rx_tlag2 %s 0", ifelse(.derived, "=", "~"))
            }
        } else if (length(.lagCenter)>1L) {
            stop("can only specify 'alag(central)' once");
        } else {
            if (.oral){
                .lines[length(.lines) + 1]  <- sprintf("rx_tlag2 %s 0", ifelse(.derived, "=", "~"))
            } else {
                .lines[length(.lines) + 1]  <- sprintf("rx_tlag %s 0", ifelse(.derived, "=", "~"))
                .lines[length(.lines) + 1]  <- sprintf("rx_tlag2 %s 0", ifelse(.derived, "=", "~"))
            }
        }
        .fDepot  <- which(regexpr(.regFdepot,.txt)!=-1)
        if (length(.fDepot)==1L){
            .tmp <- .txt[.fDepot];
            .txt <- .txt[-.fDepot];
            if (.oral){
                .lines[length(.lines)+1]  <- sub(.regFdepot,"rx_F %s \\1", ifelse(.derived, "=", "~"), .tmp);
            } else {
                stop("'f(depot)' does not exist without a 'depot' compartment, specify a 'ka' parameter");
            }
        } else if (length(.fDepot)>1L){
            stop("'f(depot)' cannot be duplicated in a model");
        } else {
            if (.oral){
                .lines[length(.lines) + 1]  <- sprintf("rx_F %s 1", ifelse(.derived, "=", "~"))
            }
        }
        .fCenter  <- which(regexpr(.regFcenter, .txt) !=-1)
        if (length(.fCenter)==1L){
            .tmp <- .txt[.fCenter];
            .txt <- .txt[-.fCenter];
            if (.oral){
                .lines[length(.lines)+1]  <- sub(.regFcenter,"rx_F2 %s \\1", ifelse(.derived, "=", "~"), .tmp);
            } else {
                .lines[length(.lines)+1]  <- sub(.regFcenter,"rx_F %s \\1", ifelse(.derived, "=", "~"), .tmp);
                .lines[length(.lines) + 1]  <- sprintf("rx_F2 %s 1", ifelse(.derived, "=", "~"))
            }
        } else if (length(.fCenter)>1L) {
            stop("can only specify 'f(central)' once");
        } else {
            if (.oral){
                .lines[length(.lines) + 1]  <- sprintf("rx_F2 %s 1", ifelse(.derived, "=", "~"))
            } else {
                .lines[length(.lines) + 1]  <- sprintf("rx_F %s 1", ifelse(.derived, "=", "~"))
                .lines[length(.lines) + 1]  <- sprintf("rx_F2 %s 1", ifelse(.derived, "=", "~"))
            }
        }
        .ncmt <- 1;
        if (any(.varsUp == "CL")){
            .hasVp  <- any(regexpr("^(VP[1-9]*)$", .varsUp)!=-1);
            .hasVt  <- any(regexpr("^(VT[1-9]*)$", .varsUp)!=-1);
            .hasVn  <- any(regexpr("^V[1-9]+$", .varsUp)!=-1);
            if (.hasVp && .hasVn){
                stop("cannot mix 'Vp' and 'V#'");
            }
            if (.hasVt && .hasVn){
                stop("cannot mix 'Vt' and 'V#'");
            }
            if (.hasVp && .hasVt){
                stop("cannot mix 'Vp' and 'Vt'");
            }
            if (any(.varsUp=="VP2") && !any(.varsUp=="VP")){
                stop("defined 'VP2' without 'VP'");
            }
            if (any(.varsUp=="VT2") && !any(.varsUp=="VT")){
                stop("defined 'VT2' without 'VT'");
            }
            if (any(.varsUp=="CLD2") && !any(.varsUp=="CLD")){
                stop("defined 'CLD2' without 'CLD'");
            }
            .hasQ  <- any(regexpr("^(Q[1-9]*)$", .varsUp)!=-1);
            .hasCld  <- any(regexpr("^(CLD[1-9]*)$", .varsUp)!=-1);
            if (.hasQ && .hasCld){
                stop("cannot mix 'Q' and 'Cld'")
            }
            .cl <- .getVar("CL");
            .v <- .getVar(c("V", "VC", .vs[1]));
            if (toupper(.v)=="V" && any(.varsUp=="VC")){
                stop(sprintf("ambiguous '%s'/'%s' specification", .v, .getVar("VC")))
            }
            if (toupper(.v) !=.vs[1]){
                .vs <- c(.v, .vs);
            }
            .trans <- 1;
            .lines[length(.lines) + 1] <- sprintf("rx_v1 %s %s", ifelse(.derived, "=", "~"), .v)
            .lines[length(.lines) + 1] <- sprintf("rx_p1 %s %s", ifelse(.derived, "=", "~"), .cl)
            .hasVss  <- FALSE
            if (any(.varsUp == "VSS")){
                .ncmt <- 2;
                .trans <- 3;
                .Q <- .getVar(c("CLD", .qs[1]));
                .vss <- .getVar("VSS");
                .hasVss <- TRUE
                if (.hasVp){
                    stop("cannot define 'Vp' with a 'VSS' defined");
                }
                if (.hasVt){
                    stop("cannot define 'Vt' with a 'VSS' defined");
                }
                .tmp <- .varsUp[.varsUp != toupper(.v)];
                if (any(regexpr("^V[0-9]*$", .tmp)!=-1)){
                    stop("cannot have volumes other than 'Vss' and central volume defined");
                }
                .lines[length(.lines) + 1] <- sprintf("rx_p2 %s %s", ifelse(.derived, "=", "~"), .Q)
                .lines[length(.lines) + 1] <- sprintf("rx_p3 %s %s", ifelse(.derived, "=", "~"), .vss)
            } else if ((any(.varsUp == .vs[2]) || any(.varsUp==.vs[3]) ||
                        any(.varsUp == "VP"))){
                .ncmt <- 2
                .trans <- 1
                .Q <- .getVar(c("CLD",.qs[1]))
                .v2 <- .getVar(c("VP", .vs[2], .vs[3]))
                if (toupper(.v2)==.vs[3]){
                    .vs  <- .vs[-2]
                }
                if (toupper(.v2) !=.vs[2]){
                    .vs  <- c(.vs[1], .v2, .vs[-1])
                }
                if (toupper(.Q) !=.qs[1]){
                    .qs  <- c(.Q,.qs)
                }
                .lines[length(.lines) + 1] <- sprintf("rx_p2 %s %s", ifelse(.derived, "=", "~"), .Q)
                .lines[length(.lines) + 1] <- sprintf("rx_p3 %s %s", ifelse(.derived, "=", "~"), .v2)
            } else if (any(.varsUp == "VT")){
                .ncmt <- 2;
                .trans <- 1
                .Q <- .getVar(c("CLD", .qs[1]));
                if (toupper(.Q)=="CLD"){
                    .qs  <- c("CLD",.qs);
                }
                .v2 <- .getVar("VT");
                if (toupper(.v2) !=.vs[2]){
                    .vs  <- c(.vs[1], .v2, .vs[-1]);
                }
                .lines[length(.lines) + 1] <- sprintf("rx_p2 %s %s", ifelse(.derived, "=", "~"), .Q)
                .lines[length(.lines) + 1] <- sprintf("rx_p3 %s %s", ifelse(.derived, "=", "~"), .v2)
            } else if (any(.varsUp==.qs[1])){
                stop(sprintf("defined '%s' without corresponding volume", .qs[1]))
            } else if (any(.varsUp=="CLD")){
                stop("defined 'CLD' without corresponding volume")
            }
            if (any(.varsUp == toupper(.vs[3]))){
                if (.hasVss){
                    stop("'Vss' only supported with 2 compartment models")
                }
                .ncmt <- 3
                .trans <- 1
                .v3 <- .getVar(.vs[3])
                .q2 <- .getVar(c("CLD2", .qs[2]))
                .lines[length(.lines) + 1] <- sprintf("rx_p4 %s %s", ifelse(.derived, "=", "~"), .q2)
                .lines[length(.lines) + 1] <- sprintf("rx_p5 %s %s", ifelse(.derived, "=", "~"), .v3)
            } else if (any(.varsUp=="VP2")){
                if (.hasVss){
                    stop("'Vss' only supported with 2 compartment models")
                }
                .ncmt <- 3
                .trans <- 1
                .v3 <- .getVar("VP2");
                .q2 <- .getVar(c(.qs[2],"CLD2"))
                .lines[length(.lines) + 1] <- sprintf("rx_p4 %s %s", ifelse(.derived, "=", "~"), .q2)
                .lines[length(.lines) + 1] <- sprintf("rx_p5 %s %s", ifelse(.derived, "=", "~"), .v3)
            } else if (any(.varsUp == "VT2")) {
                if (.hasVss){
                    stop("'Vss' only supported with 2 compartment models")
                }
                .ncmt <- 3;
                .trans <- 1
                .v3 <- .getVar("VT2");
                .q2 <- .getVar(c(.qs[2], "CLD2"));
                .lines[length(.lines) + 1] <- sprintf("rx_p4 %s %s", ifelse(.derived, "=", "~"), .q2);
                .lines[length(.lines) + 1] <- sprintf("rx_p5 %s %s", ifelse(.derived, "=", "~"), .v3);
            } else if (any(.varsUp==.qs[2])){
                if (.hasVss){
                    stop("'Vss' only supported with 2 compartment models")
                }
                stop(sprintf("Defined '%s' without corresponding volume", .qs[2]))
            } else if (any(.varsUp=="CLD2")){
                stop("defined 'CLD2' without corresponding volume")
            }
        } else if (any(.varsUp == "K") || any(.varsUp == "KE") || any(.varsUp == "KEL")) {
            .k <- .getVar(c("K", "KE", "KEL"))
            .v <- .getVar(c("V", "VC", .vs[1]))
            if (toupper(.v)=="V" && any(.varsUp=="VC")){
                stop(sprintf("ambiguous '%s'/'%s' specification", .v, .getVar("VC")))
            }
            if (toupper(.v)=="VC" && any(.varsUp==.vs[1])){
                stop(sprintf("ambiguous '%s'/'%s' specification", .v, .getVar(.vs[1])))
            }
            if (toupper(.v) !=.vs[1]){
                .vs <- c(.v, .vs)
            }
            .lines[length(.lines) + 1] <- sprintf("rx_v1 %s %s", ifelse(.derived, "=", "~"), .v)
            .lines[length(.lines) + 1] <- sprintf("rx_p1 %s %s", ifelse(.derived, "=", "~"), .k)
            .trans <- 2
            if (any(.varsUp == "K12") || any(.varsUp == "K21")){
                .ncmt <- 2
                .trans <- 2
                .k12 <- .getVar("K12")
                .k21 <- .getVar("K21")
                .lines[length(.lines) + 1] <- sprintf("rx_p2 %s %s", ifelse(.derived, "=", "~"), .k12)
                .lines[length(.lines) + 1] <- sprintf("rx_p3 %s %s", ifelse(.derived, "=", "~"), .k21)
            }
            if (any(.varsUp == "K13") || any(.varsUp=="K31")){
                if (.ncmt !=2){
                    stop("'K12' and 'K21' need to be defined for a 3 compartment model")
                }
                .ncmt <- 3
                .k13 <- .getVar("K13")
                .k31 <- .getVar("K31")
                .trans <- 2
                .lines[length(.lines) + 1] <- sprintf("rx_p4 %s %s", ifelse(.derived, "=", "~"), .k13)
                .lines[length(.lines) + 1] <- sprintf("rx_p5 %s %s", ifelse(.derived, "=", "~"), .k31)
            }
        } else if (any(.varsUp == "AOB")){
            .ncmt <- 2;
            .trans <- 5
            .v <- .getVar(c("V", "VC", .vs[1]));
            if (toupper(.v)=="V" && any(.varsUp=="VC")){
                stop(sprintf("ambiguous '%s'/'%s' specification", .v, .getVar("VC")))
            }
            if (toupper(.v)=="VC" && any(.varsUp==.vs[1])){
                stop(sprintf("ambiguous '%s'/'%s' specification", .v, .getVar(.vs[1])))
            }
            if (toupper(.v) !=.vs[1]){
                .vs <- c(.v, .vs)
            }
            .aob <- .getVar("AOB")
            .alpha <- .getVar("ALPHA")
            .beta <- .getVar("BETA")
            .lines[length(.lines) + 1] <- sprintf("rx_v1 %s %s", ifelse(.derived, "=", "~"), .v)
            .lines[length(.lines) + 1] <- sprintf("rx_p1 %s %s", ifelse(.derived, "=", "~"), .alpha)
            .lines[length(.lines) + 1] <- sprintf("rx_p2 %s %s", ifelse(.derived, "=", "~"), .beta)
            .lines[length(.lines) + 1] <- sprintf("rx_p3 %s %s", ifelse(.derived, "=", "~"), .aob)
        } else if (any(.varsUp == "ALPHA") && any(.varsUp == "BETA") && any(.varsUp == "K21")){
            .ncmt <- 2
            .trans <- 4
            .k21 <- .getVar("K21")
            .alpha <- .getVar("ALPHA")
            .beta <- .getVar("BETA")
            .v <- .getVar(c("V", "VC", .vs[1]))
            if (toupper(.v)=="V" && any(.varsUp=="VC")){
                stop(sprintf("ambiguous '%s'/'%s' specification", .v, .getVar("VC")))
            }
            if (toupper(.v)=="VC" && any(.varsUp==.vs[1])){
                stop(sprintf("ambiguous '%s'/'%s' specification", .v, .getVar(.vs[1])))
            }
            if (toupper(.v) !=.vs[1]){
                .vs <- c(.v, .vs)
            }
            .lines[length(.lines) + 1] <- sprintf("rx_v1 %s %s", ifelse(.derived, "=", "~"), .v)
            .lines[length(.lines) + 1] <- sprintf("rx_p3 %s %s", ifelse(.derived, "=", "~"), .k21)
            .lines[length(.lines) + 1] <- sprintf("rx_p1 %s %s", ifelse(.derived, "=", "~"), .alpha)
            .lines[length(.lines) + 1] <- sprintf("rx_p2 %s %s", ifelse(.derived, "=", "~"), .beta)
        } else if (any(.varsUp == "ALPHA") && any(.varsUp == "A")) {
            .trans <- 10
            .isDirect  <- TRUE
        } else if (any(.varsUp=="ALPHA")){
            .ncmt <- 1
            .trans <- 11
            .v <- .getVar(c("V", "VC", .vs[1]))
            .alpha <- .getVar("ALPHA")
            if (toupper(.v)=="V" && any(.varsUp=="VC")){
                stop(sprintf("ambiguous '%s'/'%s' specification", .v, .getVar("VC")))
            }
            if (toupper(.v)=="VC" && any(.varsUp==.vs[1])){
                stop(sprintf("ambiguous '%s'/'%s' specification", .v, .getVar(.vs[1])))
            }
            if (any(.varsUp=="GAMMA")){
                stop("A 'gamma' parameter requires 'A'/'B'/'C' and 'alpha'/'beta'")
            }
            if (any(.varsUp=="K21")){
                stop("'K21' requires a beta parameter")
            }
            if (any(.varsUp=="AOB")){
                stop("'AOB' requires a beta parameter")
            }
            if (any(.varsUp=="BETA")){
                stop("'Beta' requires 'AOB' or 'K21' parameter");
            }
            .lines[length(.lines) + 1] <- sprintf("rx_v1 %s %s", ifelse(.derived, "=", "~"), .v)
            .lines[length(.lines) + 1] <- sprintf("rx_p1 %s %s", ifelse(.derived, "=", "~"), .alpha)
        } else {
            stop("could not figure out the 'linCmt' from the defined parameters")
        }
        if (.isDirect){
            .alpha <- .getVar("ALPHA")
            .a  <- .getVar("A")
            if (any(.varsUp=="BETA") || any(.varsUp=="B")){
              .ncmt <- 2
              .beta <- .getVar("BETA")
              .b <- .getVar("B")
              if (any(.varsUp=="GAMMA") || any(.varsUp=="C")){
                .ncmt <- 3
                ## 3 cmt
                .gamma <- .getVar("GAMMA")
                .c <- .getVar("C")
                .lines[length(.lines) + 1] <- sprintf("rx_p1 %s %s;", ifelse(.derived, "=", "~"), .alpha)
                .lines[length(.lines) + 1] <- sprintf("rx_p2 %s %s;", ifelse(.derived, "=", "~"), .beta)
                .lines[length(.lines) + 1] <- sprintf("rx_p4 %s %s;", ifelse(.derived, "=", "~"), .gamma)
                .lines[length(.lines) + 1] <- sprintf("rx_v1 %s %s;", ifelse(.derived, "=", "~"), .a)
                .lines[length(.lines) + 1] <- sprintf("rx_p3 %s %s;", ifelse(.derived, "=", "~"), .b)
                .lines[length(.lines) + 1] <- sprintf("rx_p5 %s %s;", ifelse(.derived, "=", "~"), .c)
              } else {
                ## 2 cmt
                .lines[length(.lines) + 1] <- sprintf("rx_p1 %s %s", ifelse(.derived, "=", "~"), .alpha)
                .lines[length(.lines) + 1] <- sprintf("rx_p2 %s %s", ifelse(.derived, "=", "~"), .beta)
                .lines[length(.lines) + 1] <-  sprintf("rx_v1 %s %s;", ifelse(.derived, "=", "~"), .a)
                .lines[length(.lines) + 1] <- sprintf("rx_p3 %s %s;", ifelse(.derived, "=", "~"), .b)
                .lines[length(.lines) + 1] <- sprintf("rx_p4 %s 0", ifelse(.derived, "=", "~"))
                .lines[length(.lines) + 1] <- sprintf("rx_p5 %s 0", ifelse(.derived, "=", "~"))
              }
            } else {
                if (any(.varsUp=="GAMMA") || any(.varsUp=="C")){
                    stop("a three compartment model requires 'BETA'/'B'")
                }
                ## 1 cmt
                .lines[length(.lines) + 1] <- sprintf("rx_p1 %s %s", ifelse(.derived, "=", "~"), .alpha)
                .lines[length(.lines) + 1] <- sprintf("rx_v1 %s %s", ifelse(.derived, "=", "~"), .a)
                .lines[length(.lines) + 1] <- sprintf("rx_p2 %s 0", ifelse(.derived, "=", "~"))
                .lines[length(.lines) + 1] <- sprintf("rx_p3 %s 0", ifelse(.derived, "=", "~"))
                .lines[length(.lines) + 1] <- sprintf("rx_p4 %s 0", ifelse(.derived, "=", "~"))
                .lines[length(.lines) + 1] <- sprintf("rx_p5 %s 0", ifelse(.derived, "=", "~"))
            }
        } else {
            if (.ncmt == 1){
                .lines[length(.lines) + 1] <- sprintf("rx_p2 %s 0", ifelse(.derived, "=", "~"))
                .lines[length(.lines) + 1] <- sprintf("rx_p3 %s 0", ifelse(.derived, "=", "~"))
                .lines[length(.lines) + 1] <- sprintf("rx_p4 %s 0", ifelse(.derived, "=", "~"))
                .lines[length(.lines) + 1] <- sprintf("rx_p5 %s 0", ifelse(.derived, "=", "~"))
            } else if (.ncmt == 2){
                .lines[length(.lines) + 1] <- sprintf("rx_p4 %s 0", ifelse(.derived, "=", "~"))
                .lines[length(.lines) + 1] <- sprintf("rx_p5 %s 0", ifelse(.derived, "=", "~"))
            }
        }
        if (.derived){
          return(sprintf(paste0(".Call(`_calcDerived`,as.integer(%d),as.integer(%d),with(.lst,{", paste(.lines, collapse=";\n"),
                        ";\nlist(as.double(rx_p1), as.double(rx_v1), as.double(rx_p2), as.double(rx_p3), as.double(rx_p4), as.double(rx_p5), as.double(%s), as.double(rx_tlag), as.double(rx_tlag2), as.double(rx_F), as.double(rx_F2), as.double(rx_rate), as.double(rx_dur), as.double(rx_rate2), as.double(rx_dur2))}))"), .trans, .ncmt, ifelse(.oral, "rx_ka", "0.0")))
        }
        if (linCmtSens){
            .solve <- sprintf("linCmtB(rx__PTR__, t, %s, %s, %s, 0, rx_p1, rx_v1, rx_p2, rx_p3, rx_p4, rx_p5, %s, rx_tlag, rx_tlag2, rx_F, rx_F2, rx_rate, rx_dur, rx_rate2, rx_dur2)", .linCmt, .ncmt, .trans, ifelse(.oral, "rx_ka", "0.0"));
        } else {
            .solve <- sprintf("linCmtA(rx__PTR__, t, %s, %s, %s, rx_p1, rx_v1, rx_p2, rx_p3, rx_p4, rx_p5, %s, rx_tlag, rx_tlag2, rx_F, rx_F2, rx_rate, rx_dur, rx_rate2, rx_dur2)", .linCmt, .ncmt, .trans, ifelse(.oral, "rx_ka", "0.0"));
        }
        .lines <- paste(.lines, collapse="\n");
        .txt <- paste(sub(.re, sprintf("%s\n\\1%s\\3", .lines, .solve), .txt), collapse="\n");
        ## Put in extra compartment information
        return(rxGetModel(.txt))
    } else {
        stop("can only have one 'linCmt' function in the model");
    }
}

.rxDerivedReg <- rex::rex(start,
                          or(group(or("V", "Q", "VP", "VT", "CLD"), number),
                             "KA", "VP", "VT", "CLD", "V", "VC", "CL", "VSS", "K", "KE", "KEL",
                             group("K", number, number), "AOB", "ALPHA", "BETA", "GAMMA",
                             "A", "B", "C"),
                          end)

##' Calculate derived parameters
##'
##' This calculates the derived parameters based on what is provided
##'
##' @param ...
##'
##' @export
rxDerived <- function(...) {
  .lst <- list(...)
  if (inherits(.lst[[1]], "data.frame")) {
    .lst <- .lst[[1]]
  }
  .namesU <- toupper(names(.lst))
  .w <- which(regexpr(.rxDerivedReg, .namesU) != -1)
  if (length(.w) > 1L){
    .linCmt <- rxLinCmtTrans(paste0("C2=linCmt(", paste(names(.lst)[.w], collapse=","), ")"), NA)
    .env <- environment()
    return(eval(parse(text=.linCmt), envir=.env))
  } else {
    stop("cannot figure out PK parameters to convert")
  }
}
