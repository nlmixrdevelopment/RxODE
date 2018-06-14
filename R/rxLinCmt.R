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
              identical(x[[1]], quote(`=`))) &&
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
    .txt <- try({strsplit(RxODE::rxNorm(modText), "\n")[[1]]}, silent = TRUE);
    if (inherits(.txt, "try-error")){
        return(modText);
    }
    .re <- rex::rex(capture(anything), boundary, "linCmt(", capture(except_any_of(")")), ")", capture(anything));
    .w <- which(regexpr(.re, .txt, perl=TRUE) != -1)
    if (length(.w) == 0){
        return(modText);
    } else if (length(.w) == 1){
        .linCmt <- gsub(.re, "\\2", .txt[.w]);
        if (.linCmt == ""){
            .linCmt <- length(RxODE::rxState(modText));
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
                .linCmt <- length(RxODE::rxState(modText));
            } else {
                .linCmt <- .linCmt - 1;
            }
        }
        .vars <- unique(c(.vars, RxODE::rxLhs(modText), RxODE::rxParam(modText)))
        .varsUp <- toupper(.vars);
        .reg <- rex::rex(start, "V", capture(number), end);
        .w <- which(regexpr(.reg, .varsUp) != -1);
        if (any(.varsUp == "V")){
            if (length(.w) > 0){
                .minV <- min(as.numeric(gsub(.reg, "\\1", .varsUp[.w])));
                .vs <- c("V", paste0("V", seq(.minV, .minV + 1)))
            } else {
                .vs <- c("V", paste0("V", 1:2));
            }
        } else{
            if (length(.w) > 0){
                .minV <- min(as.numeric(gsub(.reg, "\\1", .varsUp[.w])));
                .vs <- paste0("V", seq(.minV, .minV + 2))
            } else {
                .vs <- paste0("V", 1:3);
            }
        }
        .reg <- rex::rex(start, "Q", capture(number), end);
        .w <- which(regexpr(.reg, .varsUp) != -1);
        if (any(.varsUp == "Q")){
            if (length(.w) > 0){
                .minQ <- min(as.numeric(gsub(.reg, "\\1", .varsUp[.w])));
                .qs <- c("Q", paste0("Q", seq(.minQ, .minQ + 1)))
            } else {
                .qs <- c("Q", paste0("Q", 1:2));
            }
        } else{
            if (length(.w) > 0){
                .minQ <- min(as.numeric(gsub(.reg, "\\1", .varsUp[.w])));
                .vs <- paste0("Q", seq(.minQ, .minQ + 2))
            } else {
                .vs <- paste0("Q", 1:3);
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
            if (any(.varsUp == "TLAG")){
                tlag <- .getVar("TLAG");
                .lines[length(.lines) + 1] <- sprintf("rx_tlag ~ %s", tlag);
            } else {
                .lines[length(.lines) + 1] <- "rx_tlag ~ 0"
            }
        } else {
            .lines[length(.lines) + 1] <- sprintf("rx_ka ~ 0");
            .lines[length(.lines) + 1] <- sprintf("rx_tlag ~ 0");
        }
        .ncmt <- 1;
        if (any(.varsUp == "CL")){
            .cl <- .getVar("CL");
            .v <- .getVar(c("V", "VC", .vs[1]));
            .lines[length(.lines) + 1] <- sprintf("rx_v ~ %s", .v);
            .lines[length(.lines) + 1] <- sprintf("rx_k ~ %s/%s", .cl, .v);
            if ((any(.varsUp == .vs[2]) || any(.varsUp == "VP"))){
                .ncmt <- 2;
                .Q <- .getVar(.qs[1]);
                .v2 <- .getVar(c(.vs[2], "VP"));
                .lines[length(.lines) + 1] <- sprintf("rx_k12 ~ %s/%s", .Q, .v);
                .lines[length(.lines) + 1] <- sprintf("rx_k21 ~ %s/%s", .Q, .v2);
            } else if (any(.varsUp == "VT")){
                .ncmt <- 2;
                .Q <- .getVar("CLD");
                .v2 <- .getVar("VT");
                .lines[length(.lines) + 1] <- sprintf("rx_k12 ~ %s/%s", .Q, .v);
                .lines[length(.lines) + 1] <- sprintf("rx_k21 ~ %s/%s", .Q, .v2);
            } else if (any(.varsUp == "VSS")){
                .ncmt <- 2;
                .Q <- .getVar(.qs[1]);
                .vss <- .getVar("VSS");
                .lines[length(.lines) + 1] <- sprintf("rx_k12 ~ %s/%s", .Q, .v);
                .lines[length(.lines) + 1] <- sprintf("rx_k21 ~ %s/(%s-%s)", .Q, .vss, .v);
            }
            if (any(.varsUp == .vs[3])){
                .ncmt <- 3;
                .v3 <- .getVar(.vs[3]);
                .q2 <- .getVar(.qs[2]);
                .lines[length(.lines) + 1] <- sprintf("rx_k13 ~ %s/%s", .q2, .v);
                .lines[length(.lines) + 1] <- sprintf("rx_k31 ~ %s/%s", .q2, .v3);
            } else if (any(.varsUp == "VT2")) {
                .ncmt <- 3;
                .v3 <- .getVar("VT2");
                .q2 <- .getVar("CLD2");
                .lines[length(.lines) + 1] <- sprintf("rx_k13 ~ %s/%s", .q2, .v);
                .lines[length(.lines) + 1] <- sprintf("rx_k31 ~ %s/%s", .q2, .v3);
            }
        } else if (any(.varsUp == "K") || any(.varsUp == "KE") || any(.varsUp == "KEL")) {
            .k <- .getVar(c("K", "KE", "KEL"))
            .v <- .getVar(c("V", "VC", .vs[1]));
            .lines[length(.lines) + 1] <- sprintf("rx_v ~ %s", .v);
            .lines[length(.lines) + 1] <- sprintf("rx_k ~ %s", .k);
            if (any(.varsUp == "K12")){
                .ncmt <- 2;
                .k12 <- .getVar("K12")
                .k12 <- .getVar("K21")
                .lines[length(.lines) + 1] <- sprintf("rx_k12 ~ %s", .k12);
                .lines[length(.lines) + 1] <- sprintf("rx_k21 ~ %s", .k21);
            }
            if (any(.varsUp == "K13")){
                .ncmt <- 3;
                .k13 <- .getVar("K13")
                .k31 <- .getVar("K31")
                .lines[length(.lines) + 1] <- sprintf("rx_k13 ~ %s", .k13);
                .lines[length(.lines) + 1] <- sprintf("rx_k31 ~ %s", .k31);
            }
        } else if (any(.varsUp == "AOB")){
            .ncmt <- 2;
            .v <- .getVar(c("V", "VC", .vs[1]));
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
            .lines[length(.lines) + 1] <- sprintf("rx_v ~ %s", .v);
            .lines[length(.lines) + 1] <- sprintf("rx_k21 ~ %s", .k21);
            .lines[length(.lines) + 1] <- sprintf("rx_k ~ (%s*%s)/rx_k21", .alpha, .beta);
            .lines[length(.lines) + 1] <- sprintf("rx_k12 ~ %s+%s - rx_k21 - rx_k", .alpha, .beta);
        }
        if (.ncmt == 1){
            .lines[length(.lines) + 1] <- "rx_alpha ~ rx_k";
            if (.oral){
                .lines[length(.lines) + 1] <- "rx_A ~ rx_ka / (rx_ka - rx_alpha) / rx_v";
            } else {
                .lines[length(.lines) + 1] <- "rx_A ~ 1.0 / rx_v";
            }
            .lines[length(.lines) + 1] <- "rx_beta ~ 0";
            .lines[length(.lines) + 1] <- "rx_B ~ 0";
            .lines[length(.lines) + 1] <- "rx_gamma ~ 0";
            .lines[length(.lines) + 1] <- "rx_C ~ 0";
        } else if (.ncmt == 2){
            .lines[length(.lines) + 1] <- "rx_beta  ~ 0.5 * (rx_k12 + rx_k21 + rx_k - sqrt((rx_k12 + rx_k21 + rx_k) * (rx_k12 + rx_k21 + rx_k) - 4.0 * rx_k21 * rx_k))"
            .lines[length(.lines) + 1] <-  "rx_alpha ~ rx_k21 * rx_k / rx_beta"
            if (.oral){
                .lines[length(.lines) + 1] <-  "rx_A ~ rx_ka / (rx_ka - rx_alpha) * (rx_alpha - rx_k21) / (rx_alpha - rx_beta) / rx_v"
                .lines[length(.lines) + 1] <- "rx_B ~ rx_ka / (rx_ka - rx_beta) * (rx_beta - rx_k21) / (rx_beta - rx_alpha) / rx_v;"
            } else {
                .lines[length(.lines) + 1] <-  "rx_A ~ (rx_alpha - rx_k21) / (rx_alpha - rx_beta) / rx_v"
                .lines[length(.lines) + 1] <- "rx_B ~ (rx_beta - rx_k21) / (rx_beta - rx_alpha) / rx_v;"
            }
            .lines[length(.lines) + 1] <- "rx_gamma ~ 0";
            .lines[length(.lines) + 1] <- "rx_C ~ 0";
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
                .lines[length(.lines) + 1] <- "rx_A ~ rx_ka / (rx_ka - rx_alpha) * rx_A";
                .lines[length(.lines) + 1] <- "rx_B ~ rx_ka / (rx_ka - rx_beta) * rx_B";
                .lines[length(.lines) + 1] <- "rx_C ~ rx_ka / (rx_ka - rx_gamma) * rx_C";
            }
        }
        .solve <- sprintf("solveLinB(rx__PTR__, t, %s, 0, 0, rx_A, rx_alpha, rx_B, rx_beta, rx_C, rx_gamma, rx_ka, rx_tlag)", .linCmt);
        .lines <- paste(.lines, collapse="\n");
        .txt <- paste(sub(.re, sprintf("%s\n\\1%s\\3", .lines, .solve), .txt), collapse="\n");
        return(.txt)
    } else {
        stop("Can only have one linCmt() function in the model.  Assign it to a variable if you need the concentrations more than once.");
    }
}
