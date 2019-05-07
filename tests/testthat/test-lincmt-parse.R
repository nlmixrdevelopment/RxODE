rxPermissive({

    tran1  <- expand.grid(Ka=NA,## c("Ka",NA),
                          Vc=c("V","Vc","V1", NA),
                          Cl=c("Cl",NA),
                          Q=c("Q","Q1", "Cld", NA),
                          Vp=c("V2","Vp", "Vt", NA),
                          Q2=c("Q2","Cld2", NA),
                          Vp2=c("Vp2","V3", "Vt2", NA))

    library(RxODE);
    library(testthat)

    .fun  <- function(x){
        assign(".x", x, globalenv())
        x  <- setNames(x,names(tran1))
        .v1  <- as.character(na.omit(c(x["Vc"], x["Vp"], x["Vp2"])))
        .v2  <- as.character(na.omit(c(x["Cl"], x["Q"], x["Q2"])))
        .rx  <- paste(c(ifelse(is.na(x["Ka"]), "", paste0(x["Ka"], "=tKa*exp(eta.ka)")),
                        ifelse(is.na(x["Vc"]), "", paste0(x["Vc"], "=tVc*exp(eta.vc)")),
                        ifelse(is.na(x["Cl"]), "", paste0(x["Cl"], "=tCl*exp(eta.cl)")),
                        ifelse(is.na(x["Q"]), "", paste0(x["Q"], "=tQ*exp(eta.q)")),
                        ifelse(is.na(x["Vp"]), "", paste0(x["Vp"], "=tVp*exp(eta.vp)")),
                        ifelse(is.na(x["Q2"]), "", paste0(x["Q2"], "=tq2*exp(eta.q)")),
                        ifelse(is.na(x["Vp2"]), "", paste0(x["Vp2"], "=tvp2*exp(eta.tvp2)")),
                        "cp=linCmt()"
                        ),collapse="\n");
        assign(".rx", .rx, globalenv());
        .good <- FALSE;
        if (length(.v1) ==length(.v2)){
            .good <- TRUE;
            if (length(.v1)==0){
                test_that(sprintf("linCmt() should error without parameters"),{
                    expect_error(RxODE(.rx))
                })
                .good <- NA;
            } else {
                .varsUp  <- toupper(c(.v1, .v2))
                if (!any(.varsUp=="CL")){
                    .good <- FALSE
                }
                .hasVc <- any(regexpr("^V(C|[1-9]|)$",.varsUp)!=-1)
                if (!.hasVc){
                    .good <- FALSE
                }
                .hasVp  <- any(regexpr("^(VP[1-9]*)$", .varsUp)!=-1);
                .hasVt  <- any(regexpr("^(VT[1-9]*)$", .varsUp)!=-1);
                .hasVn  <- any(regexpr("^V[1-9]+$", .varsUp)!=-1);
                if (.hasVp && .hasVn){
                    .good <- FALSE
                }
                if (.hasVt && .hasVn){
                    .good <- FALSE
                }
                if (.hasVp && .hasVt){
                    .good <- FALSE
                }
                .hasQ  <- any(regexpr("^(Q[1-9]*)$", .varsUp)!=-1);
                .hasCld  <- any(regexpr("^(CLD[1-9]*)$", .varsUp)!=-1);
                if (.hasQ && .hasCld){
                    .good <- FALSE
                }
                if (any("VP2"==.varsUp) && !any("VP"==.varsUp)){
                    .good <- FALSE
                }
                if (any("VT2"==.varsUp) && !any("VT"==.varsUp)){
                    .good <- FALSE
                }
                if (any("CLD2"==.varsUp) && !any("CLD"==.varsUp)){
                    .good <- FALSE
                }
            }
        }
        if (is.na(.good)){
        } else if (.good){
            test_that(sprintf("linCmt() successful with parameters: %s", paste(na.omit(x),collapse=", ")),{
                .rx <- RxODE(.rx)
                expect_true(inherits(.rx, "RxODE"));
            })
        } else {
            test_that(sprintf("linCmt() should error with parameters: %s", paste(na.omit(x),collapse=", ")),{
                expect_error(RxODE(.rx))
            })
        }
    }

    context("Cl style translations")
    apply(tran1, 1, .fun)

}, on.validate=TRUE);
