rxPermissive({
    context("Nesting tests")

    mod <- RxODE({
        eff(0) = 1
        C2 = centr/V2*(1+prop.err);
        C3 = peri/V3;
        CL =  TCl*exp(eta.Cl + iov.Cl)
        KA = TKA * exp(eta.Ka + iov.Ka)
        d/dt(depot) =-KA*depot;
        d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
        d/dt(peri)  =                    Q*C2 - Q*C3;
        d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
    })

    mod.eta <- RxODE({
        eff(0) = 1
        C2 = centr/V2*(1+prop.err);
        C3 = peri/V3;
        CL =  TCl*exp(ETA[1] + iov.Cl)
        KA = TKA * exp(ETA[2] + iov.Ka)
        d/dt(depot) =-KA*depot;
        d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
        d/dt(peri)  =                    Q*C2 - Q*C3;
        d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
    })

    library(dplyr)

    et(amountUnits="mg", timeUnits="hours") %>%
        et(amt=10000, addl=9,ii=12,cmt="depot") %>%
        et(time=120, amt=2000, addl=4, ii=14, cmt="depot") %>%
        et(seq(0, 240, by=4)) %>% # Assumes sampling when there is no dosing information
        et(seq(0, 240, by=4) + 0.1) %>% ## adds 0.1 for separate eye
        et(id=1:20) %>%
        ## Add an occasion per dose
        mutate(occ=cumsum(!is.na(amt))) %>%
        mutate(occ=ifelse(occ == 0, 1, occ)) %>%
        mutate(occ=2- occ %% 2) %>%
        mutate(eye=ifelse(round(time) == time, 1, 2)) %>%
        mutate(inv=ifelse(id < 10, 1, 2)) ->
        ev

    omega <- lotri(lotri(eta.Cl ~ 0.1,
                         eta.Ka ~ 0.1) | id,
                   lotri(eye.Cl ~ 0.05,
                         eye.Ka ~ 0.05) | eye,
                   lotri(iov.Cl ~ 0.01,
                         iov.Ka ~ 0.01) | occ,
                   lotri(inv.Cl ~ 0.02,
                         inv.Ka ~ 0.02) | inv)

    .ni <- .nestingInfo(ev$id, omega, ev)
    expect_equal(.ni$below, c(eye = 2L, occ = 2L))
    expect_equal(.ni$above, c(inv = 2L))
    expect_true(inherits(.ni$data$eye, "factor"))
    expect_equal(attr(.ni$data$eye, "nu"), 40L)
    expect_true(inherits(.ni$data$inv, "factor"))
    expect_equal(attr(.ni$data$inv, "nu"), NULL)

    expect_true(inherits(.ni$data$occ, "factor"))
    expect_equal(attr(.ni$data$occ, "nu"), 40L)

    expect_equal(.ni$extraTheta, 4)
    expect_equal(.ni$extraEta, 8)

    rxExpandNesting(mod, .ni)

    ## Tests -- Test different occ size
    ## Tests -- Eye becomes OD/OS

    ## Test edge case -- no between or above occasion variability

    .ni <- .nestingInfo(ev$id, lotri(eta.Cl ~ 0.1, eta.Ka ~ 0.1),
                        ev)

    expect_equal(.ni$above, NULL)
    expect_equal(.ni$below, NULL)
    expect_equal(.ni$idName, "id")
    expect_true(inherits(.ni$omega, "lotri"))
    expect_equal(names(.ni$omega), "id")


    ## test_that("Expanding IOV", {

##         expect_error(RxODE:::rxExpandOcc(mod,1,"iov.Cl"))

##         expect_equal(RxODE:::rxExpandOcc(mod,2,c("iov.Cl", "iov.Ka")),
##                      list(mod="iov.Ka=(rxOCC==2)*ETA[4]+(rxOCC==1)*ETA[2];\niov.Cl=(rxOCC==2)*ETA[3]+(rxOCC==1)*ETA[1];\neff(0)=1;\nC2=centr/V2*(1+prop.err);\nC3=peri/V3;\nCL=TCl*exp(eta.Cl+iov.Cl);\nKA=TKA*exp(eta.Ka+iov.Ka);\nd/dt(depot)=-KA*depot;\nd/dt(centr)=KA*depot-CL*C2-Q*C2+Q*C3;\nd/dt(peri)=Q*C2-Q*C3;\nd/dt(eff)=Kin-Kout*(1-C2/(EC50+C2))*eff;\n",
##                           iov1=paste0("ETA[", 1:4, "]"),
##                           iov2=c("iov.Cl", "iov.Ka"),
##                           iov3=c("iov.Cl(1)", "iov.Ka(1)",
##                                  "iov.Cl(2)", "iov.Ka(2)")))

##         expect_equal(RxODE:::rxExpandOcc(mod,2,c("iov.Ka", "iov.Cl")),
##                      list(mod="iov.Cl=(rxOCC==2)*ETA[4]+(rxOCC==1)*ETA[2];\niov.Ka=(rxOCC==2)*ETA[3]+(rxOCC==1)*ETA[1];\neff(0)=1;\nC2=centr/V2*(1+prop.err);\nC3=peri/V3;\nCL=TCl*exp(eta.Cl+iov.Cl);\nKA=TKA*exp(eta.Ka+iov.Ka);\nd/dt(depot)=-KA*depot;\nd/dt(centr)=KA*depot-CL*C2-Q*C2+Q*C3;\nd/dt(peri)=Q*C2-Q*C3;\nd/dt(eff)=Kin-Kout*(1-C2/(EC50+C2))*eff;\n",
##                           iov1=paste0("ETA[", 1:4, "]"),
##                           iov2=c("iov.Ka", "iov.Cl"),
##                           iov3=c("iov.Ka(1)", "iov.Cl(1)",
##                                  "iov.Ka(2)", "iov.Cl(2)")))

##         expect_equal(RxODE:::rxExpandOcc(mod.eta,3,c("iov.Cl", "iov.Ka")),
##                      list(mod="iov.Ka=(rxOCC==3)*ETA[8]+(rxOCC==2)*ETA[6]+(rxOCC==1)*ETA[4];\niov.Cl=(rxOCC==3)*ETA[7]+(rxOCC==2)*ETA[5]+(rxOCC==1)*ETA[3];\neff(0)=1;\nC2=centr/V2*(1+prop.err);\nC3=peri/V3;\nCL=TCl*exp(ETA[1]+iov.Cl);\nKA=TKA*exp(ETA[2]+iov.Ka);\nd/dt(depot)=-KA*depot;\nd/dt(centr)=KA*depot-CL*C2-Q*C2+Q*C3;\nd/dt(peri)=Q*C2-Q*C3;\nd/dt(eff)=Kin-Kout*(1-C2/(EC50+C2))*eff;\n",
##                           iov1=paste0("ETA[", 3:8, "]"),
##                           iov2=c("iov.Cl", "iov.Ka"),
##                           iov3=c("iov.Cl(1)", "iov.Ka(1)",
##                                  "iov.Cl(2)", "iov.Ka(2)",
##                                  "iov.Cl(3)", "iov.Ka(3)")))

##         .i <- RxODE:::rxExpandOcc(mod.eta,3,c("iov.Cl", "iov.Ka"), TRUE)
##         expect_true(inherits(.i[[1]], "RxODE"))

##         omega <- lotri(eta.cl+eta.v ~c(1, 0.5, 1))
##         iov <- lotri(iov.cl+iov.v ~c(1, 0.5, 1))

##         .m1 <- RxODE:::.rxExpandOmegaIov(omega, iov, paste0("ETA[",1:4, "]"))

##         expect_equal(.m1,
##                      list(m = structure(c(1, 0.5, 0, 0, 0, 0, 0.5, 1, 0, 0, 0, 0,
## 0, 0, 1, 0.5, 0, 0, 0, 0, 0.5, 1, 0, 0, 0, 0, 0, 0, 1, 0.5, 0,
## 0, 0, 0, 0.5, 1), .Dim = c(6L, 6L), .Dimnames = list(c("eta.cl",
## "eta.v", "ETA[1]", "ETA[2]", "ETA[3]", "ETA[4]"), c("eta.cl",
## "eta.v", "ETA[1]", "ETA[2]", "ETA[3]", "ETA[4]"))), lower = c(-Inf,
## -Inf, -Inf, -Inf, -Inf, -Inf), upper = c(Inf, Inf, Inf, Inf,
## Inf, Inf)))

##         .m1 <- RxODE:::.rxExpandOmegaIov(omega, iov, paste0("ETA[",1:4, "]"),
##                                          omegaLower= -1, omegaUpper=1,
##                                          iovLower= -1.5, iovUpper=1.5)
##     expect_equal(.m1,
##                  list(m = structure(c(1, 0.5, 0, 0, 0, 0, 0.5, 1, 0, 0, 0, 0,
## 0, 0, 1, 0.5, 0, 0, 0, 0, 0.5, 1, 0, 0, 0, 0, 0, 0, 1, 0.5, 0,
## 0, 0, 0, 0.5, 1), .Dim = c(6L, 6L), .Dimnames = list(c("eta.cl",
## "eta.v", "ETA[1]", "ETA[2]", "ETA[3]", "ETA[4]"), c("eta.cl",
## "eta.v", "ETA[1]", "ETA[2]", "ETA[3]", "ETA[4]"))), lower = c(-1,
## -1, -1.5, -1.5, -1.5, -1.5), upper = c(1, 1, 1.5, 1.5, 1.5, 1.5)))

##     .m1 <- RxODE:::.rxExpandOmegaIov(omega, iov, paste0("ETA[",1:4, "]"),
##                                      omegaLower= -c(1, 1.5),
##                                      omegaUpper=c(1, 1.5),
##                                      iovLower= -c(2.5, 1.5),
##                                      iovUpper=c(2.5, 1.5))

##     expect_equal(.m1, list(m = structure(c(1, 0.5, 0, 0, 0, 0, 0.5, 1, 0, 0, 0, 0,
## 0, 0, 1, 0.5, 0, 0, 0, 0, 0.5, 1, 0, 0, 0, 0, 0, 0, 1, 0.5, 0,
## 0, 0, 0, 0.5, 1), .Dim = c(6L, 6L), .Dimnames = list(c("eta.cl",
## "eta.v", "ETA[1]", "ETA[2]", "ETA[3]", "ETA[4]"), c("eta.cl",
## "eta.v", "ETA[1]", "ETA[2]", "ETA[3]", "ETA[4]"))), lower = c(-1,
## -1.5, -2.5, -1.5, -2.5, -1.5), upper = c(1, 1.5, 2.5, 1.5, 2.5,
## 1.5)))


##         expect_error(RxODE:::.rxExpandOmegaIov(omega, iov, paste0("ETA[",1:3, "]"),
##                                                omegaLower= -1, omegaUpper=1,
##                                                iovLower= -1.5, iovUpper=1.5))

##     expect_error(RxODE:::.rxExpandOmegaIov(omega, iov, paste0("ETA[",1:4, "]"),
##                                            omegaLower= -c(1, 1, 1)))

##     expect_error(RxODE:::.rxExpandOmegaIov(omega, iov, paste0("ETA[",1:4, "]"),
##                                            omegaUpper= -c(1, 1, 1)))

##     expect_error(RxODE:::.rxExpandOmegaIov(omega, iov, paste0("ETA[",1:4, "]"),
##                                            iovUpper= -c(1, 1, 1)))

##     expect_error(RxODE:::.rxExpandOmegaIov(omega, iov, paste0("ETA[",1:4, "]"),
##                                      iovLower= -c(1, 1, 1)))


##         dimnames(iov) <- NULL
##         expect_error(RxODE:::.rxExpandOmegaIov(omega, iov, paste0("ETA[",1:4, "]"),
##                                          omegaLower= -1, omegaUpper=1,
##                                          iovLower= -1.5, iovUpper=1.5))

##         RxODE:::.rxExpandOmegaIov(omega, omegaLower= -1, omegaUpper=1)

##         dimnames(omega) <- NULL
##         expect_error(RxODE:::.rxExpandOmegaIov(omega, omegaLower= -1, omegaUpper=1))

##     })

})
