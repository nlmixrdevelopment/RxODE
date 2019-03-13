
require(RxODE);

rxPermissive({

    context("etTrans checks");

    mod <- RxODE("
a = 6
b = 0.6
d/dt(intestine) = -a*intestine
d/dt(blood)     = a*intestine - b*blood
")

    et <- eventTable()
    et$add.dosing(dose=2/24,rate=2,start.time=0,
                  nbr.doses=10,dosing.interval=1)
    et <- et %>% et(0.05,evid=2) %>%
        et(amt=3,time=0.5,cmt=out) %>%
        et(amt=3,time=0.1,cmt=intestine,ss=1,ii=3) %>%
        et(amt=3,time=0.3,cmt=intestine,ss=2,ii=3) %>%
        et(time=0.2,cmt="-intestine") %>%
        as.data.frame

    ett1 <- RxODE:::etTrans(et, mod)
    tmp1 <- sort(unique(ett1$evid));

    et$cmt <- factor(et$cmt)
    ett2 <- RxODE:::etTrans(et, mod);

    tmp2 <- sort(unique(ett2$evid));

    test_that("factor and character give same compartment information",{
        expect_equal(attr(class(ett2), ".RxODE.lst")$cmtInfo, attr(class(ett1), ".RxODE.lst")$cmtInfo);
        expect_equal(attr(class(ett2), ".RxODE.lst")$cmtInfo, c("intestine", "blood", "out"))
    })

    test_that("factor and character give same evids",{
        expect_equal(tmp1,tmp2);
        expect_equal(tmp1, c(2L, 101L, 110L, 120L, 130L, 301L, 10101L))
    })


    et <- eventTable()
    et$add.dosing(dose=2/24,rate=2,start.time=0,
                  nbr.doses=10,dosing.interval=1)
    et <- et %>% et(0.05,evid=2) %>%
        et(amt=3,time=0.5,cmt="-out") %>% as.data.frame

    test_that("error for negative non ODE compartments",{
        expect_error(RxODE:::etTrans(et, mod))
        et$cmt <- factor(et$cmt)
        expect_error(RxODE:::etTrans(et, mod))
    })

    et <- eventTable()
    et$add.dosing(dose=2/24,rate=2,start.time=0,
                  nbr.doses=10,dosing.interval=1)
    et <- et %>% et(0.05,evid=2) %>%
        et(amt=3,time=0.25,cmt="out") %>%
        et(amt=3,time=0.5,cmt="-out") %>% as.data.frame

    test_that("error for negative non ODE compartments after defined compartment", {
        expect_error(RxODE:::etTrans(et, mod))
        et$cmt <- factor(et$cmt)
        expect_error(RxODE:::etTrans(et, mod))
    })

    et <- et() %>% et(amt=3,time=0.24,evid=4)

    test_that("EVID=4 makes sense", {
        expect_equal(RxODE:::etTrans(et, mod)$evid, c(3L, 101L))
    })

}, cran=TRUE, silent=TRUE)
