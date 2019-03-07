require(RxODE);
require(digest)
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
        et(amt=3,time=0.5,cmt="out") %>% as.data.frame

    tmp1 <- sort(unique(RxODE:::etTrans(et, mod)$evid));

    et$cmt <- factor(et$cmt)
    tmp2 <- sort(unique(RxODE:::etTrans(et, mod)$evid));

    test_that("factor and character give same evids",{
        expect_equal(tmp1,tmp2);
        expect_equal(tmp1,c(2L, 101L, 301L, 10101L))
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

}, cran=TRUE, silent=TRUE)
