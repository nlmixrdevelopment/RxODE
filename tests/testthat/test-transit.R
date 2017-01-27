context("Test Transit compartment model");
library(digest);
rxPermissive({

    mod <- RxODE("
## Table 3 from Savic 2007
cl = 17.2 # (L/hr)
vc = 45.1 # L
ka = 0.38 # 1/hr
mtt = 0.37 # hr
bio=1
n = 20.1
k = cl/vc
ktr = (n+1)/mtt
## note that lgammafn is the same as lgamma in R.
d/dt(depot) = exp(log(bio*podo)+log(ktr)+n*log(ktr*t)-ktr*t-lgammafn(n+1))-ka*depot
d/dt(cen) = ka*depot-k*cen
")

    et <- eventTable();
    et$add.sampling(seq(0, 10, length.out=200));
    et$add.dosing(20, start.time=0);

    transit <- rxSolve(mod, et, transit_abs=TRUE)

    test_that("Transit absorption works", {
        expect_equal(digest(round(as.data.frame(transit), 4)),
                     "19470df8fd8bcac4f87fd2ee6431791a");
    })

    transit <- rxSolve(mod, et)

    test_that("Transit absorption is not specified, but still works.", {
        expect_equal(digest(round(as.data.frame(transit), 4)),
                     "19470df8fd8bcac4f87fd2ee6431791a");
        expect_warning(rxSolve(mod, et));
    })

    no_transit <- rxSolve(mod, et, transit_abs=FALSE);

    test_that("Transit absorption is turned off, and gives other results", {
        expect_equal(digest(round(as.data.frame(no_transit), 4)),
                     "0c8ba19d3e8fa307a1f7256afd556f40");
    })

    mod <- RxODE("
## Table 3 from Savic 2007
cl = 17.2 # (L/hr)
vc = 45.1 # L
ka = 0.38 # 1/hr
mtt = 0.37 # hr
bio=1
n = 20.1
k = cl/vc
ktr = (n+1)/mtt
## note that lgamma1p is the same as lgamma(1+p) in R.
d/dt(depot) = exp(log(bio*podo)+log(ktr)+n*log(ktr*t)-ktr*t-lgamma1p(n))-ka*depot
d/dt(cen) = ka*depot-k*cen
")

    transit <- rxSolve(mod, et);

    test_that("Transit absorption using lagmma1p works.", {
        expect_equal(digest(round(as.data.frame(transit), 4)),
                     "19470df8fd8bcac4f87fd2ee6431791a");
        expect_warning(rxSolve(mod, et));
    })

    mod <- RxODE("
## Table 3 from Savic 2007
cl = 17.2 # (L/hr)
vc = 45.1 # L
ka = 0.38 # 1/hr
mtt = 0.37 # hr
bio=1
n = 20.1
k = cl/vc
ktr = (n+1)/mtt
d/dt(depot) = exp(log(bio*podo)+log(ktr)+n*log(ktr*t)-ktr*t-log(n!))-ka*depot
d/dt(cen) = ka*depot-k*cen
")

    transit <- rxSolve(mod, et);

    test_that("Transit absorption using !.", {
        expect_equal(digest(round(as.data.frame(transit), 4)),
                     "19470df8fd8bcac4f87fd2ee6431791a");
        expect_warning(rxSolve(mod, et));
    })

    mod <- RxODE("
## Table 3 from Savic 2007
cl = 17.2 # (L/hr)
vc = 45.1 # L
ka = 0.38 # 1/hr
mtt = 0.37 # hr
bio=1
n = 20.1
k = cl/vc
ktr = (n+1)/mtt
d/dt(depot) = exp(log(bio*podo)+log(ktr)+n*log(ktr*t)-ktr*t-lgamma(n+1))-ka*depot
d/dt(cen) = ka*depot-k*cen
")

    transit <- rxSolve(mod, et);

    test_that("Transit absorption using C's lgamma function.", {
        expect_equal(digest(round(as.data.frame(transit), 4)),
                     "19470df8fd8bcac4f87fd2ee6431791a");
        expect_warning(rxSolve(mod, et));
    })

    mod <- RxODE("
## Table 3 from Savic 2007
cl = 17.2 # (L/hr)
vc = 45.1 # L
ka = 0.38 # 1/hr
mtt = 0.37 # hr
bio=1
n = 20.1
k = cl/vc
ktr = (n+1)/mtt
d/dt(depot) = exp(log(bio*podo)+log(ktr)+n*log(ktr*t)-ktr*t-lfactorial(n))-ka*depot
d/dt(cen) = ka*depot-k*cen
")

    transit <- rxSolve(mod, et);

    test_that("Transit absorption using lfactorial", {
        expect_equal(digest(round(as.data.frame(transit), 4)),
                     "19470df8fd8bcac4f87fd2ee6431791a");
        expect_warning(rxSolve(mod, et));
    })

    mod <- RxODE("
## Table 3 from Savic 2007
cl = 17.2 # (L/hr)
vc = 45.1 # L
ka = 0.38 # 1/hr
mtt = 0.37 # hr
bio=1
n = 20.1
k = cl/vc
ktr = (n+1)/mtt
d/dt(depot) = exp(log(bio*podo)+log(ktr)+n*log(ktr*t)-ktr*t-log(factorial(n)))-ka*depot
d/dt(cen) = ka*depot-k*cen
")

    transit <- rxSolve(mod, et);

    test_that("Transit absorption using log(factorial)", {
        expect_equal(digest(round(as.data.frame(transit), 4)),
                     "19470df8fd8bcac4f87fd2ee6431791a");
        expect_warning(rxSolve(mod, et));
    })

    mod <- RxODE("
## Table 3 from Savic 2007
cl = 17.2 # (L/hr)
vc = 45.1 # L
ka = 0.38 # 1/hr
mtt = 0.37 # hr
n = 20.1
k = cl/vc
## note that lgammafn is the same as lgamma in R.
d/dt(depot) = transit(n, mtt)-ka*depot
d/dt(cen) = ka*depot-k*cen
")

    transit2 <- rxSolve(mod, et);

    test_that("Transit absorption is function that can take 2 arguments", {
        expect_equal(round(as.data.frame(transit)[,names(transit) != "ktr"], 4),
                     round(as.data.frame(transit2), 4));
    })

    mod <- RxODE("
## Table 3 from Savic 2007
cl = 17.2 # (L/hr)
vc = 45.1 # L
ka = 0.38 # 1/hr
mtt = 0.37 # hr
bio=1
n = 20.1
k = cl/vc
d/dt(depot) = transit(n, mtt, bio)-ka*depot
d/dt(cen) = ka*depot-k*cen
")

    transit2 <- rxSolve(mod, et);

    test_that("Transit absorption is function that can take 3 arguments", {
        expect_equal(round(as.data.frame(transit)[,names(transit) != "ktr"], 4),
                     round(as.data.frame(transit2), 4));
    })

    mod <- RxODE("
## Table 3 from Savic 2007
cl = 17.2 # (L/hr)
vc = 45.1 # L
ka = 0.38 # 1/hr
mtt = 0.37 # hr
bio=1
n = 20.1
k = cl/vc
d/dt(depot) = transit(n, mtt, 1)-ka*depot
d/dt(cen) = ka*depot-k*cen
")

    transit2 <- rxSolve(mod, et);

    test_that("Transit absorption can take numeric arguments", {
        expect_equal(round(as.data.frame(transit)[,names(transit) != "ktr"], 4),
                     round(as.data.frame(transit2), 4));
    })
}, silent=TRUE)
