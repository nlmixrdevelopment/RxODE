context("Test Transit compartment model");
library(digest);

dat1 <- structure(list(time = c(0, 0.0503, 0.1005, 0.1508, 0.201, 0.2513,  0.3015, 0.3518, 0.402, 0.4523, 0.5025, 0.5528, 0.603, 0.6533,  0.7035, 8.9447, 8.995, 9.0452, 9.0955, 9.1457, 9.196, 9.2462,  9.2965, 9.3467, 9.397, 9.4472, 9.4975, 9.5477, 9.598, 9.6482,  9.6985), depot = c(0, 0, 0, 0.0044, 0.1348, 1.096, 3.9792, 8.5762,  13.1524, 16.2823, 17.8004, 18.2531, 18.1919, 17.9355, 17.6201,  0.7694, 0.7548, 0.7405, 0.7265, 0.7128, 0.6993, 0.6861, 0.6731,  0.6604, 0.6479, 0.6356, 0.6236, 0.6118, 0.6002, 0.5889, 0.5777 ), cen = c(0, 0, 0, 0, 8e-04, 0.0102, 0.0546, 0.1709, 0.3748,  0.6489, 0.9611, 1.285, 1.6058, 1.9171, 2.217, 2.4915, 2.4586,  2.426, 2.3939, 2.362, 2.3306, 2.2994, 2.2686, 2.2382, 2.2081,  2.1783, 2.1488, 2.1197, 2.091, 2.0625, 2.0344), k = c(0.3814,  0.3814, 0.3814, 0.3814, 0.3814, 0.3814, 0.3814, 0.3814, 0.3814,  0.3814, 0.3814, 0.3814, 0.3814, 0.3814, 0.3814, 0.3814, 0.3814,  0.3814, 0.3814, 0.3814, 0.3814, 0.3814, 0.3814, 0.3814, 0.3814,  0.3814, 0.3814, 0.3814, 0.3814, 0.3814, 0.3814), ktr = c(57.027,  57.027, 57.027, 57.027, 57.027, 57.027, 57.027, 57.027, 57.027,  57.027, 57.027, 57.027, 57.027, 57.027, 57.027, 57.027, 57.027,  57.027, 57.027, 57.027, 57.027, 57.027, 57.027, 57.027, 57.027,  57.027, 57.027, 57.027, 57.027, 57.027, 57.027)), .Names = c("time",  "depot", "cen", "k", "ktr"), row.names = c(1L, 2L, 3L, 4L, 5L,  6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L, 179L, 180L, 181L,  182L, 183L, 184L, 185L, 186L, 187L, 188L, 189L, 190L, 191L, 192L,  193L, 194L), class = "data.frame");

dat2 <- structure(list(time = c(0, 0.0503, 0.1005, 0.1508, 0.201, 0.2513,  0.3015, 0.3518, 0.402, 0.4523, 0.5025, 0.5528, 0.603, 0.6533,  0.7035, 8.9447, 8.995, 9.0452, 9.0955, 9.1457, 9.196, 9.2462,  9.2965, 9.3467, 9.397, 9.4472, 9.4975, 9.5477, 9.598, 9.6482,  9.6985), depot = c(20, 19.6217, 19.2506, 18.8909, 18.6641, 19.2748,  21.8143, 26.0739, 30.3191, 33.1243, 34.3238, 34.4641, 34.0962,  33.539, 32.9284, 1.4375, 1.4103, 1.3836, 1.3575, 1.3318, 1.3066,  1.2819, 1.2576, 1.2339, 1.2105, 1.1876, 1.1652, 1.1431, 1.1215,  1.1003, 1.0795), cen = c(0, 0.3747, 0.7351, 1.0818, 1.4159, 1.7456,  2.0976, 2.5093, 2.9965, 3.5425, 4.1153, 4.6888, 5.2487, 5.7888,  6.3075, 4.7486, 4.6853, 4.6228, 4.561, 4.4999, 4.4394, 4.3796,  4.3205, 4.262, 4.2042, 4.1471, 4.0906, 4.0348, 3.9796, 3.9251,  3.8712), k = c(0.3814, 0.3814, 0.3814, 0.3814, 0.3814, 0.3814,  0.3814, 0.3814, 0.3814, 0.3814, 0.3814, 0.3814, 0.3814, 0.3814,  0.3814, 0.3814, 0.3814, 0.3814, 0.3814, 0.3814, 0.3814, 0.3814,  0.3814, 0.3814, 0.3814, 0.3814, 0.3814, 0.3814, 0.3814, 0.3814,  0.3814), ktr = c(57.027, 57.027, 57.027, 57.027, 57.027, 57.027,  57.027, 57.027, 57.027, 57.027, 57.027, 57.027, 57.027, 57.027,  57.027, 57.027, 57.027, 57.027, 57.027, 57.027, 57.027, 57.027,  57.027, 57.027, 57.027, 57.027, 57.027, 57.027, 57.027, 57.027,  57.027)), .Names = c("time", "depot", "cen", "k", "ktr"), row.names = c(1L,  2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L,  179L, 180L, 181L, 182L, 183L, 184L, 185L, 186L, 187L, 188L, 189L,  190L, 191L, 192L, 193L, 194L), class = "data.frame")

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
        expect_equal(round(as.data.frame(transit[c(1:15,seq(194-15,194)),]), 4),
                     dat1);
    })

    transit <- rxSolve(mod, et)

    test_that("Transit absorption is not specified, but still works.", {
        expect_equal(round(as.data.frame(transit[c(1:15,seq(194-15,194)),]), 4),
                     dat1);
        expect_warning(rxSolve(mod, et));
    })

    no_transit <- rxSolve(mod, et, transit_abs=FALSE);

    test_that("Transit absorption is turned off, and gives other results", {
        expect_equal(round(as.data.frame(no_transit[c(1:15,seq(194-15,194)),]), 4),
                     dat2);
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
        expect_equal(round(as.data.frame(transit[c(1:15,seq(194-15,194)),]), 4),
                     dat1);
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
        expect_equal(round(as.data.frame(transit[c(1:15,seq(194-15,194)),]), 4), dat1);
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
        expect_equal(round(as.data.frame(transit[c(1:15,seq(194-15,194)),]), 4), dat1)
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
        expect_equal(round(as.data.frame(transit[c(1:15,seq(194-15,194)),]), 4), dat1)
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
        expect_equal(round(as.data.frame(transit[c(1:15,seq(194-15,194)),]), 4), dat1);
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
