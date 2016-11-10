context("Test Transit compartment model");
library(digest);
options(RxODE.verbose=FALSE);
rxClean();
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
d/dt(abs) = exp(log(bio*podo)+log(ktr)+n*log(ktr*t)-ktr*t-lgammafn(n+1))-ka*abs
d/dt(cen) = ka*abs-k*cen
")

et <- eventTable();
et$add.sampling(seq(0, 10, length.out=200));
et$add.dosing(20, start.time=0);

transit <- rxSolve(mod, et, transit_abs=TRUE)

test_that("Transit absorption works", {
    expect_equal(digest(round(as.data.frame(transit), 4)),
                 "07ba829ef75129d50b02578f76dc123b");
})

transit <- rxSolve(mod, et)

test_that("Transit absorption is not specified, but still works.", {
    expect_equal(digest(round(as.data.frame(transit), 4)),
                 "07ba829ef75129d50b02578f76dc123b");
    expect_warning(rxSolve(mod, et));
})

no_transit <- rxSolve(mod, et, transit_abs=FALSE);

test_that("Transit absorption is turned off, and gives other results", {
    expect_equal(digest(round(as.data.frame(no_transit), 4)),
                 "89b3c21d4367dd82da08063b5fe5883b");
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
d/dt(abs) = exp(log(bio*podo)+log(ktr)+n*log(ktr*t)-ktr*t-lgamma1p(n))-ka*abs
d/dt(cen) = ka*abs-k*cen
")

transit <- rxSolve(mod, et);

test_that("Transit absorption using lagmma1p works.", {
    expect_equal(digest(round(as.data.frame(transit), 4)),
                 "07ba829ef75129d50b02578f76dc123b");
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
d/dt(abs) = exp(log(bio*podo)+log(ktr)+n*log(ktr*t)-ktr*t-log(n!))-ka*abs
d/dt(cen) = ka*abs-k*cen
")

test_that("Transit absorption using factorial.", {
    expect_equal(digest(round(as.data.frame(transit), 4)),
                 "07ba829ef75129d50b02578f76dc123b");
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
d/dt(abs) = transit(n, mtt)-ka*abs
d/dt(cen) = ka*abs-k*cen
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
d/dt(abs) = transit(n, mtt, bio)-ka*abs
d/dt(cen) = ka*abs-k*cen
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
d/dt(abs) = transit(n, mtt, 1)-ka*abs
d/dt(cen) = ka*abs-k*cen
")

transit2 <- rxSolve(mod, et);

test_that("Transit absorption can take numeric arguments", {
    expect_equal(round(as.data.frame(transit)[,names(transit) != "ktr"], 4),
                 round(as.data.frame(transit2), 4));
})
