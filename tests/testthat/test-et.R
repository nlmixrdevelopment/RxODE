context("Test event Table et(...)");
rxPermissive({

    et <- et()

    test_that("Empty event table check",{
        expect_equal(et$nobs, 0L)
        expect_equal(et$ndose, 0L)
        expect_equal(et$get.EventTable(), NULL);
        expect_equal(et$get.dosing(), NULL);
        expect_equal(et$get.sampling(), NULL);
    })

    et <- et(1:10);

    test_that("10 sample items", {
        expect_equal(et$nobs, 10L);
        expect_equal(et$ndose, 0L);
        expect_equal(class(et$get.EventTable()), "data.frame");
        expect_true(is(et, "rxEt"));
        expect_false(et$show["id"])
        expect_false(et$show["cmt"])
        expect_equal(et$get.dosing(), NULL);
        expect_equal(length(et$get.sampling()$time), 10);
    })

    et <- et(1:10, "matt");

    test_that("compartment check",{
        expect_equal(et$nobs, 10L);
        expect_equal(et$ndose, 0L);
        expect_equal(class(et$get.EventTable()), "data.frame");
        expect_true(is(et, "rxEt"));
        expect_false(et$show["id"])
        expect_true(et$show["cmt"])
        expect_equal(et$get.dosing(), NULL);
        expect_equal(length(et$get.sampling()$time), 10);
    })

    et1 <- et(1:10, id=1:10)

    test_that("Observation only table check",{
        expect_equal(et1$nobs, 100L);
        expect_equal(et1$ndose, 0L);
        expect_equal(class(et1$get.EventTable()), "data.frame");
        expect_true(is(et1, "rxEt"));
        expect_true(et1$show["id"])
        expect_false(et1$show["cmt"])
        expect_equal(et1$get.dosing(), NULL);
        expect_equal(length(et1$get.sampling()$time), 100);
    })

    ## now resize down
    et2 <- et1 %>% et(id=-(2:10))

    test_that("compartment check",{
        expect_equal(et2$nobs, 10L);
        expect_equal(et2$ndose, 0L);
        expect_equal(class(et2$get.EventTable()), "data.frame");
        expect_true(is(et, "rxEt"));
        expect_true(et2$show["id"])
        expect_false(et2$show["cmt"])
        expect_equal(et2$get.dosing(), NULL);
        expect_equal(length(et2$get.sampling()$time), 10);
    })


    ## now resize back up
    et3 <- et2 %>% et(id=1:10)

    library(units)

    et3 <- et3 %>% set_units(mg);
    ## Make sure units are right.

    test_that("Observation only table check",{
        expect_equal(et1$nobs, 100L);
        expect_equal(et1$ndose, 0L);
        expect_equal(class(et1$get.EventTable()), "data.frame");
        expect_true(is(et1, "rxEt"));
        expect_true(et1$show["id"])
        expect_false(et1$show["cmt"])
        expect_equal(et1$get.dosing(), NULL);
        expect_equal(length(et1$get.sampling()$time), 100);
    })

    ## Check adding different units of time, rate, amt work


    test_that("units tests", {

        library(units)
        e <- et(amount.units="mg", time_units="hr") %>%
            add.sampling(seq(0,24,by=3)) %>%
            add.dosing(1,3) %>%
            et(rate=3,amt=2,time=120)
        e2 <- e %>% et(amt=set_units(0.0003, "lb"), time=0.5)
        expect_equal(e2$amt[e2$time == set_units(0.5,hr)], set_units(set_units(0.0003,lb),mg))
        e2 <- e %>% et(set_units(30,min))
        expect_true(any(e2$time == set_units(0.5,hr)))
        e2 <- e %>% et(time=0.25, rate=set_units(30,ug/min),amt=set_units(4,ug))
        tmp <- e2[e2$time == set_units(0.25,hr),];
        expect_equal(set_units(1.8,mg/h), tmp$rate)
        expect_equal(set_units(0.004,mg), tmp$amt)
        e2 <- e %>% et(time=0.25,ii=set_units(30,min), amt=4,addl=4)
        expect_equal(e2$ii[e2$time == set_units(0.25,hr)],set_units(0.5,hr));

        ## Check importing wrong different ii and time units as well as different rate units work.
        e <- et(amount.units="mg", time_units="hr") %>%
            add.sampling(seq(0,24,by=3)) %>%
            add.dosing(1,3) %>%
            et(rate=3,amt=2,time=120)

        etDf <- as.data.frame(e)
        etDf$rate <- set_units(etDf$rate,ug/s)
        etDf$ii <- set_units(etDf$ii, min)

        et <- et();
        et$import.EventTable(etDf);

        expect_equal(et$ii, e$ii)
        expect_equal(et$rate, e$rate)

    })

    test_that("seq works with wait", {
        e1 <- et(amt=100, ii=24, addl=6)

        e2 <- et(amt = 200, ii = 24, addl = 6)

        e3 <- seq(e1,e2,e1)

        e4 <- seq(e1, wait=72, e2, wait=72, e1) %>% as.data.frame

        expect_equal(structure(list(id = c(1L, 1L, 1L), low = c(NA_real_, NA_real_,
NA_real_), time = c(0, 216, 432), high = c(NA_real_, NA_real_,
NA_real_), cmt = c("(default)", "(default)", "(default)"), amt = c(100,
200, 100), rate = c(0, 0, 0), ii = c(24, 24, 24), addl = c(6L,
6L, 6L), evid = c(1L, 1L, 1L), ss = c(0L, 0L, 0L)), class = "data.frame", row.names = c(NA,
-3L)), e4)

        e5 <- etSeq(e1, wait=72, e2, wait=72, e1, handleWait="alwaysAddII") %>%
            as.data.frame

        expect_equal(structure(list(id = c(1L, 1L, 1L), low = c(NA_real_, NA_real_,
NA_real_), time = c(0, 240, 480), high = c(NA_real_, NA_real_,
NA_real_), cmt = c("(default)", "(default)", "(default)"), amt = c(100,
200, 100), rate = c(0, 0, 0), ii = c(24, 24, 24), addl = c(6L,
6L, 6L), evid = c(1L, 1L, 1L), ss = c(0L, 0L, 0L)), class = "data.frame", row.names = c(NA,
-3L)), e5)

    })

})
