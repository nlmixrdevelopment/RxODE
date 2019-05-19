context("Test event Table et(...)");
rxPermissive({

    library(units)
    library(dplyr)


    et <- et()

    test_that("Empty event table check",{
        expect_equal(class(et$env),"rxHidden")
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

    et <- et(1:10, cmt=1);
    test_that("compartment check",{
        expect_equal(et$nobs, 10L);
        expect_equal(et$ndose, 0L);
        expect_equal(class(et$get.EventTable()), "data.frame");
        expect_true(is(et, "rxEt"));
        expect_false(et$show["id"])
        expect_true(et$show["cmt"])
        expect_equal(et$get.dosing(), NULL);
        expect_equal(length(et$get.sampling()$time), 10);
        expect_true(all(et$cmt==1L));
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

    test_that("Using simulate with et without windows will warn",{
        f1 <- as.data.frame(et3)
        f2 <- suppressWarnings({as.data.frame(simulate(et3))})
        expect_equal(f1$time,f2$time);
        expect_warning(simulate(et3));
    })



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

        expect_equal(structure(list(time = c(0, 216, 432),
                                    amt = c(100,200, 100),
                                    ii = c(24, 24, 24),
                                    addl = c(6L,6L, 6L),
                                    evid = c(1L, 1L, 1L)),
                               class = "data.frame",
                               row.names = c(NA,-3L)), e4)

        e5 <- etSeq(e1, wait=72, e2, wait=72, e1, waitII="+ii") %>%
            as.data.frame

        expect_equal(structure(list(time = c(0, 240, 480),
                                    amt = c(100,200, 100),
                                    ii = c(24, 24, 24),
                                    addl = c(6L,6L, 6L),
                                    evid = c(1L, 1L, 1L)),
                               class = "data.frame",
                               row.names = c(NA, -3L)), e5)

        e1 <- et(amt=500)
        e2 <- et(amt=250, ii=24, addl=4)

        expect_equal(structure(list(time = c(0, 24),
                                    amt = c(500, 250),
                                    ii = c(0, 24),
                                    addl = c(0L, 4L),
                                    evid = c(1L, 1L)),
                               class = "data.frame",
                               row.names = c(NA, -2L)),
                     c(e1, e2) %>% as.data.frame)

        expect_equal(structure(list(time = c(0, 120, 144),
                                    amt = c(250, 500, 250),
                                    ii = c(24, 0, 24),
                                    addl = c(4L, 0L, 4L),
                                    evid = c(1L, 1L, 1L)),
                               class = "data.frame",
                               row.names = c(NA, -3L)),
                     c(e2,e1,e2) %>% as.data.frame)

        e3 <- et(amt=200)

        expect_warning(c(e1,e3))

        e4 <- suppressWarnings(c(e1,e3) %>% as.data.frame)

        expect_equal(structure(list(time = c(0, 24),
                                    amt = c(500, 200),
                                    ii = c(0, 0),
                                    addl = c(0L, 0L),
                                    evid = c(1L, 1L)),
                               class = "data.frame",
                               row.names = c(NA, -2L)),
                     e4)

        e4 <- suppressWarnings(c(e1,e3,ii=12) %>% as.data.frame)

        expect_equal(structure(list(time = c(0, 12),
                                    amt = c(500, 200),
                                    ii = c(0, 0),
                                    addl = c(0L, 0L),
                                    evid = c(1L, 1L)),
                               class = "data.frame",
                               row.names = c(NA, -2L)),
                     e4)


        e1 <- et(amt=100, ii=24, addl=6) %>%
            et(seq(0, 2*168,by=0.1));
        e2 <- c(e1,e1,samples="use")
        expect_equal(range(e2$time), c(0, 672))

        ## combine without changing, use rbind
        e1 <- et(amt=100, ii=24, addl=6,  ID=1:5)
        e2 <- et(amt=50,  ii=12, addl=13, ID=1:3)
        e3 <- et(amt=200, ii=24, addl=2,  ID=1:2)

        e4 <- rbind(e1,e2,e3)
        expect_equal(e4 %>% select(id, time, amt, ii, addl) %>% as.data.frame,
                     structure(list(id = c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 3L, 4L, 5L),
                                    time = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                    amt = c(100, 50, 200, 100, 50, 200, 100, 50, 100, 100),
                                    ii = c(24, 12, 24, 24, 12, 24, 24, 12, 24, 24),
                                    addl = c(6L, 13L, 2L, 6L, 13L, 2L, 6L, 13L, 6L, 6L)),
                               class = "data.frame",
                               row.names = c(NA, -10L)))

        e4 <- rbind(e1, e2, e3, id="unique")
        expect_equal(e4 %>% select(id, time, amt, ii, addl) %>% as.data.frame,
                     structure(list(id = 1:10,
                                    time = c(0, 0, 0, 0, 0, 0, 0, 0, 0,0),
                                    amt = c(100, 100, 100, 100, 100, 50, 50, 50, 200, 200),
                                    ii = c(24, 24, 24, 24, 24, 12, 12, 12, 24, 24),
                                    addl = c(6L, 6L, 6L, 6L,6L, 13L, 13L, 13L, 2L, 2L)),
                               class = "data.frame",
                               row.names = c(NA, -10L)))

    })

    ## FIXME test windows and dur
    test_that("no ii throws an error with addl",{
        expect_error(et(c(11,13),amt=10,addl=3));
    })
    ## Test Windows
    et <- et(list(c(0,1),
                  c(1,4),
                  c(4,8),
                  c(8,12),
                  c(12,24))) %>%
        et(amt=10) %>%
        et(c(11,13),amt=10,addl=3,ii=12)

    test_that("Using simulate with et works and gives a different data frame",{
        et2 <- as.data.frame(simulate(et))
        et1 <- as.data.frame(et)
        expect_false(all(et1$time==et2$time))
        et3 <- as.data.frame(et);
        expect_true(all(et1$time==et3$time))
        et$simulate();
        et3 <- as.data.frame(et)
        expect_false(all(et1$time==et3$time))

    })

    test_that("Low/High middle tests; i.e windows",{
        et2 <- et[which(!is.na(et$low)),]
        expect_true(et$show["low"])
        expect_true(et$show["high"])
        expect_true(all(et2$time < et2$high))
        expect_true(all(et2$time > et2$low))
    })

    et <- et(list(c(0,1),
                  c(1,4),
                  c(4,8),
                  c(8,12),
                  c(12,24)), cmt=1) %>%
        et(amt=10) %>%
        et(c(11,13),amt=10,addl=3,ii=12)

    test_that("Low/High middle tests; i.e windows",{
        et2 <- et[which(!is.na(et$low)),]
        expect_true(et$show["low"])
        expect_true(et$show["high"])
        expect_true(all(et2$time < et2$high))
        expect_true(all(et2$time > et2$low))
    })

    test_that("Window Errors", {
        expect_error(et(list(c(1,0))))
        expect_error(et(list(c(0,1,2))))
    })
    context("until is inclusive")
    test_that("until is inclusive",{
        expect_equal(et(amt=1, time=50, until=57.5, ii=1.5)$addl, 5)
        expect_equal(et(amt=1, time=50, until=57.49999, ii=1.5)$addl, 4)
        expect_equal(et(amt=1, time=50, until=57.50001, ii=1.5)$addl, 5)
    })

    context("et Expected errors")
    test_that("et errors",{
        expect_error(et(list(c(2,1),c(3,4)), amt=3))
        expect_error(et(list(c(1,2),c(3), c(1,2,3)), amt=3))
        expect_error(et(list(c(1,2),c(3), TRUE), amt=3))
    })


}, silent=TRUE, cran=TRUE)
