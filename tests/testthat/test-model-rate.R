rxPermissive({

    library(dplyr);
    library(testthat);
ms <- c("liblsoda", "lsoda", "dop853")
    if (grepl('SunOS',Sys.info()['sysname'])) ms <- "lsoda"

    for (m in ms){
        context(sprintf("Modeled rate (%s)",m))

        mod.rate <- RxODE({
            a = 6
            b = 0.6
            f = 1
            ri = 10
            li = 0
            d/dt(intestine) = -a*intestine
            f(intestine) = f
            rate(intestine) = ri
            lag(intestine) = li
            d/dt(blood)     = a*intestine - b*blood
        })

        library(ggplot2)

        et <- et() %>%
            et(amt=2/24,rate=-1, 0, addl=9, ii=1) %>%
            et(seq(0,10,by=1/24))


        et2 <- et() %>%
            et(amt=2/24,rate=2, 0, addl=9, ii=1) %>%
            et(seq(0,10,by=1/24))

        f1.a <- rxSolve(mod.rate, et, c(ri=2, f=1, li=0),method=m);
        f1.b <- rxSolve(mod.rate, et2, c(ri=0.1, f=1, li=0),method=m);

        test_that("rate=2, f=1, lag=0 makes sense for modeled rate",{
            expect_equal(seq(0,10,by=1/24),f1.a$time)
            expect_equal(seq(0,10,by=1/24),f1.b$time)
            expect_equal(as.data.frame(f1.a), as.data.frame(f1.b));
        })

        f1.a <- rxSolve(mod.rate, et, c(ri=2, f=1, li=0.5),method=m);
        f1.b <- rxSolve(mod.rate, et2, c(ri=0.1, f=1, li=0.5),method=m);

        test_that("rate=2, f=1, lag=0.5 makes sense for modeled rate",{
            expect_equal(seq(0,10,by=1/24),f1.a$time)
            expect_equal(seq(0,10,by=1/24),f1.b$time)
            expect_equal(as.data.frame(f1.a), as.data.frame(f1.b));
        })

        f1.a <- rxSolve(mod.rate, et, c(ri=2, f=2, li=0.5),method=m);
        f1.b <- rxSolve(mod.rate, et2, c(ri=0.1, f=2, li=0.5),method=m);

        test_that("rate=2, f=2, lag=0.5 makes sense for modeled rate & fixed rate",{
            expect_equal(seq(0,10,by=1/24),f1.a$time)
            expect_equal(seq(0,10,by=1/24),f1.b$time)
            expect_equal(as.data.frame(f1.a), as.data.frame(f1.b));
        })

        et <- et() %>% et(amt=2/24,rate=-1) %>%
            et(amt=0.2,time=1) %>%
            et(seq(0,10,by=1/24))

        et2 <- et() %>% et(amt=2/24,rate=2) %>%
            et(amt=0.2,time=1) %>%
            et(seq(0,10,by=1/24))

        f.a <- rxSolve(mod.rate, et, c(ri=2,f=2,li=0.3),method=m)
        f.b <- rxSolve(mod.rate, et2, c(ri=1e-6,f=2,li=0.3),method=m)

        test_that("rate=2 makes sense for modeled rate",{
            expect_equal(seq(0,10,by=1/24),f.a$time)
            expect_equal(seq(0,10,by=1/24),f.b$time)
            expect_equal(as.data.frame(f.a), as.data.frame(f.b));
        })

        f.a <- rxSolve(mod.rate, et, c(ri=0.125,f=2,li=0.3),method=m)

        et2 <- et() %>% et(amt=2/24,rate=0.125) %>%
            et(amt=0.2,time=1) %>%
            et(seq(0,10,by=1/24))

        f.b <- rxSolve(mod.rate, et2, c(ri=1e-6,f=2,li=0.3),method=m)

        test_that("rate=0.125, f=2, lag=0.3 with 2 doses makes sense for modeled rate",{
            expect_equal(seq(0,10,by=1/24),f.a$time)
            expect_equal(seq(0,10,by=1/24),f.b$time)
            expect_equal(round(as.data.frame(f.a),4), round(as.data.frame(f.b),4));
        })

        f.a <- rxSolve(mod.rate, et, c(ri=0.05,f=2,li=0.3),method=m)

        et2 <- et() %>% et(amt=2/24,rate=0.05) %>%
            et(amt=0.2,time=1) %>%
            et(seq(0,10,by=1/24))

        f.b <- rxSolve(mod.rate, et2, c(ri=1,f=2,li=0.3),method=m)

        et3 <- et() %>% et(amt=2/24*2,rate=0.05) %>%
            et(amt=0.2*2,time=1) %>%
            et(seq(0,10,by=1/24))

        f.c <- rxSolve(mod.rate, et3, c(ri=1,f=1,li=0.3),method=m)

        test_that("rate=0.05, f=2, lag=0.3 makes sense for modeled rate",{
            expect_equal(seq(0, 10, by=1/24),f.a$time)
            expect_equal(seq(0, 10, by=1/24),f.b$time)
            expect_equal(as.data.frame(f.a), as.data.frame(f.b));
        })

        test_that("bad rates (zero/negative) throw errors",{
            expect_error(rxSolve(mod.rate, et, c(ri=-1, f=1, li=0.3),method=m))
            expect_error(rxSolve(mod.rate, et, c(ri=0, f=1, li=0.3),method=m))
        })

        mod.dur <- RxODE({
            a = 6
            b = 0.6
            f = 1
            di = 3
            li = 0
            d/dt(intestine) = -a*intestine
            f(intestine) = f
            dur(intestine) = di
            lag(intestine) = li
            d/dt(blood)     = a*intestine - b*blood
        })

        test_that("Error when rate is requested but not in table",{
            expect_error(rxSolve(mod.dur, et, c(ri=1,f=1,li=0.3),method=m))
        })

        context(sprintf("Modeled duration (%s)", m))

        ## Now model duration
        et <- et(amt=2/24*2,rate=-2) %>%
            et(amt=0.2*2, time=2) %>%
            et(seq(0, 10, by=1/24))

        test_that("Error is thrown without modeled duration",{
            expect_error({rxSolve(mod.rate, et, c(ri=1,f=1,li=0.3),method=m)})
        })

        f.a <- rxSolve(mod.dur, et, c(di=0.5,f=1,li=0.3),method=m)

        et2 <- et() %>%
            et(amt=2/24*2,dur=0.5) %>%
            et(amt=0.2*2,time=2) %>%
            et(seq(0,10,by=1/24))

        f.b <- rxSolve(mod.dur, et2, c(di=10,f=1,li=0.3),method=m)

        test_that("Duration 0.5, F=1, li=0.3", {
            expect_equal(as.data.frame(f.a),as.data.frame(f.b))
        })

        et2 <- et() %>%
            et(amt=2/24*2,dur=0.5) %>%
            et(amt=0.2*2,time=2) %>%
            et(seq(0,10,by=1/24))

        f.a <- rxSolve(mod.dur, et, c(di=0.5,f=1,li=0.3),method=m)

        f.b <- rxSolve(mod.dur, et2, c(di=10,f=1,li=0.3),method=m)

        test_that("Duration 0.5, F=1, li=0.3", {
            expect_equal(as.data.frame(f.a),as.data.frame(f.b))
        })

        f.a <- rxSolve(mod.dur, et, c(di=0.5,f=2,li=0.3),method=m)

        f.b <- rxSolve(mod.dur, et2, c(di=10,f=2,li=0.3),method=m)

        test_that("Duration 0.5, F=1, li=0.3", {
            expect_equal(as.data.frame(f.a),as.data.frame(f.b))
        })
    }
}, cran=TRUE)


## et <- eventTable(time.units="days")
## et$add.sampling(seq(0,10,by=1/24))
## et$add.dosing(dose=2/24,rate=2,start.time=0,
##               nbr.doses=1,dosing.interval=1)

## et1 <- et$get.EventTable() %>% filter(is.na(amt) | amt > 0) %>%
##     mutate(evid=ifelse(evid==0, 0, 90101)) %>%
##     mutate(amt=ifelse(is.na(amt),NA_real_,0.4))

## et2 <- et1 %>% filter(evid !=0) %>% mutate(evid=70101)

## et3 <- rbind(et1, et2, data.frame(time=1,evid=101,amt=0.2)) %>% arrange(time,-evid)

## f <- rxSolve(mod.rate, et3, c(ri=2,f=2,li=0.3)); ggplot(f,aes(time,intestine)) + geom_line()

## f <- rxSolve(mod.rate, et3, c(ri=0.5,f=2,li=0.3)); ggplot(f,aes(time,intestine)) + geom_line()
