rxPermissive({

    context("Test reset event EVID=3")

    ## 6.1
    library(dplyr)

    mod <- RxODE({
        a = 6
        b = 0.6
        d/dt(intestine) = -a*intestine
        d/dt(blood)     = a*intestine - b*blood
    })

    et <- eventTable(time.units="days")
    et$add.sampling(seq(0,10,by=1/24))
    et$add.dosing(dose=2/24,strt.time=0,
                  nbr.doses=10,dosing.interval=1)

    et1 <- et$get.EventTable()


    et2 <- rbind(et1,data.frame(time=7.5,evid=3, amt=NA)) %>%
        arrange(time,evid)

    for (m in c("lsoda","liblsoda", "dop853")){

        x2 <- solve(mod,et2, method=m)

        x27 <- x2 %>% filter(time>=7.5) %>% filter(time < 8)

        zeros <- rep(0,length(x27$blood));

        test_that(sprintf("EVID=3 resets the system (%s)", m),{
            expect_true(any(x27$blood==zeros))
            expect_true(any(x27$intestine==zeros))
        })
    }

    ## library(ggplot2);x2 %>% ggplot(aes(time,blood)) + geom_line()

    et <- eventTable() %>% add.dosing(dose=3, nbr.doses=6, dosing.interval=8) %>%
        add.sampling(seq(0, 48, length.out=200))

    sol.1c <- RxODE({
        V <- 20
        Cl <- 25
        blood <- linCmt();
    })

    x2 <- solve(sol.1c, et)

    x2 <- solve(sol.1c, et$get.EventTable())

    et1 <- rbind(et$get.EventTable(), c(time=9,evid=3,amt=NA_real_)) %>%
        arrange(time,-evid)

    x2 <- solve(sol.1c, et1)

    test_that("Solved Linear EVID=3",{
        expect_true(all((x2 %>% filter(time > 9) %>% filter(time < 12))$blood==0))
    })

    ## library(ggplot2);x2 %>% ggplot(aes(time,blood)) + geom_line()
})

