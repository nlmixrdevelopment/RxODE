rxPermissive({

    library(dplyr)
    context("Test reset event EVID=3")

    ## 6.1
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

        x27 <- x2 %>% filter(time==7.5)

        zeros <- rep(0,length(x27$blood));
        print(zeros)

        test_that(sprintf("EVID=3 resets the system (%s)", m),{
            expect_true(any(x27$blood==zeros))
            expect_true(any(x27$intestine==zeros))
        })
    }

    ## library(ggplot2);x2 %>% ggplot(aes(time,blood)) + geom_line()
})

