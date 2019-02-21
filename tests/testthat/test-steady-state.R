rxPermissive({

    for (m in c("liblsoda", "lsoda", "dop853")){

        et <- eventTable() %>% add.dosing(dose=3, nbr.doses=6, dosing.interval=8) %>%
            add.sampling(seq(0, 48, length.out=200))

        ode.1c <- RxODE({
            V <- 20
            Cl <- 1
            C2 = center/V;
            d/dt(center) ~ - Cl*C2
        })

        d <- 3
        et2 <- eventTable() %>% add.dosing(dose=d, nbr.doses=2, dosing.interval=8) %>%
            add.sampling(seq(0, 48, length.out=200))

        x2 <- solve(ode.1c, et2, method=m)

        c0 <- setNames(d/rxInits(ode.1c)["V"],NULL);

        ke <- with(as.list(rxInits(ode.1c)),{Cl/V})

        context(sprintf("Steady state (%s)", m))

        test_that("Non steady state dose makes sense",{
            expect_equal(x2$C2[1], c0)
        })

        sigdig <- 4;

        for (ii in seq(2,24,by=2)){
            ## now setup bolus steady state evid
            et3 <- et() %>% et(amt=d, ss=1,ii=ii) %>%
                et(amt=3, time=8) %>% et(seq(0,  48, length.out=200))
            x2 <- solve(ode.1c, et3, method=m)
            test_that(paste("Steady State dose makes sense for ii=",ii),{
                expect_equal(round(x2$C2[1],sigdig), round(c0/(1-exp(-ke*ii)),sigdig))
            })
        }

        ## Need to test:
        ## - Fixed Infusion w/rate
        ## - Fixed Infusion w/rate + bioavailability change
        ## - Fixed Infusion w/dur
        ## - Fixed Infusion w/dur + bioavailability change
        ## - Modeled rate
        ## - Modeled rate + bioavailability change
        ## - Modeled duration
        ## - Modeled duration + bioavailability change
        for (dur in c(0.5, 1)){
            for (ii in seq(2,24,by=2)){
                et3 <- et() %>% et(amt=d, ss=1,ii=ii, rate=d/dur) %>%
                    et(c(dur,seq(0,  24, length.out=19)))
                x2 <- solve(ode.1c, et3, method=m,maxsteps=10000)
                test_that(paste("Infusion Steady State dose makes sense for ii=",ii," dur=",dur),{
                    expect_equal(round(x2$C2[x2 == dur], sigdig),
                                 round(with(as.list(rxInit(ode.1c)),
                                            d/(Cl*dur)*(1-exp(-ke*dur))/(1-exp(-ke*ii))),sigdig))
                })
            }
        }
    }
})
