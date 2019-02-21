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

        ode.1cR <- RxODE({
            V <- 20
            Cl <- 1
            C2 = center/V;
            d/dt(center) ~ - Cl*C2
            rate(center) = rateIn
        })

        ode.1cD <- RxODE({
            V <- 20
            Cl <- 1
            C2 = center/V;
            d/dt(center) ~ - Cl*C2
            dur(center) = durIn
        })

        d <- 3
        et2 <- eventTable() %>% add.dosing(dose=d, nbr.doses=2, dosing.interval=8) %>%
            add.sampling(seq(0, 48, length.out=200))

        x2 <- solve(ode.1c, et2, method=m)

        c0 <- setNames(d/rxInits(ode.1c)["V"],NULL);

        ke <- with(as.list(rxInits(ode.1c)),{Cl/V})

        context(sprintf("Steady state IV Bolus (%s)", m))

        test_that("Non steady state dose makes sense",{
            expect_equal(x2$C2[1], c0)
        })

        tol <- 1e-3

        for (ii in seq(2,24,by=2)){
            ## now setup bolus steady state evid
            et3 <- et() %>% et(amt=d, ss=1,ii=ii) %>%
                et(amt=d, time=ii) %>% et(seq(0, 48, length.out=200))
            x2 <- solve(ode.1c, et3, method=m)
            test_that(paste("Steady State dose makes sense for ii=",ii),{
                expect_equal(x2$C2[1], c0/(1-exp(-ke*ii)), tolerance=tol)
            })
        }

        context(sprintf("Steady state Infusions (%s)", m))

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
                ## Fixed rate
                et3 <- et() %>% et(amt=d, ss=1,ii=ii, rate=d/dur) %>%
                    et(time=ii,amt=d,ii=ii,addl=floor(24/ii),rate=d/dur) %>%
                    et(c(dur,seq(0,  24, length.out=200)))
                x2 <- solve(ode.1c, et3, method=m,maxsteps=10000)
                infMax <- with(as.list(rxInit(ode.1c)),
                                     d/(Cl*dur)*(1-exp(-ke*dur))/(1-exp(-ke*ii)))
                inf0 <- with(as.list(rxInit(ode.1c)),
                             infMax*exp(-ke*(ii-dur)))
                test_that(paste("Infusion Steady State dose makes sense for ii=",ii," dur=",dur, "(rate)"),{
                    expect_equal(x2$C2[x2 == dur], infMax, tolerance=tol);
                    expect_equal(x2$C2[1], inf0, tolerance=tol)
                })
                ## Fixed duration
                et3 <- et() %>% et(amt=d, ss=1,ii=ii, dur=dur) %>%
                    et(c(dur,seq(0,  24, length.out=19)))
                x2 <- solve(ode.1c, et3, method=m,maxsteps=10000)
                test_that(paste("Infusion Steady State dose makes sense for ii=",ii," dur=",dur, "(dur)"),{
                    expect_equal(x2$C2[x2 == dur], infMax, tolerance=tol);
                    expect_equal(x2$C2[1], inf0, tolerance=tol)
                })
                ## rate modeled
                et3 <- et() %>% et(amt=d, ss=1,ii=ii, rate=-1) %>%
                    et(c(dur,seq(0,  24, length.out=19)))
                x2 <- rxSolve(ode.1cR, et3, c(rateIn=d/dur), method=m, maxsteps=10000)
                test_that(paste("Infusion Steady State dose makes sense for ii=",ii," dur=",dur, "(rate modeled)"),{
                    expect_equal(x2$C2[x2 == dur], infMax, tolerance=tol);
                    expect_equal(x2$C2[1], inf0, tolerance=tol)
                })
                ## duration modeled
                et3 <- et() %>% et(amt=d, ss=1,ii=ii, rate=-2) %>%
                    et(c(dur,seq(0,  24, length.out=19)))
                x2 <- rxSolve(ode.1cD, et3, c(durIn=dur), method=m, maxsteps=10000)
                test_that(paste("Infusion Steady State dose makes sense for ii=",ii," dur=",dur, "(rate modeled)"),{
                    expect_equal(x2$C2[x2 == dur], infMax, tolerance=tol);
                    expect_equal(x2$C2[1], inf0, tolerance=tol)
                })
            }
        }


    }
})
