rxPermissive({

    ## Test case for SS=2
    ## https://github.com/metrumresearchgroup/mrgsolve/issues/221
    pk2cmt <- RxODE({
        ## - Two-compartment PK model
        ##      - Dual first-order absorption
        ##      - Optional nonlinear clearance from `CENT`
        ##  - Source: Translated from mrgsolve to RxODE
        CL <- 1  # Clearance (volume/time)
        VC <- 3  # Central volume (volume)
        Q <- 2  # Inter-compartmental clearance (volume/time)
        VP <- 100  # Peripheral volume of distribution (volume)
        KA1 <- 1  # Absorption rate constant 1 (1/time)
        KA2 <- 1  # Absorption rate constant 2 (1/time)
        VMAX <- 0  # Maximum velocity (mass/time)
        KM <- 2  # Michaelis Constant (mass/volume)
        CP <- (CENT/VC)
        CT <- (PERIPH/VP)
        CLNL <- (VMAX/(KM+CP))
        d/dt(EV1)= -KA1*EV1;
        d/dt(EV2)= -KA2*EV2;
        d/dt(CENT)= KA1*EV1 + KA2*EV2 - (CL+CLNL+Q)*CP  + Q*CT;
        d/dt(PERIPH)= Q*CP - Q*CT;
    })

    m <- "liblsoda"

    ## m <- "lsoda"

    ## m <- "dop853"

    et <- eventTable() %>% add.dosing(dose=3, nbr.doses=6, dosing.interval=8) %>%
        add.sampling(seq(0, 48, length.out=200))

    ode.1c <- RxODE({
        V <- 20
        Cl <- 1
        C2 = center/V;
        d/dt(center) ~ - Cl*C2
    })

    tmp <- et$get.EventTable()
    tmp <- tmp[order(tmp$evid),]

    x2 <- solve(ode.1c, et, method=m)
    library(ggplot2); x2 %>% ggplot(aes(time, C2)) + geom_line()

    et2 <- eventTable() %>% add.dosing(dose=3, nbr.doses=2, dosing.interval=8) %>%
        add.sampling(seq(0, 48, length.out=200))

    x2 <- solve(ode.1c, et2, method=m)

    library(ggplot2); x2 %>% ggplot(aes(time, C2)) + geom_line()

    ## now setup bolus steady state evid

    et3 <- et2$get.EventTable();

    et3$evid[1] <- 110 # SS=1
    et3$ii <- 0.0
    et3$ii[1] <- 8

    x2 <- solve(ode.1c, et3, method=m)

    library(ggplot2); x2 %>% ggplot(aes(time, C2)) + geom_line()


    et4 <- et2$get.EventTable();

    et4$evid[et4$time == 8] <- 110
    et4$ii <- 0.0
    et4$ii[et4$time == 8] <- 8

    x2 <- solve(ode.1c, et4, method=m)


    library(ggplot2); x2 %>% ggplot(aes(time, C2)) + geom_line()

    et <- eventTable(time.units="days")
    et$add.sampling(seq(0,10,by=1/24))
    et$add.dosing(dose=2/24,rate=2,start.time=0,
                  nbr.doses=10,dosing.interval=0.4)

    mod <- RxODE({
        a = 6
        b = 0.6
        d/dt(intestine) = -a*intestine
        d/dt(blood)     = a*intestine - b*blood
    })

    et <- et$get.EventTable()

    x2 <- solve(mod, et, method=m)
    library(ggplot2); x2 %>% ggplot(aes(time, intestine)) + geom_line()

    library(ggplot2); x2 %>% ggplot(aes(time, blood)) + geom_line()


    ## w <- with(et,which(evid !=0 & amt==2))[2]
    w <- 1;
    et$evid[w] <- 10110
    et$ii <- 0;
    et$ii[w] <- 0.4;

    x2 <- solve(mod, et, method=m)
    library(ggplot2); x2 %>% ggplot(aes(time, intestine)) + geom_line()


})
