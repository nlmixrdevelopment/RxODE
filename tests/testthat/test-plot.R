rxPermissive({
    context("plot tests")

    ## Model from RxODE tutorial
    m1 <-RxODE({
        KA=2.94E-01;
        CL=1.86E+01;
        V2=4.02E+01;
        Q=1.05E+01;
        V3=2.97E+02;
        Kin=1;
        Kout=1;
        EC50=200;
        ## Added modeled bioavaiblity, duration and rate
        fdepot = 1;
        durDepot = 8;
        rateDepot = 1250;
        C2 = centr/V2;
        C3 = peri/V3;
        d/dt(depot) =-KA*depot;
        f(depot) = fdepot
        dur(depot) = durDepot
        rate(depot) = rateDepot
        d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
        d/dt(peri)  =                    Q*C2 - Q*C3;
        d/dt(eff)   = Kin - Kout*(1-C2/(EC50+C2))*eff;
        eff(0) = 1
    });

    ev <- et(timeUnits="hr") %>%
        et(amt=10000, ii=12,until=24) %>%
        et(seq(0, 24, length.out=100))

    s <- rxSolve(m1, ev);

    vdiffr::expect_doppelganger("plot-C2", s %>% plot(C2))

    vdiffr::expect_doppelganger("plot-C2-log-x", s %>% plot(C2, log="x"))

    vdiffr::expect_doppelganger("plot-C2-log-y", s %>% plot(C2, log="y"))

    vdiffr::expect_doppelganger("plot-C2-all-xy", s %>% plot(C2, log="xy"))

    vdiffr::expect_doppelganger("plot-C2-all-yx", s %>% plot(C2, log="yx"))

    vdiffr::expect_doppelganger("plot-all", s %>% plot())

    vdiffr::expect_doppelganger("plot-all-log-x", s %>% plot(log="x"))

    vdiffr::expect_doppelganger("plot-all-log-y", s %>% plot(log="y"))

    vdiffr::expect_doppelganger("plot-all-all-xy", s %>% plot(log="xy"))

    vdiffr::expect_doppelganger("plot-all-all-yx", s %>% plot(log="yx"))



}, silent=TRUE, test="plot")
