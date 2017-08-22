rxPermissive({

    context("Test the solved equations")

    et <- eventTable() %>% add.dosing(dose=3, nbr.doses=6, dosing.interval=8) %>%
        add.sampling(seq(0, 48, length.out=200))

    rxClean()

    ode.1c <- RxODE({
        C2 = center/V;
        d/dt(center) = - CL*C2
    })

    ## Solved systems can check the variables in the RxODE statement
    ## to figure out what type of solved system is being requested
    ode.1cs <- RxODE({
        V <- theta[1];
        CL <- theta[2];
        C2 = linCmt();
    })

    ## Instead of specifying parameters in the solved system, you can
    ## specify them in the linCmt variable.
    ode.1cs2 <- RxODE({
        C2 = linCmt(CL, V);
    })

    ## The solved systems can be mixed with ODE solving routines (to
    ## speed them up a bit...?)

    o.1c <- ode.1c %>% solve(params=c(V=20, CL=25), events=et)

    s.1c <- ode.1cs2 %>% solve(params=c(V=20, CL=25), events=et)

    s.2c <- ode.1cs %>% solve(theta=c(20, 25), events=et)

    test_that("1 compartment solved models and ODEs same.", {
        expect_equal(round(o.1c$C2,4), round(s.1c$C2,4))
        expect_equal(round(o.1c$C2,4), round(s.2c$C2,4))
    })

    ode.1c.ka <- RxODE({
        C2 = center/V;
        d / dt(depot) = -KA * depot
        d/dt(center) = KA * depot - CL*C2
    })


    sol.1c.ka <- RxODE({
        C2 = linCmt(V, CL, KA);
    })

    o.1c <- ode.1c.ka %>% solve(params=c(V=20, CL=25, KA=2), events=et)

    s.1c <- sol.1c.ka %>% solve(params=c(V=20, CL=25, KA=2), events=et)

    test_that("1 compartment oral solved models and ODEs same.", {
        expect_equal(round(o.1c$C2,4), round(s.1c$C2,4))
    })

    ode.2c <- RxODE({
        C2 = centr/V;
        C3 = peri/V2;
        d/dt(centr) = - CL*C2 - Q*C2 + Q*C3;
        d/dt(peri)  = Q*C2 - Q*C3;
    })

    sol.2c <- RxODE({
        C2=linCmt(V, CL, V2, Q);
    })

    o.2c <- ode.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=et)

    s.2c <- sol.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=et)


    test_that("2 compartment solved models and ODEs same.", {
        expect_equal(round(o.2c$C2,4), round(s.2c$C2,4))
    })

    ode.2c.ka <- RxODE({
        C2 = centr/V;
        C3 = peri/V2;
        d/dt(depot) =-KA*depot;
        d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
        d/dt(peri)  =                    Q*C2 - Q*C3;
    })

    sol.2c.ka <- RxODE({
        C2=linCmt(V, CL, V2, Q, KA);
    })

    o.2c <- ode.2c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, KA= 0.3), events=et)

    s.2c <- sol.2c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, KA=0.3), events=et)

    test_that("2 compartment oral solved models and ODEs same.", {
        expect_equal(round(o.2c$C2,4), round(s.2c$C2,4))
    })


    ode.3c <- RxODE({
        C2 = centr/V;
        C3 = peri/V2;
        C4 = peri2 / V3
        d/dt(centr) = - CL*C2 - Q*C2 + Q*C3  - Q2*C2 + Q2*C4;
        d/dt(peri)  = Q*C2 - Q*C3;
        d / dt(peri2) = Q2 * C2 - Q2 * C4
    })

    sol.3c <- RxODE({
        ## double solvedC(double t, int parameterization, int cmt, unsigned int col, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8);
        C2=linCmt(V, CL, V2, Q, Q2, V3);
    })

    o.3c <- ode.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et)

    s.3c <- sol.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et)

    test_that("3 compartment solved models and ODEs same.", {
        expect_equal(round(o.3c$C2,4), round(s.3c$C2,4))
    })


    ode.3c.ka <- RxODE({
        C2 = centr/V;
        C3 = peri/V2;
        C4 = peri2 / V3
        d/dt(depot) =-KA*depot;
        d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3  - Q2*C2 + Q2*C4;
        d/dt(peri)  = Q*C2 - Q*C3;
        d / dt(peri2) = Q2 * C2 - Q2 * C4
    })

    sol.3c.ka <- RxODE({
        ## double solvedC(double t, int parameterization, int cmt, unsigned int col, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8);
        C2=linCmt(V, CL, V2, Q, Q2, V3, KA);
    })

    o.3c <- ode.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3), events=et)

    s.3c <- sol.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3), events=et)


    test_that("3 compartment oral solved models and ODEs same.", {
        expect_equal(round(o.3c$C2,4), round(s.3c$C2,4))
    })

    context("Infusion Models")

    et <- eventTable() %>% add.dosing(dose=3, rate=1.5, nbr.doses=6, dosing.interval=8) %>%
        add.sampling(seq(0, 48, length.out=200))

    ode.1c <- RxODE({
        C2 = center/V;
        d/dt(center) = - CL*C2
    })

    ## Solved systems can check the variables in the RxODE statement
    ## to figure out what type of solved system is being requested
    ode.1cs <- RxODE({
        V <- theta[1];
        CL <- theta[2];
        C2 = linCmt();
    })

    ## Instead of specifying parameters in the solved system, you can
    ## specify them in the linCmt variable.
    ode.1cs2 <- RxODE({
        C2 = linCmt(CL, V);
    })

    ## The solved systems can be mixed with ODE solving routines (to
    ## speed them up a bit...?)

    o.1c <- ode.1c %>% solve(params=c(V=20, CL=25), events=et)

    s.1c <- ode.1cs2 %>% solve(params=c(V=20, CL=25), events=et)

    s.2c <- ode.1cs %>% solve(theta=c(20, 25), events=et)

    test_that("1 compartment solved models and ODEs same.", {
        expect_equal(round(o.1c$C2,4), round(s.1c$C2,4))
        expect_equal(round(o.1c$C2,4), round(s.2c$C2,4))
    })

    ode.2c <- RxODE({
        C2 = centr/V;
        C3 = peri/V2;
        d/dt(centr) = - CL*C2 - Q*C2 + Q*C3;
        d/dt(peri)  = Q*C2 - Q*C3;
    })

    sol.2c <- RxODE({
        C2=linCmt(V, CL, V2, Q);
    })

    o.2c <- ode.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=et)

    s.2c <- sol.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=et)

    test_that("2 compartment solved models and ODEs same.", {
        expect_equal(round(o.2c$C2,4), round(s.2c$C2,4))
    })

    ode.3c <- RxODE({
        C2 = centr/V;
        C3 = peri/V2;
        C4 = peri2 / V3
        d/dt(centr) = - CL*C2 - Q*C2 + Q*C3  - Q2*C2 + Q2*C4;
        d/dt(peri)  = Q*C2 - Q*C3;
        d / dt(peri2) = Q2 * C2 - Q2 * C4
    })

    sol.3c <- RxODE({
        ## double solvedC(double t, int parameterization, int cmt, unsigned int col, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8);
        C2=linCmt(V, CL, V2, Q, Q2, V3);
    })

    o.3c <- ode.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et)

    s.3c <- sol.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et)

    test_that("3 compartment solved models and ODEs same.", {
        expect_equal(round(o.3c$C2,4), round(s.3c$C2,4))
    })

    context("Infusion + Bolus")

    et <- eventTable() %>% add.dosing(dose=3, rate=1.5, nbr.doses=6, dosing.interval=8) %>%
        add.dosing(dose=1.5, nbr.doses=6, dosing.interval=8) %>%
        add.sampling(seq(0, 48, length.out=200))

    ode.1c <- RxODE({
        C2 = center/V;
        d/dt(center) = - CL*C2
    })

    ## Solved systems can check the variables in the RxODE statement
    ## to figure out what type of solved system is being requested
    ode.1cs <- RxODE({
        V <- theta[1];
        CL <- theta[2];
        C2 = linCmt();
    })

    ## Instead of specifying parameters in the solved system, you can
    ## specify them in the linCmt variable.
    ode.1cs2 <- RxODE({
        C2 = linCmt(CL, V);
    })

    ## The solved systems can be mixed with ODE solving routines (to
    ## speed them up a bit...?)

    o.1c <- ode.1c %>% solve(params=c(V=20, CL=25), events=et)

    s.1c <- ode.1cs2 %>% solve(params=c(V=20, CL=25), events=et)

    s.2c <- ode.1cs %>% solve(theta=c(20, 25), events=et)

    test_that("1 compartment solved models and ODEs same.", {
        expect_equal(round(o.1c$C2,4), round(s.1c$C2,4))
        expect_equal(round(o.1c$C2,4), round(s.2c$C2,4))
    })

    ode.2c <- RxODE({
        C2 = centr/V;
        C3 = peri/V2;
        d/dt(centr) = - CL*C2 - Q*C2 + Q*C3;
        d/dt(peri)  = Q*C2 - Q*C3;
    })

    sol.2c <- RxODE({
        C2=linCmt(V, CL, V2, Q);
    })

    o.2c <- ode.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=et)

    s.2c <- sol.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=et)

    test_that("2 compartment solved models and ODEs same.", {
        expect_equal(round(o.2c$C2,4), round(s.2c$C2,4))
    })

    ode.3c <- RxODE({
        C2 = centr/V;
        C3 = peri/V2;
        C4 = peri2 / V3
        d/dt(centr) = - CL*C2 - Q*C2 + Q*C3  - Q2*C2 + Q2*C4;
        d/dt(peri)  = Q*C2 - Q*C3;
        d / dt(peri2) = Q2 * C2 - Q2 * C4
    })

    sol.3c <- RxODE({
        ## double solvedC(double t, int parameterization, int cmt, unsigned int col, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8);
        C2=linCmt(V, CL, V2, Q, Q2, V3);
    })

    o.3c <- ode.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et)

    s.3c <- sol.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et)

    test_that("3 compartment solved models and ODEs same.", {
        expect_equal(round(o.3c$C2,4), round(s.3c$C2,4))
    })

})
