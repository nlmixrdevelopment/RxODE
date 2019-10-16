rxPermissive({
    tol  <- 5e-6 ## Current difference for all equations
    rxClean()

    for (ll in c(TRUE, FALSE)){
        context(sprintf("Test the solved equations (%s)", ifelse(ll, "linlog", "linear")))

        et <- eventTable() %>% add.dosing(dose=3, nbr.doses=6, dosing.interval=8) %>%
            add.sampling(seq(0, 48, length.out=200))


        ode.1c <- RxODE({
            C2 = center/V;
            d/dt(center) = - CL*C2
        })

        test_that("ode model gives extraCmt=0",{
            expect_equal(rxModelVars(ode.1c)$extraCmt,0L);
        })

        ## Solved systems can check the variables in the RxODE statement
        ## to figure out what type of solved system is being requested
        ode.1cs <- RxODE({
            V <- theta[1];
            CL <- theta[2];
            C2 = linCmt();
        })

        ode.2cK <- RxODE({
            V <- theta[1];
            CLx <- theta[2];
            K <- CLx/V
            C2 = linCmt();
        })

        ode.2cA1 <- RxODE({
            V <- theta[1];
            CLx <- theta[2];
            alpha <- CLx/V
            C2 = linCmt();
        })

        ode.2cA2 <- RxODE({
            A <- 1/theta[1];
            CLx <- theta[2];
            alpha <- CLx*A
            C2 = linCmt();
        })

        ## Instead of specifying parameters in the solved system, you can
        ## specify them in the linCmt variable.
        ode.1cs2 <- RxODE({
            C2 = linCmt(CL, V);
        })

        test_that("linear compartment model gives extraCmt=1",{
            expect_equal(rxModelVars(ode.1cs2)$extraCmt,1L);
        })

        ## The solved systems can be mixed with ODE solving routines (to
        ## speed them up a bit...?)

        o.1c <- ode.1c %>% solve(params=c(V=20, CL=25), events=et,linLog=ll)

        s.1c <- ode.1cs2 %>% solve(params=c(V=20, CL=25), events=et,linLog=ll)

        s.2c <- ode.1cs %>% solve(theta=c(20, 25), events=et,linLog=ll)
        s.2cK <- ode.2cK %>% solve(theta=c(20, 25), events=et,linLog=ll)
        s.2cA1 <- ode.2cA1 %>% solve(theta=c(20, 25), events=et,linLog=ll)
        s.2cA2 <- ode.2cA2 %>% solve(theta=c(20, 25), events=et,linLog=ll)

        test_that("Gives the correct parameters for THETAs",{
            expect_equal(s.2c$params,
                         structure(list("THETA[1]" = 20, "THETA[2]" = 25), class = "data.frame",
                                   row.names = c(NA, -1L)))
        })

        test_that("1 compartment solved models and ODEs same.", {
            expect_equal(o.1c$C2, s.1c$C2, tolerance=tol)
            expect_equal(o.1c$C2, s.2c$C2, tolerance=tol)
            expect_equal(o.1c$C2, s.2cK$C2, tolerance=tol)
            expect_equal(o.1c$C2, s.2cA1$C2, tolerance=tol)
            expect_equal(o.1c$C2, s.2cA2$C2, tolerance=tol)
        })

        ## Test steady state doses.
        etSs  <- et() %>% et(amt=3) %>%
            et(time=4,amt=3, ss=1, ii=24) %>%
            et(amt=3, ss=2, ii=24, time=8) %>%
            et(seq(0,24,length.out=200))

        o.1c <- ode.1c %>% solve(params=c(V=20, CL=1), events=etSs,linLog=ll)

        s.1c <- ode.1cs2 %>% solve(params=c(V=20, CL=1), events=etSs,linLog=ll)

        s.2c <- ode.1cs %>% solve(theta=c(20, 1), events=etSs,linLog=ll)

        test_that("1 compartment steady-state solved models and ODEs same.", {
            expect_equal(o.1c$C2, s.1c$C2, tolerance=tol)
            expect_equal(o.1c$C2, s.2c$C2, tolerance=tol)
        })


        ode.1c.ka <- RxODE({
            C2 = center/V;
            d / dt(depot) = -KA * depot
            d/dt(center) = KA * depot - CL*C2
        })


        sol.1c.ka <- RxODE({
            C2 = linCmt(V, CL, KA);
        })

        ode.2cK <- RxODE({
            V <- theta[1];
            CLx <- theta[2];
            Ka <- theta[3]
            K <- CLx/V
            C2 = linCmt();
        })

        ode.2cA1 <- RxODE({
            V <- theta[1];
            CLx <- theta[2];
            Ka <- theta[3]
            alpha <- CLx/V
            C2 = linCmt();
        })

        ode.2cA2 <- RxODE({
            A <- 1/theta[1];
            CLx <- theta[2];
            Ka <- theta[3];
            alpha <- CLx*A
            C2 = linCmt();
        })


        test_that("linear oral model gives extraCmt=2",{
            expect_equal(rxModelVars(sol.1c.ka)$extraCmt,2L);
        })

        o.1c <- ode.1c.ka %>% solve(params=c(V=20, CL=25, KA=2), events=et,linLog=ll)

        s.1c <- sol.1c.ka %>% solve(params=c(V=20, CL=25, KA=2), events=et,linLog=ll)

        s.2cK <- ode.2cK %>% solve(theta=c(20, 25, KA=2), events=et,linLog=ll)
        s.2cA1 <- ode.2cA1 %>% solve(theta=c(20, 25, KA=2), events=et,linLog=ll)
        s.2cA2 <- ode.2cA2 %>% solve(theta=c(20, 25, KA=2), events=et,linLog=ll)

        test_that("1 compartment oral solved models and ODEs same.", {
            expect_equal(o.1c$C2, s.1c$C2, tolerance=tol)
            expect_equal(o.1c$C2, s.2cK$C2, tolerance=tol)
            expect_equal(o.1c$C2, s.2cA1$C2, tolerance=tol)
            expect_equal(o.1c$C2, s.2cA2$C2, tolerance=tol)
        })

        ## Note the strange-looking dip at 4 hours.  This is because ss=1 resets the system first.
        o.1c <- ode.1c.ka %>% solve(params=c(V=20, CL=2, KA=2), events=etSs,linLog=ll)

        s.1c <- sol.1c.ka %>% solve(params=c(V=20, CL=2, KA=2), events=etSs,linLog=ll)

        test_that("1 compartment oral solved models steady state ODEs same.", {
            expect_equal(o.1c$C2, s.1c$C2, tolerance=tol)
        })

        ode.2c <- RxODE({
            C2 = centr/V;
            C3 = peri/V2;
            d/dt(centr) = - CL*C2 - Q*C2 + Q*C3;
            d/dt(peri)  = Q*C2 - Q*C3;
        })

        sol.2c <- RxODE({
            C2=linCmt(V, CL, V2, Q1);
        })

        sol.2cK <- RxODE({
            V <- theta[1]
            CLx <- theta[2]
            V2x <- theta[3]
            Q <- theta[4]
            K <- CLx/V
            K12 <- Q/V
            K21 <- Q/V2x
            C2=linCmt();
        })

        ## A1 in terms of A, alpha, B, beta

        sol.2cA1 <- RxODE({
            Vx <- theta[1]
            CLx <- theta[2]
            V2x <- theta[3]
            Q <- theta[4]
            Kx <- CLx/Vx
            K12x <- Q/Vx
            K21x <- Q/V2x
            beta <- 0.5 * (K12x + K21x + Kx -
                           sqrt((K12x +
                                   K21x + Kx) * (K12x + K21x + Kx) - 4 * K21x *
                                  Kx))
            alpha <- K21x * Kx/beta
            A <- (alpha - K21x)/(alpha - beta)/Vx
            B <- (beta - K21x)/(beta - alpha)/Vx
            C2=linCmt();
        })

        ## A2 V, alpha, beta, k21
        sol.2cA2 <- RxODE({
            V <- theta[1]
            CLx <- theta[2]
            V2x <- theta[3]
            Q <- theta[4]
            Kx <- CLx/V
            K12x <- Q/V
            K21 <- Q/V2x
            beta <- 0.5 * (K12x + K21 + Kx -
                           sqrt((K12x +
                                   K21 + Kx) * (K12x + K21 + Kx) - 4 * K21 *
                                  Kx))
            alpha <- K21 * Kx/beta
            C2=linCmt();
        })

        ## A3 alpha, beta, aob
        sol.2cA3 <- RxODE({
            V <- theta[1]
            CLx <- theta[2]
            V2x <- theta[3]
            Q <- theta[4]
            Kx <- CLx/V
            K12x <- Q/V
            K21x <- Q/V2x
            beta <- 0.5 * (K12x + K21x + Kx -
                           sqrt((K12x +
                                   K21x + Kx) * (K12x + K21x + Kx) - 4 * K21x *
                                  Kx))
            alpha <- K21x * Kx/beta
            Ax <- (alpha - K21x)/(alpha - beta)/V
            Bx <- (beta - K21x)/(beta - alpha)/V
            aob <- Ax/Bx
            C2=linCmt();
        })

        o.2c <- ode.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=et,linLog=ll)
        s.2c <- sol.2c %>% solve(params=c(V=40, CL=18, V2=297, Q1=10), events=et,linLog=ll)
        s.2cK <- sol.2cK %>% solve(theta=c(V=40, CL=18, V2=297, Q=10), events=et,linLog=ll)
        s.2cA1 <- sol.2cA1 %>% solve(theta=c(V=40, CL=18, V2=297, Q=10), events=et,linLog=ll)
        s.2cA2 <- sol.2cA2 %>% solve(theta=c(V=40, CL=18, V2=297, Q=10), events=et,linLog=ll)
        s.2cA3 <- sol.2cA3 %>% solve(theta=c(V=40, CL=18, V2=297, Q=10), events=et,linLog=ll)


        test_that("2 compartment solved models and ODEs same.", {
            expect_equal(o.2c$C2, s.2c$C2, tolerance=tol)
            expect_equal(o.2c$C2, s.2cK$C2, tolerance=tol)
            expect_equal(o.2c$C2, s.2cA1$C2, tolerance=tol)
            expect_equal(o.2c$C2, s.2cA2$C2, tolerance=tol)
            expect_equal(o.2c$C2, s.2cA3$C2, tolerance=tol)
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

        sol.2cK <- RxODE({
            V <- theta[1]
            CLx <- theta[2]
            V2x <- theta[3]
            Q <- theta[4]
            ka <- theta[5]
            K <- CLx/V
            K12 <- Q/V
            K21 <- Q/V2x
            C2=linCmt();
        })

        ## A1 in terms of A, alpha, B, beta

        sol.2cA1 <- RxODE({
            Vx <- theta[1]
            CLx <- theta[2]
            V2x <- theta[3]
            Q <- theta[4]
            ka <- theta[5]
            Kx <- CLx/Vx
            K12x <- Q/Vx
            K21x <- Q/V2x
            beta <- 0.5 * (K12x + K21x + Kx -
                           sqrt((K12x +
                                   K21x + Kx) * (K12x + K21x + Kx) - 4 * K21x *
                                  Kx))
            alpha <- K21x * Kx/beta
            A <- (alpha - K21x)/(alpha - beta)/Vx
            B <- (beta - K21x)/(beta - alpha)/Vx
            C2=linCmt();
        })

        ## A2 V, alpha, beta, k21
        sol.2cA2 <- RxODE({
            V <- theta[1]
            CLx <- theta[2]
            V2x <- theta[3]
            Q <- theta[4]
            ka <- theta[5]
            Kx <- CLx/V
            K12x <- Q/V
            K21 <- Q/V2x
            beta <- 0.5 * (K12x + K21 + Kx -
                           sqrt((K12x +
                                   K21 + Kx) * (K12x + K21 + Kx) - 4 * K21 *
                                  Kx))
            alpha <- K21 * Kx/beta
            C2=linCmt();
        })

        ## A3 alpha, beta, aob
        sol.2cA3 <- RxODE({
            V <- theta[1]
            CLx <- theta[2]
            V2x <- theta[3]
            Q <- theta[4]
            ka <- theta[5]
            Kx <- CLx/V
            K12x <- Q/V
            K21x <- Q/V2x
            beta <- 0.5 * (K12x + K21x + Kx -
                           sqrt((K12x +
                                   K21x + Kx) * (K12x + K21x + Kx) - 4 * K21x *
                                  Kx))
            alpha <- K21x * Kx/beta
            Ax <- (alpha - K21x)/(alpha - beta)/V
            Bx <- (beta - K21x)/(beta - alpha)/V
            aob <- Ax/Bx
            C2=linCmt();
        })

        sol.2cSS <- RxODE({
            V <- theta[1]
            CL <- theta[2]
            V2x <- theta[3]
            Q <- theta[4]
            Ka <- theta[5]
            Vss <- V+V2x
            C2=linCmt();
        })

        sol.2cT <- RxODE({
            V <- theta[1]
            CL <- theta[2]
            VT <- theta[3]
            Q <- theta[4]
            Ka <- theta[5]
            C2=linCmt();
        })



        o.2c <- ode.2c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, KA= 0.3), events=et,linLog=ll)
        s.2c <- sol.2c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, KA=0.3), events=et,linLog=ll)
        s.2cK <- sol.2cK %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, ka=0.3), events=et,linLog=ll)
        s.2cA1 <- sol.2cA1 %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, ka=0.3), events=et,linLog=ll)
        s.2cA2 <- sol.2cA2 %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, ka=0.3), events=et,linLog=ll)
        s.2cA3 <- sol.2cA3 %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, ka=0.3), events=et,linLog=ll)
        s.2cSS <- sol.2cSS %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, ka=0.3), events=et,linLog=ll)
        s.2cT <- sol.2cT %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, ka=0.3), events=et,linLog=ll)

        test_that("2 compartment oral solved models and ODEs same.", {
            expect_equal(o.2c$C2, s.2c$C2, tolerance=tol)
            expect_equal(o.2c$C2, s.2cK$C2, tolerance=tol)
            expect_equal(o.2c$C2, s.2cA1$C2, tolerance=tol)
            expect_equal(o.2c$C2, s.2cA2$C2, tolerance=tol)
            expect_equal(o.2c$C2, s.2cA3$C2, tolerance=tol)
            expect_equal(o.2c$C2, s.2cSS$C2, tolerance=tol)
            expect_equal(o.2c$C2, s.2cT$C2, tolerance=tol)
        })

        o.2c <- ode.2c.ka %>% solve(params=c(V=40, CL=1, V2=297, Q=10, KA= 0.3), events=etSs,linLog=ll)

        s.2c <- sol.2c.ka %>% solve(params=c(V=40, CL=1, V2=297, Q=10, KA=0.3), events=etSs,linLog=ll)

        test_that("2 compartment oral steady-state solved models and ODEs same.", {
            expect_equal(o.2c$C2, s.2c$C2,tolerance=tol)
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

        sol.3cK <- RxODE({
            V <- theta[1]
            CLx <- theta[2]
            V2x <- theta[3]
            Q <- theta[4]
            Q2x <- theta[5]
            V3x <- theta[6]
            K <- CLx/V
            K12 <- Q/V
            K21 <- Q/V2x
            k13 <- Q2x/V
            k31 <- Q2x/V3x
            C2=linCmt();
        })

        sol.3cA1 <- RxODE({
            Vx <- theta[1]
            CLx <- theta[2]
            V2x <- theta[3]
            Q <- theta[4]
            Q2x <- theta[5]
            V3x <- theta[6]
            Kx <- CLx/Vx
            K12x <- Q/Vx
            K21x <- Q/V2x
            K13x <- Q2x/Vx
            K31x <- Q2x/V3x
            a0 <- Kx * K21x * K31x
            a1 <- Kx * K31x + K21x * K31x + K21x * K13x +
                Kx * K21x + K31x * K12x
            a2 <- Kx + K12x + K13x + K21x + K31x
            p <- a1 - a2 * a2/3
            q <- 2 * a2 * a2 * a2/27 - a1 * a2/3 + a0
            r1 <- sqrt(-p * p * p/27)
            r2 <- 2 * r1^(1/3)
            theta <- acos(-q/(2 * r1))/3
            alpha <- -(cos(theta) * r2 - a2/3)
            beta <- -(cos(theta + 2/3 * pi) * r2 - a2/3)
            gamma <- -(cos(theta + 4/3 * pi) * r2 - a2/3)
            A <- (K21x - alpha) * (K31x - alpha)/(alpha - beta)/(alpha - gamma)/Vx
            B <- (K21x - beta) * (K31x - beta)/(beta - alpha)/(beta - gamma)/Vx
            C <- (K21x - gamma) * (K31x - gamma)/(gamma - alpha)/(gamma - beta)/Vx
            C2=linCmt();
        })

        sol.3cVp <- RxODE({
            V <- theta[1]
            CL <- theta[2]
            Vp <- theta[3]
            Q <- theta[4]
            Q2 <- theta[5]
            Vp2 <- theta[6]
            C2=linCmt();
        })

        sol.3cVt <- RxODE({
            V <- theta[1]
            CL <- theta[2]
            Vt <- theta[3]
            Q <- theta[4]
            Q2 <- theta[5]
            Vt2 <- theta[6]
            C2=linCmt();
        })

        o.3c <- ode.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et,linLog=ll)

        s.3c <- sol.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et,linLog=ll)
        s.3cK <- sol.3cK %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et,linLog=ll)
        s.3cA1 <- sol.3cA1 %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et,linLog=ll)
        s.3cVp <- sol.3cVp %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et,linLog=ll)
        s.3cVt <- sol.3cVt %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et,linLog=ll)

        test_that("3 compartment solved models and ODEs same.", {
            expect_equal(o.3c$C2, s.3c$C2, tolerance=tol)
            expect_equal(o.3c$C2, s.3cK$C2, tolerance=tol)
            expect_equal(o.3c$C2, s.3cA1$C2, tolerance=tol)
            expect_equal(o.3c$C2, s.3cVp$C2, tolerance=tol)
            expect_equal(o.3c$C2, s.3cVt$C2, tolerance=tol)
        })

        o.3c <- ode.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=etSs,linLog=ll)

        s.3c <- sol.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=etSs,linLog=ll)

        test_that("3 compartment solved models and ODEs same with steady state.", {
            expect_equal(o.3c$C2, s.3c$C2, tolerance=tol)
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

        sol.3cK <- RxODE({
            V <- theta[1]
            CLx <- theta[2]
            V2x <- theta[3]
            Q <- theta[4]
            Q2x <- theta[5]
            V3x <- theta[6]
            ka <- theta[7]
            K <- CLx/V
            K12 <- Q/V
            K21 <- Q/V2x
            k13 <- Q2x/V
            k31 <- Q2x/V3x
            C2=linCmt();
        })

        sol.3cA1 <- RxODE({
            Vx <- theta[1]
            CLx <- theta[2]
            V2x <- theta[3]
            Q <- theta[4]
            Q2x <- theta[5]
            V3x <- theta[6]
            ka <- theta[7]
            Kx <- CLx/Vx
            K12x <- Q/Vx
            K21x <- Q/V2x
            K13x <- Q2x/Vx
            K31x <- Q2x/V3x
            a0 <- Kx * K21x * K31x
            a1 <- Kx * K31x + K21x * K31x + K21x * K13x +
                Kx * K21x + K31x * K12x
            a2 <- Kx + K12x + K13x + K21x + K31x
            p <- a1 - a2 * a2/3
            q <- 2 * a2 * a2 * a2/27 - a1 * a2/3 + a0
            r1 <- sqrt(-p * p * p/27)
            r2 <- 2 * r1^(1/3)
            theta <- acos(-q/(2 * r1))/3
            alpha <- -(cos(theta) * r2 - a2/3)
            beta <- -(cos(theta + 2/3 * pi) * r2 - a2/3)
            gamma <- -(cos(theta + 4/3 * pi) * r2 - a2/3)
            A <- (K21x - alpha) * (K31x - alpha)/(alpha - beta)/(alpha - gamma)/Vx
            B <- (K21x - beta) * (K31x - beta)/(beta - alpha)/(beta - gamma)/Vx
            C <- (K21x - gamma) * (K31x - gamma)/(gamma - alpha)/(gamma - beta)/Vx
            C2=linCmt();
        })

        o.3c <- ode.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3), events=et,linLog=ll)
        s.3c <- sol.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3), events=et,linLog=ll)
        s.3cK <- sol.3cK %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3), events=et,linLog=ll)
        s.3cA1 <- sol.3cA1 %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3), events=et,linLog=ll)



        test_that("3 compartment oral solved models and ODEs same.", {
            expect_equal(o.3c$C2, s.3c$C2,tolerance=tol)
            expect_equal(o.3c$C2, s.3cK$C2,tolerance=tol)
            expect_equal(o.3c$C2, s.3cA1$C2,tolerance=tol)
        })

        o.3c <- ode.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3), events=etSs,linLog=ll)

        s.3c <- sol.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3), events=etSs,linLog=ll)

        ## Again the 4 hour strange discontinuity because ss=1
        test_that("3 compartment oral solved models and ODEs same for steady state.", {
            expect_equal(o.3c$C2, s.3c$C2, tolerance=tol)
        })

        context(sprintf("Infusion Models (%s)", ifelse(ll, "linlog", "linear")))

        et <- eventTable() %>% add.dosing(dose=3, rate=1.5, nbr.doses=6, dosing.interval=8) %>%
            add.sampling(seq(0, 48, length.out=200))

        etSs <- et() %>% et(amt=3, rate=1.5) %>%
            et(time=4,amt=3, rate=1.5, ss=1, ii=24) %>%
            et(time=8, amt=3, rate=1.5, ss=2, ii=24) %>%
            et(seq(0,24,length.out=200))

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

        o.1c <- ode.1c %>% solve(params=c(V=20, CL=25), events=et,linLog=ll)

        s.1c <- ode.1cs2 %>% solve(params=c(V=20, CL=25), events=et,linLog=ll)

        s.2c <- ode.1cs %>% solve(theta=c(20, 25), events=et,linLog=ll)

        test_that("1 compartment solved models and ODEs same.", {
            expect_equal(o.1c$C2, s.1c$C2, tolerance=tol)
            expect_equal(o.1c$C2, s.2c$C2, tolerance=tol)
        })

        o.1c <- ode.1c %>% solve(params=c(V=20, CL=10), events=etSs,linLog=ll)

        s.1c <- ode.1cs2 %>% solve(params=c(V=20, CL=10), events=etSs,linLog=ll)

        s.2c <- ode.1cs %>% solve(theta=c(20, 10), events=etSs,linLog=ll)

        test_that("1 compartment solved models and ODEs same; Steady State", {
            expect_equal(o.1c$C2, s.1c$C2,tolerance=tol)
            expect_equal(o.1c$C2, s.2c$C2,tolerance=tol)
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

        o.2c <- ode.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=et,linLog=ll)

        s.2c <- sol.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=et,linLog=ll)

        test_that("2 compartment solved models and ODEs same.", {
            expect_equal(o.2c$C2, s.2c$C2,tolerance=tol)
        })

        o.2c <- ode.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=etSs,linLog=ll)

        s.2c <- sol.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=etSs,linLog=ll)

        test_that("2 compartment steady state solved models and ODEs same.", {
            expect_equal(o.2c$C2, s.2c$C2,tolerance=tol)
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

        o.3c <- ode.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et,linLog=ll)

        s.3c <- sol.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et,linLog=ll)

        test_that("3 compartment solved models and ODEs same.", {
            expect_equal(o.3c$C2, s.3c$C2, tolerance=tol)
        })

        o.3c <- ode.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=etSs,linLog=ll)

        s.3c <- sol.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=etSs,linLog=ll)

        test_that("3 compartment steady state solved models and ODEs same.", {
            expect_equal(o.3c$C2, s.3c$C2,tolerance=tol)
        })

        context(sprintf("Infusion + Bolus (%s)", ifelse(ll, "linlog", "linear")))

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

        o.1c <- ode.1c %>% solve(params=c(V=20, CL=25), events=et,linLog=ll)

        s.1c <- ode.1cs2 %>% solve(params=c(V=20, CL=25), events=et,linLog=ll)

        s.2c <- ode.1cs %>% solve(theta=c(20, 25), events=et,linLog=ll)

        test_that("1 compartment solved models and ODEs same.", {
            expect_equal(o.1c$C2, s.1c$C2, tolerance=tol)
            expect_equal(o.1c$C2, s.2c$C2, tolerance=tol)
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

        o.2c <- ode.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=et,linLog=ll)

        s.2c <- sol.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=et,linLog=ll)

        test_that("2 compartment solved models and ODEs same.", {
            expect_equal(o.2c$C2, s.2c$C2, tolerance=tol)
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

        o.3c <- ode.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et,linLog=ll)

        s.3c <- sol.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et,linLog=ll)

        test_that("3 compartment solved models and ODEs same.", {
            expect_equal(o.3c$C2, s.3c$C2, tolerance=tol)
        })

        context(sprintf("Oral + Infusion + Bolus Models (%s)", ifelse(ll, "linlog", "linear")))

        et <- eventTable() %>% add.dosing(dose=3, rate=1.5, nbr.doses=3, dosing.interval=16,cmt=2) %>%
            add.dosing(dose=1.5, nbr.doses=3, dosing.interval=16,cmt=2) %>%
            add.dosing(dose=1.5, nbr.doses=3, dosing.interval=16,cmt=1,start.time=8) %>%
            add.sampling(seq(0, 48, length.out=200))

        ode.1c.ka <- RxODE({
            C2 = center/V;
            d / dt(depot) = -KA * depot
            d/dt(center) = KA * depot - CL*C2
        })

        sol.1c.ka <- RxODE({
            C2 = linCmt(V, CL, KA);
        })

        o.1c <- ode.1c.ka %>% solve(params=c(V=20, CL=25, KA=2), events=et,linLog=ll)
        s.1c <- sol.1c.ka %>% solve(params=c(V=20, CL=25, KA=2), events=et,linLog=ll)

        test_that("1 compartment solved models and ODEs same for mixed oral, iv and infusion.", {
            expect_equal(o.1c$C2, s.1c$C2, tolerance=tol)
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

        o.2c <- ode.2c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, KA= 0.3), events=et,linLog=ll)

        s.2c <- sol.2c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, KA=0.3), events=et,linLog=ll)

        test_that("2 compartment solved models and ODEs same for mixed oral, iv and infusion.", {
            expect_equal(o.2c$C2, s.2c$C2, tolerance=tol)
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

        o.3c <- ode.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3), events=et,linLog=ll)

        s.3c <- sol.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3), events=et,linLog=ll)

        test_that("3 compartment solved models and ODEs same for mixed oral, iv and infusion.", {
            expect_equal(o.3c$C2, s.3c$C2,tolerance=tol)
        })

        context(sprintf("Modeled bio-availability (%s)", ifelse(ll, "linlog", "linear")))

        et <- eventTable() %>% add.dosing(dose=3, rate=1.5, nbr.doses=3, dosing.interval=16,cmt=2) %>%
            add.dosing(dose=1.5, nbr.doses=3, dosing.interval=16,cmt=2) %>%
            add.dosing(dose=1.5, nbr.doses=3, dosing.interval=16,cmt=1,start.time=8) %>%
            add.sampling(seq(0, 48, length.out=200))

        ode.1c.ka <- RxODE({
            C2 = center/V;
            d / dt(depot) = -KA * depot
            d/dt(center) = KA * depot - CL*C2
            f(depot) = fDepot
            f(center) = fCenter
        })

        sol.1c.ka <- RxODE({
            C2 = linCmt(V, CL, KA);
            f(depot) = fDepot
            f(central) = fCenter
        })

        for (fd in c(0.5,1,2)){
            for (fc in c(0.5,1,2)){
                o.1c <- ode.1c.ka %>% solve(params=c(V=20, CL=25, KA=2,fDepot=fd,fCenter=fc), events=et,linLog=ll)
                s.1c <- sol.1c.ka %>% solve(params=c(V=20, CL=25, KA=2, fDepot=fd, fCenter=fc), events=et,linLog=ll)
                test_that(sprintf("1 compartment solved models and ODEs same for mixed oral, iv and infusion + Fd=%f,Fc=%f", fd,fc), {
                    expect_equal(o.1c$C2, s.1c$C2, tolerance=tol)
                })
            }
        }

        ode.2c.ka <- RxODE({
            C2 = centr/V;
            C3 = peri/V2;
            d/dt(depot) =-KA*depot;
            d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
            d/dt(peri)  =                    Q*C2 - Q*C3;
            f(depot) = fDepot
            f(centr) = fCenter
            ## FIXME:
            ## f(central) should throw an error
        })

        sol.2c.ka <- RxODE({
            C2=linCmt(V, CL, V2, Q, KA);
            f(depot) = fDepot
            f(central) = fCenter
        })

        for (fd in c(0.5,1,2)){
            for (fc in c(0.5,1,2)){
                o.2c <- ode.2c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, KA= 0.3, fDepot=fd, fCenter=fc), events=et,linLog=ll)
                s.2c <- sol.2c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, KA=0.3, fDepot=fd, fCenter=fc), events=et,linLog=ll)
                test_that(sprintf("2 compartment solved models and ODEs same for mixed oral, iv and infusion + Fd=%f,Fc=%f", fd,fc), {
                    expect_equal(o.2c$C2, s.2c$C2,tolerance=tol)
                })
            }
        }

        ode.3c.ka <- RxODE({
            C2 = centr/V;
            C3 = peri/V2;
            C4 = peri2/V3
            d/dt(depot) =-KA*depot;
            d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3  - Q2*C2 + Q2*C4;
            d/dt(peri)  = Q*C2 - Q*C3;
            d/dt(peri2) = Q2 * C2 - Q2 * C4
            f(depot) = fDepot
            f(centr) = fCenter
        })

        sol.3c.ka <- RxODE({
            ## double solvedC(double t, int parameterization, int cmt, unsigned int col, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8);
            C2=linCmt(V, CL, V2, Q, Q2, V3, KA);
            f(depot) = fDepot
            f(central) = fCenter
        })

        for (fd in c(0.5,1,2)){
            for (fc in c(0.5,1,2)){
                o.3c <- ode.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7,
                                                     V3=400, KA=0.3, fDepot=fd, fCenter=fc), events=et,linLog=ll)
                s.3c <- sol.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3,
                                                     fDepot=fd, fCenter=fc), events=et,linLog=ll)
                test_that(sprintf("3 compartment solved models and ODEs same for mixed oral, iv and infusion + Fd=%f,Fc=%f", fd,fc), {
                    expect_equal(o.3c$C2, s.3c$C2,tolerance=tol)
                })
            }
        }

        context(sprintf("Modeled lag time (%s)", ifelse(ll, "linlog", "linear")))

        et <- eventTable() %>% add.dosing(dose=3, rate=1.5, nbr.doses=3, dosing.interval=16,cmt=2) %>%
            add.dosing(dose=1.5, nbr.doses=3, dosing.interval=16,cmt=2) %>%
            add.dosing(dose=1.5, nbr.doses=3, dosing.interval=16,cmt=1,start.time=8) %>%
            add.sampling(seq(0, 48, length.out=200))

        ode.1c.ka <- RxODE({
            C2 = center/V;
            d / dt(depot) = -KA * depot
            d/dt(center) = KA * depot - CL*C2
            alag(depot) = lagDepot
            alag(center) = lagCenter
        })

        sol.1c.ka <- RxODE({
            C2 = linCmt(V, CL, KA);
            alag(depot) = lagDepot
            alag(central) = lagCenter
        })

        for (fd in c(1,2,10)){
            for (fc in c(1,2,10)){
                o.1c <- ode.1c.ka %>% solve(params=c(V=20, CL=25, KA=2,lagDepot=fd,lagCenter=fc), events=et,linLog=ll)
                s.1c <- sol.1c.ka %>% solve(params=c(V=20, CL=25, KA=2, lagDepot=fd, lagCenter=fc), events=et,linLog=ll)
                test_that(sprintf("1 compartment solved models and ODEs same for mixed oral, iv and infusion + Fd=%f,Fc=%f", fd,fc), {
                    expect_equal(o.1c$C2, s.1c$C2, tolerance=tol)
                })
            }
        }

        ode.2c.ka <- RxODE({
            C2 = centr/V;
            C3 = peri/V2;
            d/dt(depot) =-KA*depot;
            d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
            d/dt(peri)  =                    Q*C2 - Q*C3;
            alag(depot) = lagDepot
            alag(centr) = lagCenter
        })

        sol.2c.ka <- RxODE({
            C2=linCmt(V, CL, V2, Q, KA);
            alag(depot) = lagDepot
            alag(central) = lagCenter
        })

        for (fd in c(1,2,10)){
            for (fc in c(1,2,10)){
                o.2c <- ode.2c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, KA= 0.3, lagDepot=fd, lagCenter=fc), events=et,linLog=ll)
                s.2c <- sol.2c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, KA=0.3, lagDepot=fd, lagCenter=fc), events=et,linLog=ll)
                test_that(sprintf("2 compartment solved models and ODEs same for mixed oral, iv and infusion + Fd=%f,Fc=%f", fd,fc), {
                    expect_equal(o.2c$C2, s.2c$C2,tolerance=tol)
                })
            }
        }

        ode.3c.ka <- RxODE({
            C2 = centr/V;
            C3 = peri/V2;
            C4 = peri2/V3
            d/dt(depot) =-KA*depot;
            d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3  - Q2*C2 + Q2*C4;
            d/dt(peri)  = Q*C2 - Q*C3;
            d/dt(peri2) = Q2 * C2 - Q2 * C4
            alag(depot) = lagDepot
            alag(centr) = lagCenter
        })

        sol.3c.ka <- RxODE({
            ## double solvedC(double t, int parameterization, int cmt, unsigned int col, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8);
            C2=linCmt(V, CL, V2, Q, Q2, V3, KA);
            alag(depot) = lagDepot
            alag(central) = lagCenter
        })

        for (fd in c(1,2,10)){
            for (fc in c(1,2,10)){
                o.3c <- ode.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7,
                                                     V3=400, KA=0.3, lagDepot=fd, lagCenter=fc), events=et,linLog=ll)
                s.3c <- sol.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3,
                                                     lagDepot=fd, lagCenter=fc), events=et,linLog=ll)
                test_that(sprintf("3 compartment solved models and ODEs same for mixed oral, iv and infusion + Fd=%f,Fc=%f", fd,fc), {
                    expect_equal(o.3c$C2, s.3c$C2,tolerance=tol)
                })
            }
        }

        context(sprintf("Modeled rate (%s)", ifelse(ll, "linlog", "linear")))

        ode.1c <- RxODE({
            C2 = center/V;
            d/dt(center) = - CL*C2
            rate(center) = rt
        })

        sol.1c <- RxODE({
            C2 = linCmt(CL, V);
            rate(central) = rt
        })

        et <- eventTable() %>% add.dosing(dose=3, rate=-1, nbr.doses=3, cmt=1,dosing.interval=12) %>%
            add.sampling(seq(0, 36, length.out=200))

        for (rt in seq(0.5, 1, 1.5)){
            o.1c <- ode.1c %>% solve(params=c(V=20, CL=25,rt=rt), events=et,linLog=ll)
            s.1c <- sol.1c %>% solve(params=c(V=20, CL=25,rt=rt), events=et,linLog=ll)
            test_that(sprintf("1 compartment solved models and ODEs same for rate-modeled infusion: %s", rt), {
                expect_equal(o.1c$C2, s.1c$C2,tolerance=tol)
            })
        }

        ode.2c <- RxODE({
            C2 = centr/V;
            C3 = peri/V2;
            d/dt(centr) = - CL*C2 - Q*C2 + Q*C3;
            d/dt(peri)  = Q*C2 - Q*C3;
            rate(centr) = rt
        })

        sol.2c <- RxODE({
            C2=linCmt(V, CL, V2, Q);
            rate(central) = rt
        })

        for (rt in seq(0.5, 1, 1.5)){
            o.2c <- ode.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, rt=rt), events=et,linLog=ll)
            s.2c <- sol.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, rt=rt), events=et,linLog=ll)
            test_that(sprintf("2 compartment solved models and ODEs same for rate-modeled infusion: %s", rt), {
                expect_equal(o.2c$C2, s.2c$C2,tolerance=tol)
            })
        }

        ode.3c <- RxODE({
            C2 = centr/V;
            C3 = peri/V2;
            C4 = peri2 / V3
            d/dt(centr) = - CL*C2 - Q*C2 + Q*C3  - Q2*C2 + Q2*C4;
            d/dt(peri)  = Q*C2 - Q*C3;
            d / dt(peri2) = Q2 * C2 - Q2 * C4
            rate(centr) = rt
        })

        sol.3c <- RxODE({
            ## double solvedC(double t, int parameterization, int cmt, unsigned int col, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8);
            C2=linCmt(V, CL, V2, Q, Q2, V3);
            rate(central) = rt
        })

        for (rt in seq(0.5, 1, 1.5)){
            s.3c <- sol.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, rt=rt), events=et,linLog=ll)
            o.3c <- ode.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, rt=rt), events=et,linLog=ll)
            test_that(sprintf("3 compartment solved models and ODEs same for rate-modeled infusion: %s", rt), {
                expect_equal(o.3c$C2, s.3c$C2,tolerance=tol)
            })
        }

        context(sprintf("Modeled duration (%s)", ifelse(ll, "linlog", "linear")))

        ode.1c <- RxODE({
            C2 = center/V;
            d/dt(center) = - CL*C2
            dur(center) = dr
        })

        sol.1c <- RxODE({
            C2 = linCmt(CL, V);
            dur(central) = dr
        })

        et <- eventTable() %>% add.dosing(dose=3, rate=-2, nbr.doses=3, cmt=1,dosing.interval=12) %>%
            add.sampling(seq(0, 36, length.out=200))

        for (dur in seq(0.5, 1, 1.5)){
            o.1c <- ode.1c %>% solve(params=c(V=20, CL=25,dr=dur), events=et,linLog=ll)
            s.1c <- sol.1c %>% solve(params=c(V=20, CL=25,dr=dur), events=et,linLog=ll)
            test_that(sprintf("1 compartment solved models and ODEs same for dur-modeled infusion: %s", dur), {
                expect_equal(o.1c$C2, s.1c$C2,tolerance=tol)
            })
        }

        ode.2c <- RxODE({
            C2 = centr/V;
            C3 = peri/V2;
            d/dt(centr) = - CL*C2 - Q*C2 + Q*C3;
            d/dt(peri)  = Q*C2 - Q*C3;
            dur(centr) = dr
        })

        sol.2c <- RxODE({
            C2=linCmt(V, CL, V2, Q);
            dur(central) = dr
        })

        for (dur in seq(0.5, 1, 1.5)){
            o.2c <- ode.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, dr=dur), events=et,linLog=ll)
            s.2c <- sol.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, dr=dur), events=et,linLog=ll)
            test_that(sprintf("2 compartment solved models and ODEs same for dur-modeled infusion: %s", dur), {
                expect_equal(o.2c$C2, s.2c$C2,tolerance=tol)
            })
        }

        ode.3c <- RxODE({
            C2 = centr/V;
            C3 = peri/V2;
            C4 = peri2 / V3
            d/dt(centr) = - CL*C2 - Q*C2 + Q*C3  - Q2*C2 + Q2*C4;
            d/dt(peri)  = Q*C2 - Q*C3;
            d / dt(peri2) = Q2 * C2 - Q2 * C4
            dur(centr) = dr
        })

        sol.3c <- RxODE({
            ## double solvedC(double t, int parameterization, int cmt, unsigned int col, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8);
            C2=linCmt(V, CL, V2, Q, Q2, V3);
            dur(central) = dr
        })

        for (dur in seq(0.5, 1, 1.5)){
            o.3c <- ode.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, dr=dur), events=et,linLog=ll)
            s.3c <- sol.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, dr=dur), events=et,linLog=ll)
            test_that(sprintf("3 compartment solved models and ODEs same for dur-modeled infusion: %s", dur), {
                expect_equal(o.3c$C2, s.3c$C2,tolerance=tol)
            })
        }

        test_that("central should throw error",{
            expect_error(RxODE({
                C2 = centr/V;
                C3 = peri/V2;
                d/dt(depot) =-KA*depot;
                d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
                d/dt(peri)  =                    Q*C2 - Q*C3;
                alag(depot) = lagDepot
                alag(centr) = lagCenter
                alag(central) = matt
            }))
        })

        test_that("depot should throw error",{
            expect_error(RxODE({
                C2 = centr/V;
                C3 = peri/V2;
                d/dt(dep) =-KA*dep;
                d/dt(centr) = KA*dep - CL*C2 - Q*C2 + Q*C3;
                d/dt(peri)  =                    Q*C2 - Q*C3;
                alag(dep) = lagDepot
                alag(centr) = lagCenter
                alag(depot) = matt
            }))
        })

        ode.1cs2 <- RxODE({
            C2 = linCmt(CL, V);
            mtime(t1) = mt1
            mtime(t2) = mt2
        })

        et <- eventTable() %>%
            add.dosing(dose=3, nbr.doses=6, dosing.interval=8) %>%
            add.sampling(0:48)

        s.1c <- ode.1cs2 %>% solve(params=c(V=20, CL=25,mt1=0.5, mt2=1.75),
                                   events=et,linLog=ll)

        test_that("mtime with solved systems work",{
            expect_equal(s.1c$time[1:4], c(0, 0.5, 1, 1.75))
        })

    }

    test_that("double linCmt has error",
              expect_error(RxODE({
                  C2 = linCmt(CL, V);
                  C2 = linCmt(CL, V);
              })))

    test_that("Steady state IV infusion", {

        ev <- et(amt=0, ss=1,rate=10000/8)

        ode.1c <- RxODE({
            C2 = center/V;
            d/dt(center) = - CL*C2
        })

        ode.1cs2 <- RxODE({
            C2 = linCmt(CL, V);
        })

        o.1c <- ode.1c %>% solve(params=c(V=20, CL=25), events=ev)

        s.1c <- ode.1cs2 %>% solve(params=c(V=20, CL=25), events=ev)

        expect_equal(o.1c$C2, s.1c$C2, tolerance=tol)

        ode.2c <- RxODE({
            C2 = centr/V;
            C3 = peri/V2;
            d/dt(centr) = - CL*C2 - Q*C2 + Q*C3;
            d/dt(peri)  = Q*C2 - Q*C3;
        })

        sol.2c <- RxODE({
            C2=linCmt(V, CL, V2, Q1);
        })

        o.2c <- ode.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=ev)

        s.2c <- sol.2c %>% solve(params=c(V=40, CL=18, V2=297, Q1=10), events=ev)

        expect_equal(o.2c$C2, s.2c$C2, tolerance=tol)

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

        o.3c <- ode.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=ev)

        s.3c <- sol.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=ev)

        expect_equal(o.3c$C2, s.3c$C2, tolerance=tol)

    })

}, silent=TRUE)
