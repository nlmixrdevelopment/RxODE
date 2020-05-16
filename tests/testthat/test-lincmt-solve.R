if (FALSE){
rxPermissive({

    et <- eventTable() %>% add.dosing(dose=3, rate=1.5, nbr.doses=6, dosing.interval=8) %>%
    add.dosing(dose=1.5, nbr.doses=6, dosing.interval=8) %>%
    add.sampling(seq(0, 48, length.out=200))

    ode.1cs <- RxODE({
      V <- theta[1];
      CL <- theta[2];
      C2 = linCmt();
    })

    s.2c <- ode.1cs %>% solve(theta=c(20, 25), events=et)


dfadvan <- structure(list(ID = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L), TIME = c(0L, 1L, 2L,
3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 12L, 13L, 14L, 15L,
16L, 17L, 18L, 19L, 20L, 21L, 22L, 23L, 24L, 0L, 1L, 2L, 3L,
4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 12L, 13L, 14L, 15L, 16L,
17L, 18L, 19L, 20L, 21L, 22L, 23L, 24L), AMT = c(100L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 100L, 0L, 0L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 100L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 0L, 100L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 0L, 0L, 0L), MDV = c(1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L,
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L), CLCR = c(120L, 120L,
120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L,
120L, 120L, 120L, 30L, 30L, 30L, 30L, 30L, 30L, 30L, 30L, 30L,
30L, 30L, 30L, 30L, 30L, 30L, 120L, 120L, 120L, 120L, 120L, 120L,
120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L,
120L, 120L, 120L, 120L)), row.names = c(NA, -52L), class = "data.frame")

dfadvanR <-structure(list(ID = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L,
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L
), TIME = c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 20.5, 20.5,
22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52,
54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84,
86, 88, 90, 92, 94, 96, 98, 100, 0, 2, 4, 6, 8, 10, 12, 14, 16,
18, 20, 20.5, 20.5, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42,
44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74,
76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100), AMT = c(100L,
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 100L, 0L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 100L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 0L, 100L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L), MDV = c(1L,
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L), RATE = c(2L, 0L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 2L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 2L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 0L), CLCR = c(120L, 120L, 120L, 120L, 120L,
120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L,
120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L,
120L, 30L, 30L, 30L, 30L, 30L, 30L, 30L, 30L, 30L, 30L, 30L,
30L, 30L, 30L, 30L, 30L, 30L, 30L, 30L, 30L, 30L, 30L, 30L, 30L,
30L, 60L, 60L, 60L, 60L, 60L, 60L, 60L, 60L, 60L, 60L, 60L, 60L,
60L, 60L, 60L, 60L, 60L, 60L, 60L, 60L, 60L, 60L, 120L, 120L,
120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L,
120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L, 120L,
120L, 120L, 120L, 120L, 120L, 120L, 120L), DV = c(NA, NA, NA,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
NA, NA, NA, NA, NA, NA, NA)), class = "data.frame", row.names = c(NA,
-106L))

    ## context("time-varying linCmt advan tests")

    ## test_that("time varing advan", {

        ## advan


        ## mod1iv <- RxODE({
        ##     CLpop ~ 2       # clearance
        ##     Vpop ~ 10      # central volume of distribution
        ##     CL <- CLpop*(CLCR/100)
        ##     V <- Vpop
        ##     cp = linCmt()
        ## })

        ## tmp2 <- rxSolve(mod1iv, events=dfadvan, advanLinCmt=TRUE)

        ## plot(tmp2, cp) + geom_vline(xintercept=14) +
        ##     geom_vline(xintercept=4)

        ## mod2iv <- RxODE({
        ##     CLpop <- 2       # clearance
        ##     V1pop <- 10      # central volume of distribution
        ##     Qpop  <- 1       # inter-compartmental clearance
        ##     V2pop <- 25      # peripheral volume of distribution
        ##     CL <- CLpop*(CLCR/100)           #creatinine clearance (CLCR) added as a covariate on CL
        ##     V1 <- V1pop
        ##     Q  <- Qpop
        ##     V2 <- V2pop
        ##     cp <- linCmt()
        ## })

        ## tmp2 <- rxSolve(mod2iv, events=dfadvan, advanLinCmt=TRUE)

        ## plot(tmp2, cp) + geom_vline(xintercept=14) +
        ##     geom_vline(xintercept=4)

        ## mod3iv <- RxODE({
        ##     CLpop  <- 2          # clearance
        ##     V1pop  <- 10         # central volume of distribution
        ##     Q12pop <- 0.5        # inter-compartmental clearance (1)
        ##     V2pop  <- 30         # peripheral volume of distribution (1)
        ##     Q13pop <- 0.3        # inter-compartmental clearance (2)
        ##     V3pop  <- 40         # peripheral volume of distribution (2)
        ##     CL  <- CLpop*(CLCR/100)    #creatinine clearance (CLCR) added as a covariate on CL
        ##     V1  <- V1pop
        ##     Q <- Q12pop
        ##     V2  <- V2pop
        ##     Q2 <- Q13pop
        ##     V3  <- V3pop
        ##     cp <- linCmt()
        ## })

        ## tmp3 <- rxSolve(mod3iv, events=dfadvan, advanLinCmt=TRUE)

        ## plot(tmp2, cp) + geom_vline(xintercept=14) +
        ##     geom_vline(xintercept=4)

        ## tmp1 <- rxSolve(mod1iv, events=dfadvanR, advanLinCmt=FALSE)

        ## plot(tmp1, cp) + geom_vline(xintercept=52) +
        ##     geom_vline(xintercept=40) + facet_wrap( ~ id)

    ## })

    ## stop()

    tol  <- 5e-6 ## Current difference for all equations
    type <- 1

    for (type in 1:2){

      .txt <- switch(type, "linear", "sensitivitiy");
      sens <- switch(type, FALSE, TRUE);
      context(sprintf("Test the solved equations (%s)", .txt))

      et <- eventTable() %>% add.dosing(dose=3, nbr.doses=6, dosing.interval=8) %>%
        add.sampling(seq(0, 48, length.out=200))

      ode.1c <- RxODE({
        C2 = center/V;
        d/dt(center) = - CL*C2
      })

      test_that("ode model gives extraCmt=0",{
        expect_equal(rxModelVars(ode.1c)$extraCmt,0L);
      })

      goodP <- function(model, cmt=1L, ka=0L){
        test_that(sprintf("model '%s' parses to cmt=%d, ka=%d", as.character(substitute(model)), cmt, ka),{
          .flags <- rxModelVars(model)$flags
          expect_equal(setNames(.flags["ncmt"], NULL), cmt)
          expect_equal(setNames(.flags["ka"], NULL), ka)
        })
      }

      ## Solved systems can check the variables in the RxODE statement
      ## to figure out what type of solved system is being requested
      ode.1cs <- RxODE({
        V <- theta[1];
        CL <- theta[2];
        C2 = linCmt();
      }, linCmtSens=sens)

      goodP(ode.1cs)

      ode.2cK <- RxODE({
        V <- theta[1];
        CLx <- theta[2];
        K <- CLx/V
        C2 = linCmt();
      }, linCmtSens=sens)

      goodP(ode.2cK)

      ode.2cA1 <- RxODE({
        V <- theta[1];
        CLx <- theta[2];
        alpha <- CLx/V
        C2 = linCmt();
      }, linCmtSens=sens)

      goodP(ode.2cA1)

      ode.2cA2 <- RxODE({
        A <- 1/theta[1];
        CLx <- theta[2];
        alpha <- CLx*A
        C2 = linCmt();
      }, linCmtSens=sens)

      goodP(ode.2cA2)

      ## Instead of specifying parameters in the solved system, you can
      ## specify them in the linCmt variable.
      ode.1cs2 <- RxODE({
        C2 = linCmt(CL, V);
      }, linCmtSens=sens)

      goodP(ode.1cs2)

      test_that("linear compartment model gives extraCmt=1",{
        expect_equal(rxModelVars(ode.1cs2)$extraCmt,1L);
      })

      ## The solved systems can be mixed with ODE solving routines (to
      ## speed them up a bit...?)

      o.1c <- ode.1c %>% solve(params=c(V=20, CL=25), events=et)
      
      s.1c <- ode.1cs2 %>% solve(params=c(V=20, CL=25), events=et)

      s.2c <- ode.1cs %>% solve(theta=c(20, 25), events=et)
      s.2cK <- ode.2cK %>% solve(theta=c(20, 25), events=et)
      s.2cA1 <- ode.2cA1 %>% solve(theta=c(20, 25), events=et)
      s.2cA2 <- ode.2cA2 %>% solve(theta=c(20, 25), events=et)

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

      o.1c <- ode.1c %>% solve(params=c(V=20, CL=1), events=etSs)

      s.1c <- ode.1cs2 %>% solve(params=c(V=20, CL=1), events=etSs)

      s.2c <- ode.1cs %>% solve(theta=c(20, 1), events=etSs)

      test_that("1 compartment steady-state solved models and ODEs same.", {
        expect_equal(o.1c$C2, s.1c$C2, tolerance=tol)
        expect_equal(o.1c$C2, s.2c$C2, tolerance=tol)
      })

      ode.1c.ka <- RxODE({
        C2 = center/V;
        d/dt(depot) = -KA * depot
        d/dt(center) = KA * depot - CL*C2
      })

      sol.1c.ka <- RxODE({
        C2 = linCmt(V, CL, KA);
      }, linCmtSens=sens)

      goodP(sol.1c.ka, ka=1L)

      ode.2cK <- RxODE({
        V <- theta[1];
        CLx <- theta[2];
        Ka <- theta[3]
        K <- CLx/V
        C2 = linCmt();
      }, linCmtSens=sens)

      goodP(ode.2cK, ka=1L)

      ode.2cA1 <- RxODE({
        V <- theta[1];
        CLx <- theta[2];
        Ka <- theta[3]
        alpha <- CLx/V
        C2 = linCmt();
      }, linCmtSens=sens)

      goodP(ode.2cA1, ka=1L)

      ode.2cA2 <- RxODE({
        A <- 1/theta[1];
        CLx <- theta[2];
        Ka <- theta[3];
        alpha <- CLx*A
        C2 = linCmt();
      }, linCmtSens=sens)

      goodP(ode.2cA2, ka=1L)

      test_that("linear oral model gives extraCmt=2",{
        expect_equal(rxModelVars(sol.1c.ka)$extraCmt,2L);
      })

      o.1c <- ode.1c.ka %>% solve(params=c(V=20, CL=25, KA=2), events=et)

      s.1c <- sol.1c.ka %>% solve(params=c(V=20, CL=25, KA=2), events=et)

      s.2cK <- ode.2cK %>% solve(theta=c(20, 25, KA=2), events=et)
      s.2cA1 <- ode.2cA1 %>% solve(theta=c(20, 25, KA=2), events=et)
      s.2cA2 <- ode.2cA2 %>% solve(theta=c(20, 25, KA=2), events=et)

      test_that("1 compartment oral solved models and ODEs same.", {
        expect_equal(o.1c$C2, s.1c$C2, tolerance=tol)
        expect_equal(o.1c$C2, s.2cK$C2, tolerance=tol)
        expect_equal(o.1c$C2, s.2cA1$C2, tolerance=tol)
        expect_equal(o.1c$C2, s.2cA2$C2, tolerance=tol)
      })

      ## Note the strange-looking dip at 4 hours.  This is because ss=1 resets the system first.
      o.1c <- ode.1c.ka %>% solve(params=c(V=20, CL=2, KA=2), events=etSs)

      s.1c <- sol.1c.ka %>% solve(params=c(V=20, CL=2, KA=2), events=etSs)

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
      }, linCmtSens=sens)

      goodP(sol.2c, cmt=2L)

      sol.2cK <- RxODE({
        V <- theta[1]
        CLx <- theta[2]
        V2x <- theta[3]
        Q <- theta[4]
        K <- CLx/V
        K12 <- Q/V
        K21 <- Q/V2x
        C2=linCmt();
      }, linCmtSens=sens)

      goodP(sol.2cK, cmt=2L)

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
      }, linCmtSens=sens)

      goodP(sol.2cA1, cmt=2L)

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
      }, linCmtSens=sens)

      goodP(sol.2cA2, cmt=2L)

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
      }, linCmtSens=sens)

      goodP(sol.2cA3, cmt=2L)
      
      o.2c <- ode.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=et)
      
      s.2c <- sol.2c %>% solve(params=c(V=40, CL=18, V2=297, Q1=10), events=et)
      
      s.2cK <- sol.2cK %>% solve(theta=c(V=40, CL=18, V2=297, Q=10), events=et)
      
      s.2cA1 <- sol.2cA1 %>% solve(theta=c(V=40, CL=18, V2=297, Q=10), events=et)
      
      s.2cA2 <- sol.2cA2 %>% solve(theta=c(V=40, CL=18, V2=297, Q=10), events=et)
      s.2cA3 <- sol.2cA3 %>% solve(theta=c(V=40, CL=18, V2=297, Q=10), events=et)
      
      test_that("2 compartment solved models and ODEs same.", {
        expect_equal(s.2cK$C2, s.2c$C2, tolerance=tol)
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
      }, linCmtSens=sens)

      goodP(sol.2c.ka, cmt=2, ka=1)

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
      }, linCmtSens=sens)

      goodP(sol.2cK, cmt=2, ka=1)

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
      }, linCmtSens=sens)

      goodP(sol.2cA1, cmt=2, ka=1)

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
      }, linCmtSens=sens)

      goodP(sol.2cA2, cmt=2, ka=1)

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
      }, linCmtSens=sens)

      goodP(sol.2cA3, cmt=2, ka=1)

      sol.2cSS <- RxODE({
        V <- theta[1]
        CL <- theta[2]
        V2x <- theta[3]
        Q <- theta[4]
        Ka <- theta[5]
        Vss <- V+V2x
        C2=linCmt();
      }, linCmtSens=sens)

      goodP(sol.2cSS, cmt=2, ka=1)

      sol.2cT <- RxODE({
        V <- theta[1]
        CL <- theta[2]
        VT <- theta[3]
        Q <- theta[4]
        Ka <- theta[5]
        C2=linCmt();
      }, linCmtSens=sens)

      goodP(sol.2cT, cmt=2, ka=1)

      o.2c <- ode.2c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, KA= 0.3), events=et)
      s.2c <- sol.2c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, KA=0.3), events=et)
      s.2cK <- sol.2cK %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, ka=0.3), events=et)
      s.2cA1 <- sol.2cA1 %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, ka=0.3), events=et)
      s.2cA2 <- sol.2cA2 %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, ka=0.3), events=et)
      s.2cA3 <- sol.2cA3 %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, ka=0.3), events=et)
      s.2cSS <- sol.2cSS %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, ka=0.3), events=et)
      s.2cT <- sol.2cT %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, ka=0.3), events=et)

      test_that("2 compartment oral solved models and ODEs same.", {
        expect_equal(o.2c$C2, s.2c$C2, tolerance=tol)
        expect_equal(o.2c$C2, s.2cK$C2, tolerance=tol)
        expect_equal(o.2c$C2, s.2cA1$C2, tolerance=tol)
        expect_equal(o.2c$C2, s.2cA2$C2, tolerance=tol)
        expect_equal(o.2c$C2, s.2cA3$C2, tolerance=tol)
        expect_equal(o.2c$C2, s.2cSS$C2, tolerance=tol)
        expect_equal(o.2c$C2, s.2cT$C2, tolerance=tol)
      })

      o.2c <- ode.2c.ka %>% solve(params=c(V=40, CL=1, V2=297, Q=10, KA= 0.3), events=etSs)

      s.2c <- sol.2c.ka %>% solve(params=c(V=40, CL=1, V2=297, Q=10, KA=0.3), events=etSs)

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
      }, linCmtSens=sens)

      goodP(sol.3c, 3)

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
      }, linCmtSens=sens)

      goodP(sol.3cK, 3)

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
      }, linCmtSens=sens)

      goodP(sol.3cA1, 3)

      sol.3cVp <- RxODE({
        V <- theta[1]
        CL <- theta[2]
        Vp <- theta[3]
        Q <- theta[4]
        Q2 <- theta[5]
        Vp2 <- theta[6]
        C2=linCmt();
      }, linCmtSens=sens)

      goodP(sol.3cVp, 3)

      sol.3cVt <- RxODE({
        V <- theta[1]
        CL <- theta[2]
        Vt <- theta[3]
        Q <- theta[4]
        Q2 <- theta[5]
        Vt2 <- theta[6]
        C2=linCmt();
      }, linCmtSens=sens)

      goodP(sol.3cVt, 3)

      o.3c <- ode.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et)

      s.3c <- sol.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et)

      s.3cK <- sol.3cK %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et)
      s.3cA1 <- sol.3cA1 %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et)
      s.3cVp <- sol.3cVp %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et)
      s.3cVt <- sol.3cVt %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et)

      test_that("3 compartment solved models and ODEs same.", {
        expect_equal(o.3c$C2, s.3c$C2, tolerance=tol)
        expect_equal(o.3c$C2, s.3cK$C2, tolerance=tol)
        expect_equal(o.3c$C2, s.3cA1$C2, tolerance=tol)
        expect_equal(o.3c$C2, s.3cVp$C2, tolerance=tol)
        expect_equal(o.3c$C2, s.3cVt$C2, tolerance=tol)
      })

      o.3c <- ode.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=etSs)

      s.3c <- sol.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=etSs)

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
        }, linCmtSens=sens)

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
        }, linCmtSens=sens)

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
        }, linCmtSens=sens)

        o.3c <- ode.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3), events=et)
        s.3c <- sol.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3), events=et)
        s.3cK <- sol.3cK %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3), events=et)
        s.3cA1 <- sol.3cA1 %>% solve(theta=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3), events=et)

        test_that("3 compartment oral solved models and ODEs same.", {
            expect_equal(o.3c$C2, s.3c$C2,tolerance=tol)
            expect_equal(o.3c$C2, s.3cK$C2,tolerance=tol)
            expect_equal(o.3c$C2, s.3cA1$C2,tolerance=tol)
        })

        o.3c <- ode.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3), events=etSs)

        s.3c <- sol.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3), events=etSs)

        ## Again the 4 hour strange discontinuity because ss=1
        test_that("3 compartment oral solved models and ODEs same for steady state.", {
            expect_equal(o.3c$C2, s.3c$C2, tolerance=tol)
        })

        context(sprintf("Infusion Models (%s)", .txt))

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
        }, linCmtSens=sens)

        ## Instead of specifying parameters in the solved system, you can
        ## specify them in the linCmt variable.
        ode.1cs2 <- RxODE({
            C2 = linCmt(CL, V);
        }, linCmtSens=sens)

        ## The solved systems can be mixed with ODE solving routines (to
        ## speed them up a bit...?)

        o.1c <- ode.1c %>% solve(params=c(V=20, CL=25), events=et)

        s.1c <- ode.1cs2 %>% solve(params=c(V=20, CL=25), events=et)

        s.2c <- ode.1cs %>% solve(theta=c(20, 25), events=et)

        test_that("1 compartment solved models and ODEs same.", {
            expect_equal(o.1c$C2, s.1c$C2, tolerance=tol)
            expect_equal(o.1c$C2, s.2c$C2, tolerance=tol)
        })

        o.1c <- ode.1c %>% solve(params=c(V=20, CL=10), events=etSs)

        s.1c <- ode.1cs2 %>% solve(params=c(V=20, CL=10), events=etSs)

        s.2c <- ode.1cs %>% solve(theta=c(20, 10), events=etSs)

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
        }, linCmtSens=sens)

        o.2c <- ode.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=et)

        s.2c <- sol.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=et)

        test_that("2 compartment solved models and ODEs same.", {
            expect_equal(o.2c$C2, s.2c$C2,tolerance=tol)
        })

        o.2c <- ode.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=etSs)

        s.2c <- sol.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=etSs)

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
        }, linCmtSens=sens)

        o.3c <- ode.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et)

        s.3c <- sol.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et)

        test_that("3 compartment solved models and ODEs same.", {
            expect_equal(o.3c$C2, s.3c$C2, tolerance=tol)
        })

        o.3c <- ode.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=etSs)

        s.3c <- sol.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=etSs)

        test_that("3 compartment steady state solved models and ODEs same.", {
            expect_equal(o.3c$C2, s.3c$C2,tolerance=tol)
        })

        context(sprintf("Infusion + Bolus (%s)", .txt))

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
        }, linCmtSens=sens)

        ## Instead of specifying parameters in the solved system, you can
        ## specify them in the linCmt variable.
        ode.1cs2 <- RxODE({
            C2 = linCmt(CL, V);
        }, linCmtSens=sens)

        ## The solved systems can be mixed with ODE solving routines (to
        ## speed them up a bit...?)

        o.1c <- ode.1c %>% solve(params=c(V=20, CL=25), events=et)

        s.1c <- ode.1cs2 %>% solve(params=c(V=20, CL=25), events=et)

        s.2c <- ode.1cs %>% solve(theta=c(20, 25), events=et)

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
        }, linCmtSens=sens)

        o.2c <- ode.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=et)

        s.2c <- sol.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10), events=et)

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
        }, linCmtSens=sens)

        o.3c <- ode.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et)

        s.3c <- sol.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400), events=et)

        test_that("3 compartment solved models and ODEs same.", {
            expect_equal(o.3c$C2, s.3c$C2, tolerance=tol)
        })

        if(!advan){

            context(sprintf("Oral + Infusion + Bolus Models (%s)", .txt))

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
            }, linCmtSens=sens)

            o.1c <- ode.1c.ka %>% solve(params=c(V=20, CL=25, KA=2), events=et)
            s.1c <- sol.1c.ka %>% solve(params=c(V=20, CL=25, KA=2), events=et)

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
            }, linCmtSens=sens)

            o.2c <- ode.2c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, KA= 0.3), events=et)

            s.2c <- sol.2c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, KA=0.3), events=et)

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
            }, linCmtSens=sens)

            o.3c <- ode.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3), events=et)

            s.3c <- sol.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3), events=et)

            test_that("3 compartment solved models and ODEs same for mixed oral, iv and infusion.", {
                expect_equal(o.3c$C2, s.3c$C2,tolerance=tol)
            })

        }
        context(sprintf("Modeled bio-availability (%s)", .txt))

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
        }, linCmtSens=sens)

        for (fd in c(0.5,1,2)){
            for (fc in c(0.5,1,2)){
                o.1c <- ode.1c.ka %>% solve(params=c(V=20, CL=25, KA=2,fDepot=fd,fCenter=fc), events=et)
                s.1c <- sol.1c.ka %>% solve(params=c(V=20, CL=25, KA=2, fDepot=fd, fCenter=fc), events=et)
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
        }, linCmtSens=sens)

        sol.2c.ka <- RxODE({
            C2=linCmt(V, CL, V2, Q, KA);
            f(depot) = fDepot
            f(central) = fCenter
        }, linCmtSens=sens)

        for (fd in c(0.5,1,2)){
            for (fc in c(0.5,1,2)){
                o.2c <- ode.2c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, KA= 0.3, fDepot=fd, fCenter=fc), events=et)
                s.2c <- sol.2c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, KA=0.3, fDepot=fd, fCenter=fc), events=et)
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
        }, linCmtSens=sens)

        for (fd in c(0.5,1,2)){
            for (fc in c(0.5,1,2)){
                o.3c <- ode.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7,
                                                     V3=400, KA=0.3, fDepot=fd, fCenter=fc), events=et)
                s.3c <- sol.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3,
                                                     fDepot=fd, fCenter=fc), events=et)
                test_that(sprintf("3 compartment solved models and ODEs same for mixed oral, iv and infusion + Fd=%f,Fc=%f", fd,fc), {
                    expect_equal(o.3c$C2, s.3c$C2,tolerance=tol)
                })
            }
        }

        context(sprintf("Modeled lag time (%s)", .txt))

        et <- eventTable() %>%
            add.dosing(dose=3, rate=1.5, nbr.doses=3, dosing.interval=16,cmt=2) %>%
            add.dosing(dose=1.5, nbr.doses=3, dosing.interval=16,cmt=2) %>%
            add.dosing(dose=1.5, nbr.doses=3, dosing.interval=16,cmt=1,start.time=8) %>%
            add.sampling(seq(0, 48, length.out=200))

        ode.1c.ka <- RxODE({
            C2 = center/V;
            d / dt(depot) = -KA * depot
            d/dt(center) = KA * depot - CL*C2
            alag(depot) = lagDepot
            alag(center) = lagCenter
        }, linCmtSens=sens)

        sol.1c.ka <- RxODE({
            C2 = linCmt(V, CL, KA);
            alag(depot) = lagDepot
            alag(central) = lagCenter
        }, linCmtSens=sens)

        for (fd in c(1,2,10)){
            for (fc in c(1,2,10)){
                o.1c <- ode.1c.ka %>% solve(params=c(V=20, CL=25, KA=2,lagDepot=fd,lagCenter=fc), events=et)
                s.1c <- sol.1c.ka %>% solve(params=c(V=20, CL=25, KA=2, lagDepot=fd, lagCenter=fc), events=et)
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
        }, linCmtSens=sens)

        for (fd in c(1,2,10)){
            for (fc in c(1,2,10)){
                o.2c <- ode.2c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, KA= 0.3, lagDepot=fd, lagCenter=fc), events=et)
                s.2c <- sol.2c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, KA=0.3, lagDepot=fd, lagCenter=fc), events=et)
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
        }, linCmtSens=sens)

        sol.3c.ka <- RxODE({
            ## double solvedC(double t, int parameterization, int cmt, unsigned int col, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8);
            C2=linCmt(V, CL, V2, Q, Q2, V3, KA);
            alag(depot) = lagDepot
            alag(central) = lagCenter
        }, linCmtSens=sens)

        for (fd in c(1,2,10)){
            for (fc in c(1,2,10)){
                o.3c <- ode.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7,
                                                     V3=400, KA=0.3, lagDepot=fd, lagCenter=fc), events=et)
                s.3c <- sol.3c.ka %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, KA=0.3,
                                                     lagDepot=fd, lagCenter=fc), events=et)
                test_that(sprintf("3 compartment solved models and ODEs same for mixed oral, iv and infusion + Fd=%f,Fc=%f", fd,fc), {
                    expect_equal(o.3c$C2, s.3c$C2,tolerance=tol)
                })
            }
        }

        context(sprintf("Modeled rate (%s)", .txt))

        ode.1c <- RxODE({
            C2 = center/V;
            d/dt(center) = - CL*C2
            rate(center) = rt
        })

        sol.1c <- RxODE({
            C2 = linCmt(CL, V);
            rate(central) = rt
        }, linCmtSens=sens)

        et <- eventTable() %>%
            add.dosing(dose=3, rate=-1, nbr.doses=3, cmt=1,dosing.interval=12) %>%
            add.sampling(seq(0, 36, length.out=200))

        for (rt in seq(0.5, 1, 1.5)){
            o.1c <- ode.1c %>% solve(params=c(V=20, CL=25,rt=rt), events=et)
            s.1c <- sol.1c %>% solve(params=c(V=20, CL=25,rt=rt), events=et)
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
        }, linCmtSens=sens)

        sol.2c <- RxODE({
            C2=linCmt(V, CL, V2, Q);
            rate(central) = rt
        }, linCmtSens=sens)

        for (rt in seq(0.5, 1, 1.5)){
            o.2c <- ode.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, rt=rt), events=et)
            s.2c <- sol.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, rt=rt), events=et)
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
        }, linCmtSens=sens)

        for (rt in seq(0.5, 1, 1.5)){
            s.3c <- sol.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, rt=rt), events=et)
            o.3c <- ode.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, rt=rt), events=et)
            test_that(sprintf("3 compartment solved models and ODEs same for rate-modeled infusion: %s", rt), {
                expect_equal(o.3c$C2, s.3c$C2,tolerance=tol)
            })
        }

        context(sprintf("Modeled duration (%s)", .txt))

        ode.1c <- RxODE({
            C2 = center/V;
            d/dt(center) = - CL*C2
            dur(center) = dr
        })

        sol.1c <- RxODE({
            C2 = linCmt(CL, V);
            dur(central) = dr
        }, linCmtSens=sens)

        et <- eventTable() %>% add.dosing(dose=3, rate=-2, nbr.doses=3, cmt=1,dosing.interval=12) %>%
            add.sampling(seq(0, 36, length.out=200))

        for (dur in seq(0.5, 1, 1.5)){
            o.1c <- ode.1c %>% solve(params=c(V=20, CL=25,dr=dur), events=et)
            s.1c <- sol.1c %>% solve(params=c(V=20, CL=25,dr=dur), events=et)
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
        }, linCmtSens=sens)

        for (dur in seq(0.5, 1, 1.5)){
            o.2c <- ode.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, dr=dur), events=et)
            s.2c <- sol.2c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, dr=dur), events=et)
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
        }, linCmtSens=sens)

        for (dur in seq(0.5, 1, 1.5)){
            o.3c <- ode.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, dr=dur), events=et)
            s.3c <- sol.3c %>% solve(params=c(V=40, CL=18, V2=297, Q=10, Q2=7, V3=400, dr=dur), events=et)
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
        }, linCmtSens=sens)

        et <- eventTable() %>%
            add.dosing(dose=3, nbr.doses=6, dosing.interval=8) %>%
            add.sampling(0:48)

        s.1c <- ode.1cs2 %>% solve(params=c(V=20, CL=25,mt1=0.5, mt2=1.75),
                                   events=et)

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

    tol <<- 1e-5 ## Current difference for all equations

    context("1 cmt sensitivites")
    test_that("1 compartment sensitivities; IV bolus, Cl, V", {

        pred <- function () {return(Central)}

        pk <- function ()
        {
            lCl = THETA[1]
            lVc = THETA[2]
            prop.err = THETA[3]
            eta.Vc = ETA[1]
            eta.Cl = ETA[2]
            Vc <- exp(lVc + eta.Vc)
            Cl <- exp(lCl + eta.Cl)
        }

        err <- function ()
        {
            return(prop(prop.err))
        }

        mod <- RxODE({
            Central= linCmt(Vc, Cl);
        })

        pk1s <- rxSymPySetupPred(mod, predfn=pred, pkpars=pk, err=err)

        mod2 <- RxODE({
            Central = center/Vc;
            d/dt(center) = - Cl*Central
        })

        pk1o <- rxSymPySetupPred(mod2, predfn=pred, pkpars=pk, err=err)

        et <- eventTable() %>% add.dosing(dose=3, nbr.doses=6, dosing.interval=8) %>%
            add.sampling(seq(0, 48, length.out=200))


        parms <- c("THETA[1]" = 20, "THETA[2]" = 25, "ETA[1]"=1, "ETA[2]"=1,
                   "THETA[3]"=0.2)

        s1 <- rxSolve(pk1s$inner, parms, et)
        o1 <- rxSolve(pk1o$inner, parms, et)

        expect_equal(s1$rx_pred_, o1$rx_pred_)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_1___, o1$rx__sens_rx_pred__BY_ETA_1___)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_2___, o1$rx__sens_rx_pred__BY_ETA_2___)
        expect_equal(s1$rx_r_, o1$rx_r_)

        expect_equal(s1$rx__sens_rx_r__BY_ETA_1___, o1$rx__sens_rx_r__BY_ETA_1___)

        expect_equal(s1$rx__sens_rx_r__BY_ETA_2___, o1$rx__sens_rx_r__BY_ETA_2___)

        etSs  <- et() %>% et(amt=3) %>%
            et(time=4,amt=3, ss=1, ii=24) %>%
            et(amt=3, ss=2, ii=24, time=8) %>%
            et(seq(0,24,length.out=200))

        s1 <- rxSolve(pk1s$inner, parms, etSs)
        o1 <- rxSolve(pk1o$inner, parms, etSs)

        expect_equal(s1$rx_pred_, o1$rx_pred_)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_1___, o1$rx__sens_rx_pred__BY_ETA_1___)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_2___, o1$rx__sens_rx_pred__BY_ETA_2___)
        expect_equal(s1$rx_r_, o1$rx_r_)

        expect_equal(s1$rx__sens_rx_r__BY_ETA_1___, o1$rx__sens_rx_r__BY_ETA_1___)

        expect_equal(s1$rx__sens_rx_r__BY_ETA_2___, o1$rx__sens_rx_r__BY_ETA_2___)

        etInf <- eventTable() %>% add.dosing(dose=3, rate=1.5, nbr.doses=6, dosing.interval=8) %>%
            add.sampling(seq(0, 48, length.out=200))

        s1 <- rxSolve(pk1s$inner, parms, etInf)
        o1 <- rxSolve(pk1o$inner, parms, etInf)

        expect_equal(s1$rx_pred_, o1$rx_pred_)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_1___, o1$rx__sens_rx_pred__BY_ETA_1___)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_2___, o1$rx__sens_rx_pred__BY_ETA_2___)
        expect_equal(s1$rx_r_, o1$rx_r_)

        expect_equal(s1$rx__sens_rx_r__BY_ETA_1___, o1$rx__sens_rx_r__BY_ETA_1___)

        expect_equal(s1$rx__sens_rx_r__BY_ETA_2___, o1$rx__sens_rx_r__BY_ETA_2___)

        etInfSs <- et() %>% et(amt=3, rate=1.5) %>%
            et(time=4,amt=3, rate=1.5, ss=1, ii=24) %>%
            et(time=8, amt=3, rate=1.5, ss=2, ii=24) %>%
            et(seq(0,24,length.out=200))

        s1 <- rxSolve(pk1s$inner, parms, etInfSs)
        o1 <- rxSolve(pk1o$inner, parms, etInfSs)

        expect_equal(s1$rx_pred_, o1$rx_pred_)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_1___, o1$rx__sens_rx_pred__BY_ETA_1___)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_2___, o1$rx__sens_rx_pred__BY_ETA_2___)
        expect_equal(s1$rx_r_, o1$rx_r_)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_1___, o1$rx__sens_rx_r__BY_ETA_1___)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_2___, o1$rx__sens_rx_r__BY_ETA_2___)

    })

    ## Now oral

    test_that("1 compartment sensitivities; Oral Cl, V, Ka", {

        pred <- function () {return(Central)}

        pk <- function ()
        {
            lCl = THETA[1]
            lVc = THETA[2]
            lKa = THETA[3]
            prop.err = THETA[4]
            eta.Vc = ETA[1]
            eta.Cl = ETA[2]
            eta.Ka = ETA[3]
            Vc <- exp(lVc + eta.Vc)
            Cl <- exp(lCl + eta.Cl)
            Ka <- exp(lKa + eta.Ka)
        }

        err <- function ()
        {
            return(prop(prop.err))
        }

        mod <- RxODE({
            Central= linCmt(Vc, Cl, Ka);
        })

        pk1s <- rxSymPySetupPred(mod, predfn=pred, pkpars=pk, err=err)

        mod2 <- RxODE({
            Central = center/Vc;
            d/dt(depot) = -Ka * depot
            d/dt(center) = Ka * depot - Cl*Central
        })

        pk1o <- rxSymPySetupPred(mod2, predfn=pred, pkpars=pk, err=err)

        et <- eventTable() %>% add.dosing(dose=3, nbr.doses=6, dosing.interval=8) %>%
            add.sampling(seq(0, 48, length.out=200))

        parms <- c("THETA[1]" = log(20), "THETA[2]" = log(25), "THETA[3]"=log(2),
                   "ETA[1]"=0, "ETA[2]"=0, "ETA[3]"=0,
                   "THETA[4]"=0.2)

        s1 <- rxSolve(pk1s$inner, parms, et)

        o1 <- rxSolve(pk1o$inner, parms, et)

        expect_equal(s1$rx_pred_, o1$rx_pred_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_1___, o1$rx__sens_rx_pred__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_2___, o1$rx__sens_rx_pred__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_3___, o1$rx__sens_rx_pred__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx_r_, o1$rx_r_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_1___, o1$rx__sens_rx_r__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_2___, o1$rx__sens_rx_r__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_3___, o1$rx__sens_rx_r__BY_ETA_3___, tolerance=tol)

        etSs  <- et() %>% et(amt=3) %>%
            et(time=4,amt=3, ss=1, ii=24) %>%
            et(amt=3, ss=2, ii=24, time=8) %>%
            et(seq(0,24,length.out=200))

        s1 <- rxSolve(pk1s$inner, parms, etSs)
        o1 <- rxSolve(pk1o$inner, parms, etSs)

        expect_equal(s1$rx_pred_, o1$rx_pred_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_1___, o1$rx__sens_rx_pred__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_2___, o1$rx__sens_rx_pred__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_3___, o1$rx__sens_rx_pred__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx_r_, o1$rx_r_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_1___, o1$rx__sens_rx_r__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_2___, o1$rx__sens_rx_r__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_3___, o1$rx__sens_rx_r__BY_ETA_3___, tolerance=tol)

        etInf <- eventTable() %>%
            add.dosing(dose=3, rate=1.5, nbr.doses=6, dosing.interval=8, cmt=2) %>%
            add.sampling(seq(0, 48, length.out=200))

        s1 <- rxSolve(pk1s$inner, parms, etInf)
        o1 <- rxSolve(pk1o$inner, parms, etInf)

        expect_equal(s1$rx_pred_, o1$rx_pred_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_1___, o1$rx__sens_rx_pred__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_2___, o1$rx__sens_rx_pred__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_3___, o1$rx__sens_rx_pred__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx_r_, o1$rx_r_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_1___, o1$rx__sens_rx_r__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_2___, o1$rx__sens_rx_r__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_3___, o1$rx__sens_rx_r__BY_ETA_3___, tolerance=tol)

        etInfSs <- et() %>% et(amt=3, rate=1.5, cmt=2) %>%
            et(time=4,amt=3, rate=1.5, ss=1, ii=24, cmt=2) %>%
            et(time=8, amt=3, rate=1.5, ss=2, ii=24, cmt=2) %>%
            et(seq(0,24,length.out=200))

        s1 <- rxSolve(pk1s$inner, parms, etInfSs)

        o1 <- rxSolve(pk1o$inner, parms, etInfSs)

        expect_equal(s1$rx_pred_, o1$rx_pred_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_1___, o1$rx__sens_rx_pred__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_2___, o1$rx__sens_rx_pred__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_3___, o1$rx__sens_rx_pred__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx_r_, o1$rx_r_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_1___, o1$rx__sens_rx_r__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_2___, o1$rx__sens_rx_r__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_3___, o1$rx__sens_rx_r__BY_ETA_3___, tolerance=tol)

        etMix <- eventTable() %>%
            add.dosing(dose=3, rate=1.5, nbr.doses=3, dosing.interval=16,cmt=2) %>%
            add.dosing(dose=1.5, nbr.doses=3, dosing.interval=16,cmt=2) %>%
            add.dosing(dose=1.5, nbr.doses=3, dosing.interval=16,cmt=1,start.time=8) %>%
            add.sampling(seq(0, 48, length.out=200))

        s1 <- rxSolve(pk1s$inner, parms, etMix)
        o1 <- rxSolve(pk1o$inner, parms, etMix)

        expect_equal(s1$rx_pred_, o1$rx_pred_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_1___, o1$rx__sens_rx_pred__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_2___, o1$rx__sens_rx_pred__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_3___, o1$rx__sens_rx_pred__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx_r_, o1$rx_r_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_1___, o1$rx__sens_rx_r__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_2___, o1$rx__sens_rx_r__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_3___, o1$rx__sens_rx_r__BY_ETA_3___, tolerance=tol)
    })

    ##
    context("2 cmt sensitivites")
    test_that("2 compartment sensitivities; IV bolus Cl, Vc, Q, Vp", {

        pred <- function () {return(Central)}

        pk <- function () {
            lCl = THETA[1]
            lVc = THETA[2]
            lQ = THETA[3]
            lVp = THETA[4]
            prop.err = THETA[5]
            eta.Vc = ETA[1]
            eta.Cl = ETA[2]
            eta.Vp = ETA[3]
            eta.Q = ETA[4]
            Vc <- exp(lVc + eta.Vc)
            Cl <- exp(lCl + eta.Cl)
            Vp <- exp(lVp + eta.Vp)
            Q <- exp(lQ + eta.Q)
        }

        err <- function ()
        {
            return(prop(prop.err))
        }

        mod <- RxODE({
            Central=linCmt(Vc, Cl, Vp, Q);
        })

        pk2s <- rxSymPySetupPred(mod, predfn=pred, pkpars=pk, err=err)

        mod2 <- RxODE({
            Central = centr/Vc;
            C3 = peri/Vp;
            d/dt(centr) = - Cl*Central - Q*Central + Q*C3;
            d/dt(peri)  = Q*Central - Q*C3;
        })

        pk2o <- rxSymPySetupPred(mod2, predfn=pred, pkpars=pk, err=err)

        parms <- c("THETA[1]"=log(18), "THETA[2]"=log(40),
                   "THETA[3]"=log(10), "THETA[4]"=log(297),
                   "ETA[1]"=0, "ETA[2]"=0,
                   "ETA[3]"=0, "ETA[4]"=0,
                   "THETA[5]"=0.2)

        et <- eventTable() %>% add.dosing(dose=3, nbr.doses=6, dosing.interval=8) %>%
            add.sampling(seq(0, 48, length.out=200))

        s1 <- rxSolve(pk2s$inner, parms, et)
        o1 <- rxSolve(pk2o$inner, parms, et)

        expect_equal(s1$rx_pred_, o1$rx_pred_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_1___, o1$rx__sens_rx_pred__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_2___, o1$rx__sens_rx_pred__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_3___, o1$rx__sens_rx_pred__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_4___, o1$rx__sens_rx_pred__BY_ETA_4___, tolerance=tol)
        expect_equal(s1$rx_r_, o1$rx_r_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_1___, o1$rx__sens_rx_r__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_2___, o1$rx__sens_rx_r__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_3___, o1$rx__sens_rx_r__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_4___, o1$rx__sens_rx_r__BY_ETA_4___, tolerance=tol)

        etSs  <- et() %>% et(amt=3) %>%
            et(time=4,amt=3, ss=1, ii=24) %>%
            et(amt=3, ss=2, ii=24, time=8) %>%
            et(seq(0,24,length.out=200))

        s1 <- rxSolve(pk2s$inner, parms, etSs)
        o1 <- rxSolve(pk2o$inner, parms, etSs)

        expect_equal(s1$rx_pred_, o1$rx_pred_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_1___, o1$rx__sens_rx_pred__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_2___, o1$rx__sens_rx_pred__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_3___, o1$rx__sens_rx_pred__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_4___, o1$rx__sens_rx_pred__BY_ETA_4___, tolerance=tol)
        expect_equal(s1$rx_r_, o1$rx_r_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_1___, o1$rx__sens_rx_r__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_2___, o1$rx__sens_rx_r__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_3___, o1$rx__sens_rx_r__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_4___, o1$rx__sens_rx_r__BY_ETA_4___, tolerance=tol)

        etInf <- eventTable() %>%
            add.dosing(dose=3, rate=1.5, nbr.doses=6, dosing.interval=8, cmt=2) %>%
            add.sampling(seq(0, 48, length.out=200))

        s1 <- rxSolve(pk2s$inner, parms, etInf)
        o1 <- rxSolve(pk2o$inner, parms, etInf)

        ## Close to o1
        expect_true(all(s1$rx_pred_ == 0))
        expect_true(all(s1$rx__sens_rx_pred__BY_ETA_1___ == 0))
        expect_true(all(s1$rx__sens_rx_pred__BY_ETA_2___ == 0))
        expect_true(all(s1$rx__sens_rx_pred__BY_ETA_3___ == 0))
        expect_true(all(s1$rx__sens_rx_pred__BY_ETA_4___ == 0))

        expect_equal(s1$rx_r_, o1$rx_r_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_1___, o1$rx__sens_rx_r__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_2___, o1$rx__sens_rx_r__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_3___, o1$rx__sens_rx_r__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_4___, o1$rx__sens_rx_r__BY_ETA_4___, tolerance=tol)

        etInfSs <- et() %>% et(amt=3, rate=1.5) %>%
            et(time=4,amt=3, rate=1.5, ss=1, ii=24) %>%
            et(time=8, amt=3, rate=1.5, ss=2, ii=24) %>%
            et(seq(0,24,length.out=200))

        s1 <- rxSolve(pk2s$inner, parms, etInfSs)
        o1 <- rxSolve(pk2o$inner, parms, etInfSs)

        expect_equal(s1$rx_pred_, o1$rx_pred_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_1___, o1$rx__sens_rx_pred__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_2___, o1$rx__sens_rx_pred__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_3___, o1$rx__sens_rx_pred__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_4___, o1$rx__sens_rx_pred__BY_ETA_4___, tolerance=tol)
        expect_equal(s1$rx_r_, o1$rx_r_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_1___, o1$rx__sens_rx_r__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_2___, o1$rx__sens_rx_r__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_3___, o1$rx__sens_rx_r__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_4___, o1$rx__sens_rx_r__BY_ETA_4___, tolerance=tol)

    })


    test_that("2 compartment sensitivities; Oral Cl, Vc, Q, Vp, Ka", {

        pred <- function () {return(Central)}

        pk <- function () {
            lCl = THETA[1]
            lVc = THETA[2]
            lQ = THETA[3]
            lVp = THETA[4]
            lKa = THETA[5]
            prop.err = THETA[6]
            eta.Vc = ETA[1]
            eta.Cl = ETA[2]
            eta.Vp = ETA[3]
            eta.Q = ETA[4]
            eta.Ka = ETA[5]
            Vc <- exp(lVc + eta.Vc)
            Cl <- exp(lCl + eta.Cl)
            Vp <- exp(lVp + eta.Vp)
            Q <- exp(lQ + eta.Q)
            Ka <- exp(lKa + eta.Ka)
        }

        err <- function ()
        {
            return(prop(prop.err))
        }

        mod <- RxODE({
            Central=linCmt(Vc, Cl, Vp, Q, Ka);
        })

        pk2s <- rxSymPySetupPred(mod, predfn=pred, pkpars=pk, err=err)

        mod2 <- RxODE({
            Central = centr/Vc;
            C3 = peri/Vp;
            d/dt(depot) =-Ka*depot;
            d/dt(centr) = Ka*depot - Cl*Central - Q*Central + Q*C3;
            d/dt(peri)  = Q*Central - Q*C3;
        })

        pk2o <- rxSymPySetupPred(mod2, predfn=pred, pkpars=pk, err=err)

        parms <- c("THETA[1]"=log(18), ## Cl
                   "THETA[2]"=log(40), ## Vc
                   "THETA[3]"=log(10), ## Q
                   "THETA[4]"=log(297),## Vp
                   "THETA[5]"=log(0.3), ## Ka
                   "ETA[1]"=0, "ETA[2]"=0,
                   "ETA[3]"=0, "ETA[4]"=0,
                   "ETA[5]"=0,
                   "THETA[6]"=0.2)

        et <- eventTable() %>% add.dosing(dose=3, nbr.doses=6, dosing.interval=8) %>%
            add.sampling(seq(0, 48, length.out=200))

        s1 <- rxSolve(pk2s$inner, parms, et)
        o1 <- rxSolve(pk2o$inner, parms, et)

        expect_equal(s1$rx_pred_, o1$rx_pred_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_1___, o1$rx__sens_rx_pred__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_2___, o1$rx__sens_rx_pred__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_3___, o1$rx__sens_rx_pred__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_4___, o1$rx__sens_rx_pred__BY_ETA_4___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_5___, o1$rx__sens_rx_pred__BY_ETA_5___, tolerance=tol)
        expect_equal(s1$rx_r_, o1$rx_r_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_1___, o1$rx__sens_rx_r__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_2___, o1$rx__sens_rx_r__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_3___, o1$rx__sens_rx_r__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_4___, o1$rx__sens_rx_r__BY_ETA_4___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_5___, o1$rx__sens_rx_r__BY_ETA_5___, tolerance=tol)

        etSs  <- et() %>% et(amt=3) %>%
            et(time=4,amt=3, ss=1, ii=24) %>%
            et(amt=3, ss=2, ii=24, time=8) %>%
            et(seq(0,24,length.out=200))

        s1 <- rxSolve(pk2s$inner, parms, etSs)
        o1 <- rxSolve(pk2o$inner, parms, etSs)

        expect_equal(s1$rx_pred_, o1$rx_pred_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_1___, o1$rx__sens_rx_pred__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_2___, o1$rx__sens_rx_pred__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_3___, o1$rx__sens_rx_pred__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_4___, o1$rx__sens_rx_pred__BY_ETA_4___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_5___, o1$rx__sens_rx_pred__BY_ETA_5___, tolerance=tol)
        expect_equal(s1$rx_r_, o1$rx_r_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_1___, o1$rx__sens_rx_r__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_2___, o1$rx__sens_rx_r__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_3___, o1$rx__sens_rx_r__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_4___, o1$rx__sens_rx_r__BY_ETA_4___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_5___, o1$rx__sens_rx_r__BY_ETA_5___, tolerance=tol)

        etInf <- eventTable() %>%
            add.dosing(dose=3, rate=1.5, nbr.doses=6, dosing.interval=8, cmt=2) %>%
            add.sampling(seq(0, 48, length.out=200))

        s1 <- rxSolve(pk2s$inner, parms, etInf)
        o1 <- rxSolve(pk2o$inner, parms, etInf)

        expect_equal(s1$rx_pred_, o1$rx_pred_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_1___, o1$rx__sens_rx_pred__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_2___, o1$rx__sens_rx_pred__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_3___, o1$rx__sens_rx_pred__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_4___, o1$rx__sens_rx_pred__BY_ETA_4___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_5___, o1$rx__sens_rx_pred__BY_ETA_5___, tolerance=tol)
        expect_equal(s1$rx_r_, o1$rx_r_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_1___, o1$rx__sens_rx_r__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_2___, o1$rx__sens_rx_r__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_3___, o1$rx__sens_rx_r__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_4___, o1$rx__sens_rx_r__BY_ETA_4___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_5___, o1$rx__sens_rx_r__BY_ETA_5___, tolerance=tol)

        etInfSs <- et() %>% et(amt=3, rate=1.5, cmt=2) %>%
            et(time=4,amt=3, rate=1.5, ss=1, ii=24, cmt=2) %>%
            et(time=8, amt=3, rate=1.5, ss=2, ii=24, cmt=2) %>%
            et(seq(0,24,length.out=200))

        s1 <- rxSolve(pk2s$inner, parms, etInfSs)
        o1 <- rxSolve(pk2o$inner, parms, etInfSs)

        expect_equal(s1$rx_pred_, o1$rx_pred_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_1___, o1$rx__sens_rx_pred__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_2___, o1$rx__sens_rx_pred__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_3___, o1$rx__sens_rx_pred__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_4___, o1$rx__sens_rx_pred__BY_ETA_4___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_5___, o1$rx__sens_rx_pred__BY_ETA_5___, tolerance=tol)
        expect_equal(s1$rx_r_, o1$rx_r_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_1___, o1$rx__sens_rx_r__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_2___, o1$rx__sens_rx_r__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_3___, o1$rx__sens_rx_r__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_4___, o1$rx__sens_rx_r__BY_ETA_4___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_5___, o1$rx__sens_rx_r__BY_ETA_5___, tolerance=tol)

        etMix <- eventTable() %>%
            add.dosing(dose=3, rate=1.5, nbr.doses=3, dosing.interval=16,cmt=2) %>%
            add.dosing(dose=1.5, nbr.doses=3, dosing.interval=16,cmt=2) %>%
            add.dosing(dose=1.5, nbr.doses=3, dosing.interval=16,cmt=1,start.time=8) %>%
            add.sampling(seq(0, 48, length.out=200))

        s1 <- rxSolve(pk2s$inner, parms, etMix)
        o1 <- rxSolve(pk2o$inner, parms, etMix)

        expect_equal(s1$rx_pred_, o1$rx_pred_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_1___, o1$rx__sens_rx_pred__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_2___, o1$rx__sens_rx_pred__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_3___, o1$rx__sens_rx_pred__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_4___, o1$rx__sens_rx_pred__BY_ETA_4___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_pred__BY_ETA_5___, o1$rx__sens_rx_pred__BY_ETA_5___, tolerance=tol)
        expect_equal(s1$rx_r_, o1$rx_r_, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_1___, o1$rx__sens_rx_r__BY_ETA_1___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_2___, o1$rx__sens_rx_r__BY_ETA_2___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_3___, o1$rx__sens_rx_r__BY_ETA_3___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_4___, o1$rx__sens_rx_r__BY_ETA_4___, tolerance=tol)
        expect_equal(s1$rx__sens_rx_r__BY_ETA_5___, o1$rx__sens_rx_r__BY_ETA_5___, tolerance=tol)

    })

}, silent=TRUE, test="lincmt")


}
