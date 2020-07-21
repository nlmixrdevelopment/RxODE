rxPermissive({

  context("tad family of functions with odes")

  test_that("tad family works with ode", {

    mod1 <- RxODE({
      KA <- 2.94E-01
      CL <- 1.86E+01
      V2 <- 4.02E+01
      Q <- 1.05E+01
      V3 <- 2.97E+02
      Kin <- 1
      Kout <- 1
      EC50 <- 200
      C2 <- centr / V2
      C3 <- peri / V3
      d/dt(depot) <- -KA * depot
      d/dt(centr) <- KA * depot - CL * C2 - Q * C2 + Q * C3
      d/dt(peri) <- Q * C2 - Q * C3
      d/dt(eff) <- Kin - Kout * (1 - C2 / (EC50 + C2)) * eff
      ## TAD tests
      tad <- tad()
      tl <- tlast()
      tafd <- tafd()
      tf <- tfirst()
      tadd <- tad(depot)
      tfirstd <- tfirst(depot)
      tlastd <- tlast(depot)
      tafdd <- tafd(depot)
      tadp <- tad(peri)
      tafdp <- tafd(peri)
      tfirstp <- tfirst(peri)
      tlastp <- tlast(peri)
    })

    ev <- et(amountUnits = "mg", timeUnits = "hours") %>%
      et(time=1, amt = 10000, addl = 9, ii = 12, cmt = "depot") %>%
      et(time = 120, amt = 2000, addl = 4, ii = 14, cmt = "depot") %>%
      et(time = 122, amt = 2000, addl = 4, ii = 14, cmt = "peri") %>%
      et(0, 240, by=3)

    r1 <- rxSolve(mod1, ev, addDosing = TRUE)

    expect_equal(r1$tad,
                 c(NA, 0, 2, 5, 8, 11, 0, 2, 5, 8, 11, 0, 2, 5, 8, 11, 0, 2, 5,
                   8, 11, 0, 2, 5, 8, 11, 0, 2, 5, 8, 11, 0, 2, 5, 8, 11, 0, 2,
                   5, 8, 11, 0, 2, 5, 8, 11, 0, 2, 5, 8, 0, 0, 0, 1, 4, 7, 10, 0,
                   1, 0, 2, 5, 8, 11, 0, 0, 0, 3, 6, 9, 0, 0, 0, 1, 4, 7, 10, 0,
                   1, 0, 2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44,
                   47, 50, 53, 56, 59, 62))

    expect_equal(r1$tafd, c(NA, 0, 2, 5, 8, 11, 12, 14, 17, 20, 23, 24, 26, 29, 32, 35,
                            36, 38, 41, 44, 47, 48, 50, 53, 56, 59, 60, 62, 65, 68, 71, 72,
                            74, 77, 80, 83, 84, 86, 89, 92, 95, 96, 98, 101, 104, 107, 108,
                            110, 113, 116, 119, 119, 121, 122, 125, 128, 131, 133, 134, 135,
                            137, 140, 143, 146, 147, 149, 149, 152, 155, 158, 161, 161, 163,
                            164, 167, 170, 173, 175, 176, 177, 179, 182, 185, 188, 191, 194,
                            197, 200, 203, 206, 209, 212, 215, 218, 221, 224, 227, 230, 233,
                            236, 239))

    expect_equal(r1$tadd,
                 c(NA, 0, 2, 5, 8, 11, 0, 2, 5, 8, 11, 0, 2, 5, 8, 11, 0, 2, 5,
                   8, 11, 0, 2, 5, 8, 11, 0, 2, 5, 8, 11, 0, 2, 5, 8, 11, 0, 2,
                   5, 8, 11, 0, 2, 5, 8, 11, 0, 2, 5, 8, 0, 0, 2, 3, 6, 9, 12, 0,
                   1, 2, 4, 7, 10, 13, 0, 2, 2, 5, 8, 11, 0, 0, 2, 3, 6, 9, 12,
                   0, 1, 2, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43,
                   46, 49, 52, 55, 58, 61, 64))

    expect_equal(r1$tafdd, c(NA, 0, 2, 5, 8, 11, 12, 14, 17, 20, 23, 24, 26, 29, 32, 35,
                             36, 38, 41, 44, 47, 48, 50, 53, 56, 59, 60, 62, 65, 68, 71, 72,
                             74, 77, 80, 83, 84, 86, 89, 92, 95, 96, 98, 101, 104, 107, 108,
                             110, 113, 116, 119, 119, 121, 122, 125, 128, 131, 133, 134, 135,
                             137, 140, 143, 146, 147, 149, 149, 152, 155, 158, 161, 161, 163,
                             164, 167, 170, 173, 175, 176, 177, 179, 182, 185, 188, 191, 194,
                             197, 200, 203, 206, 209, 212, 215, 218, 221, 224, 227, 230, 233,
                             236, 239))

    expect_equal(r1$tadp,
                 c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                   NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                   NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                   NA, NA, NA, NA, NA, 0, 1, 4, 7, 10, 12, 13, 0, 2, 5, 8, 11, 12,
                   0, 0, 3, 6, 9, 12, 12, 0, 1, 4, 7, 10, 12, 13, 0, 2, 5, 8, 11,
                   14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53, 56, 59,
                   62))

    expect_equal(r1$tfirst[1], NA_real_)
    expect_true(all(r1$tfirst[-1] == 1))

    expect_equal(r1$tlast,
                 c(NA, 1, 1, 1, 1, 1, 13, 13, 13, 13, 13, 25, 25, 25, 25, 25,
                   37, 37, 37, 37, 37, 49, 49, 49, 49, 49, 61, 61, 61, 61, 61, 73,
                   73, 73, 73, 73, 85, 85, 85, 85, 85, 97, 97, 97, 97, 97, 109,
                   109, 109, 109, 120, 120, 120, 120, 120, 120, 120, 134, 134, 134,
                   134, 134, 134, 134, 148, 148, 148, 148, 148, 148, 162, 162, 162,
                   162, 162, 162, 162, 176, 176, 176, 176, 176, 176, 176, 176, 176,
                   176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176,
                   176, 176))

    expect_equal(r1$tfirstd,
                 c(NA, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                   1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))


    expect_equal(r1$tlastd,
                 c(NA, 1, 1, 1, 1, 1, 13, 13, 13, 13, 13, 25, 25, 25, 25, 25,
                   37, 37, 37, 37, 37, 49, 49, 49, 49, 49, 61, 61, 61, 61, 61, 73,
                   73, 73, 73, 73, 85, 85, 85, 85, 85, 97, 97, 97, 97, 97, 109,
                   109, 109, 109, 120, 120, 120, 120, 120, 120, 120, 134, 134, 134,
                   134, 134, 134, 134, 148, 148, 148, 148, 148, 148, 162, 162, 162,
                   162, 162, 162, 162, 176, 176, 176, 176, 176, 176, 176, 176, 176,
                   176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176,
                   176, 176))


    expect_equal(r1$tlastp,
                 c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                   NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                   NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                   NA, NA, NA, NA, NA, 122, 122, 122, 122, 122, 122, 122, 136, 136,
                   136, 136, 136, 136, 150, 150, 150, 150, 150, 150, 150, 164, 164,
                   164, 164, 164, 164, 164, 178, 178, 178, 178, 178, 178, 178, 178,
                   178, 178, 178, 178, 178, 178, 178, 178, 178, 178, 178, 178, 178,
                   178))

    expect_equal(r1$tfirstp,
                 c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                   NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                   NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                   NA, NA, NA, NA, NA, 122, 122, 122, 122, 122, 122, 122, 122, 122,
                   122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122,
                   122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122,
                   122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122,
                   122))

  })

  test_that("test parsing of ode", {

    expect_error(RxODE({
      KA <- 2.94E-01
      CL <- 1.86E+01
      V2 <- 4.02E+01
      Q <- 1.05E+01
      V3 <- 2.97E+02
      Kin <- 1
      Kout <- 1
      EC50 <- 200
      C2 <- centr / V2
      C3 <- peri / V3
      ## error in capturing early
      f <- tad(eff)
      d/dt(depot) <- -KA * depot
      d/dt(centr) <- KA * depot - CL * C2 - Q * C2 + Q * C3
      d/dt(peri) <- Q * C2 - Q * C3
      d/dt(eff) <- Kin - Kout * (1 - C2 / (EC50 + C2)) * eff
    }))

    ## now change the option
    op <- options()
    options(RxODE.syntax.require.ode.first=FALSE)
    expect_error(RxODE({
      KA <- 2.94E-01
      CL <- 1.86E+01
      V2 <- 4.02E+01
      Q <- 1.05E+01
      V3 <- 2.97E+02
      Kin <- 1
      Kout <- 1
      EC50 <- 200
      C2 <- centr / V2
      C3 <- peri / V3
      ## error in capturing early
      f <- tad(eff)
      d/dt(depot) <- -KA * depot
      d/dt(centr) <- KA * depot - CL * C2 - Q * C2 + Q * C3
      d/dt(peri) <- Q * C2 - Q * C3
      d/dt(eff) <- Kin - Kout * (1 - C2 / (EC50 + C2)) * eff
    }), NA)

    options(op)

  })

  context("tad family of functions with linCmt()")


  test_that("lincmt solution tad family", {

    sol.1c.ka <- RxODE({
      KA = 2
      V = 20
      CL = 25
      C2 <- linCmt(V, CL, KA)
      tad <- tad()
      tafd <- tafd()
      tadd <- tad(depot)
      tafdd <- tafd(depot)
      tadc <- tad(central)
      tafdc <- tafd(central)
    })

    et <- eventTable() %>%
      add.dosing(dose = 3, nbr.doses = 6, dosing.interval = 8) %>%
      add.dosing(dose=6, nbr.doses=6, dosing.interval = 8,
                 start.time = 2, dosing.to = "central") %>%
      add.sampling(seq(0, 48, length.out = 200))

    s1 <- rxSolve(sol.1c.ka, et)

    sol.1c.ka <- RxODE({
      KA = 2
      V = 20
      CL = 25
      C2 <- linCmt(V, CL, KA)
      tad <- tad()
      tafd <- tafd()
      tadc <- tad(central)
      tafdc <- tafd(central)
    })

    s2 <- rxSolve(sol.1c.ka, et)

    expect_equal(s2$tad, s1$tad)
    expect_equal(s2$tafd, s1$tafd)
    expect_equal(s2$tadc, s1$tadc)
    expect_equal(s2$tafdc, s1$tafdc)


  })

  context("tad family of functions with linCmt()/ode mix")





})
