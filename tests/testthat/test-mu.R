rxodeTest({

  test_that("bounded functions needs numeric bounds", {

    lmat <- lotri({
      tka <- 0.2
      tcl <- 0.2
      tv <- 0.1
      eta.ka ~ 0.1
      eta.cl ~ 0.1
      eta.v ~ 0.1
      add.sd <- 0.1
    })

    testBounded <- function(type="expit") {

      expect_error(rxMuRef(paste0("a=", type, "(tka + eta.ka, a, b)"), lmat))

      expect_error(rxMuRef(paste0("a=", type, "(tka + eta.ka, 1, b)"), lmat))


      expect_error(rxMuRef(paste0("a=", type, "(tka + eta.ka, 1, b)"), lmat))

      expect_error(rxMuRef(paste0("a=", type, "(tka + eta.ka, 1, 2)"), lmat), NA)

      expect_error(rxMuRef(paste0("a=", type, "(tka + eta.ka, 2, 1)"), lmat))

      expect_error(rxMuRef(paste0("a=", type, "(tka + eta.ka, 0.5)"), lmat), NA)

      expect_error(rxMuRef(paste0("a=", type, "(tka + eta.ka, a)"), lmat))

      expect_error(rxMuRef(paste0("a=", type, "(tka + eta.ka, 4)"), lmat))
    }

    testBounded("logit")
    testBounded("expit")
    testBounded("probit")
    testBounded("probitInv")

  })

  lmat <- lotri({
    theta1 <- 1
    theta2 <- 1
    theta3 <- 1
    eta1 ~ 0.1
  })

  test_that("bad mu referencing examples (throw error)", {

    expect_error(rxMuRef("a=theta1+theta2+theta3*wt+eta1", lmat))

    expect_error(rxMuRef("a=theta1+theta2*wt+theta3*wt+eta1", lmat))
  })

  test_that("simple mu referencing", {

    lmat <- lotri({
      ## You may label each parameter with a comment
      tka <- 0.45 # Log Ka
      tcl <- log(c(0, 2.7, 100)) # Log Cl
      ## This works with interactive models
      ## You may also label the preceding line with label("label text")
      tv <- 3.45; label("log V")
      ## the label("Label name") works with all models
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })

    env <- rxMuRef(RxODE({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- tv + eta.v
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl/v * center
      cp = center/v
      ## cp ~ add(add.sd)
    }), lmat)

    expect_equal(env$muRefDataFrame, structure(list(theta = c("tka", "tcl", "tv"),
                                                    eta = c("eta.ka", "eta.cl", "eta.v"),
                                                    curEval = c("exp", "exp", "")),
                                               row.names = c(NA, -3L),
                                               class = "data.frame"))

  })

    test_that("shared eta is not mu referencing", {

      lmat <- lotri({
        t.EmaxA = 1
        t.EmaxB = 2
        t.EmaxC = 3
        eta.emax ~ 0.1
        add.sd <- 0.7
      })

      ## ## Test a duplicated eta; It shouldn't be counted as mu-referenced
      env <- rxMuRef(RxODE({
        EmaxA <- exp(t.EmaxA + eta.emax)
        EmaxB <- exp(t.EmaxB + eta.emax)
        EmaxC <- exp(t.EmaxC + eta.emax)
      }), lmat)

      expect_equal(length(env$muRefDataFrame[,1]), 0L)
      expect_equal(env$nonMuEtas,
                   structure(list(eta = "eta.emax", curEval = "exp"),
                             row.names = c(NA, -1L), class = "data.frame"))


      env <- rxMuRef(RxODE({
        EmaxA <- exp(t.EmaxA + eta.emax)
        EmaxB <- exp(t.EmaxB + eta.emax)
        EmaxC <- t.EmaxC + eta.emax
      }), lmat)

      expect_equal(length(env$muRefDataFrame[,1]), 0L)
      expect_equal(env$nonMuEtas,
                   structure(list(eta = "eta.emax", curEval = ""),
                             row.names = c(NA, -1L), class = "data.frame"))

      env <- rxMuRef(RxODE({
        EmaxB <- t.EmaxB + eta.emax
        EmaxA <- exp(t.EmaxA + eta.emax)
      }), lmat)

      expect_equal(length(env$muRefDataFrame[,1]), 0L)
      expect_equal(env$nonMuEtas,
                   structure(list(eta = "eta.emax", curEval = ""),
                             row.names = c(NA, -1L), class = "data.frame"))


    })


    test_that("composite ode expressions", {

      lmat <- lotri({
        ## You may label each parameter with a comment
        tka <- 0.45 # Log Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })


      env <- rxMuRef(RxODE({
        d/dt(depot) = -exp(tka + eta.ka) * depot
        d/dt(center) = exp(tka + eta.ka) * depot - exp(tcl + eta.cl)/exp(tv + eta.v) * center
        cp = center/exp(tv + eta.v)
        #cp ~ add(add.sd)
      }), lmat)

      expect_equal(length(env$nonMuEtas[, 1]), 0L)
      expect_equal(env$muRefDataFrame,
                   structure(list(theta = c("tka", "tcl", "tv"),
                                  eta = c("eta.ka", "eta.cl", "eta.v"),
                                  curEval = c("exp", "exp", "exp")),
                             row.names = c(NA, -3L), class = "data.frame"))

    })


    test_that("old style tka*eta(eta.ka)", {

      lmat <- lotri({
        ## You may label each parameter with a comment
        tka <- 0.45 # Log Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })

      env <- rxMuRef(RxODE({
        ka <- tka * exp(eta.ka + 0)
        cl <- tcl * exp(eta.cl + 0)
        v <- tv * exp(eta.v + 0)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl/v * center
        cp = center/v
        ## cp ~ add(add.sd)
      }), lmat)

      expect_equal(env$nonMuEtas,
                   structure(list(eta = c("eta.ka", "eta.cl", "eta.v"),
                                  curEval = c("exp", "exp", "exp")),
                             row.names = c(NA_integer_, -3L), class = "data.frame"))
      expect_equal(length(env$muRefDataFrame[, 1]), 0L)

      env <- rxMuRef(RxODE({
        ka <- tka * exp(eta.ka + 0)
        cl <- tcl * exp(eta.cl + 0)
        v <- tv * exp(eta.v + 0)
        v2 <- tv + eta.v
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl/v * center
        cp = center/v
        ## cp ~ add(add.sd)
      }), lmat)

      expect_equal(env$nonMuEtas,
                   structure(list(eta = c("eta.ka", "eta.cl", "eta.v"),
                                  curEval = c("exp", "exp", "")),
                             row.names = c(NA_integer_, -3L), class = "data.frame"))
      expect_equal(length(env$muRefDataFrame[, 1]), 0L)

      env <- rxMuRef(RxODE({
        ka <- tka * exp(eta.ka)
        cl <- tcl * exp(eta.cl)
        v <- tv * exp(eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl/v * center
        cp = center/v
        ## cp ~ add(add.sd)
      }), lmat)

      expect_equal(env$nonMuEtas,
                   structure(list(eta = c("eta.ka", "eta.cl", "eta.v"),
                                  curEval = c("exp", "exp", "exp")),
                             row.names = c(NA_integer_, -3L), class = "data.frame"))
      expect_equal(length(env$muRefDataFrame[, 1]), 0L)

      env <- rxMuRef(RxODE({
        ka <- tka * exp(eta.ka)
        cl <- tcl * exp(eta.cl)
        v <- tv * exp(eta.v)
        etv <- eta.v
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl/v * center
        cp = center/v
        ## cp ~ add(add.sd)
      }), lmat)

      expect_equal(env$nonMuEtas,
                   structure(list(eta = c("eta.ka", "eta.cl", "eta.v"),
                                  curEval = c("exp", "exp", "")),
                             row.names = c(NA_integer_, -3L), class = "data.frame"))
      expect_equal(length(env$muRefDataFrame[, 1]), 0L)



    })





}, test="lvl2")


## rxMuRef(RxODE({
##   ka <- exp(tka + eta.ka)
##   cl <- exp(tcl + eta.cl)
##   v <- exp(tv + eta.v)
##   d/dt(depot) = -ka * depot
##   d/dt(center) = ka * depot - cl/v * center
##   cp = center/v
##   ## cp ~ add(add.sd)
## }), theta=c("tka", "tcl", "tv", "add.sd"),
## eta=c("eta.ka", "eta.cl", "eta.v"))



## ## Composite expressions should be extracted to their own lines
## rxMuRef(RxODE({
##   ratio <- exp(t.EmaxA + eta.emaxA) / exp(t.EmaxB + eta.emaxB)
## }), theta=c("tka", "tcl", "tv", "add.sd"),
## eta=c("eta.ka", "eta.cl", "eta.v"))

## ## This composite RxODE model should be extracted to something like above
## rxMuRef(RxODE({
##   d/dt(depot) = -exp(tka + eta.ka) * depot
##   d/dt(center) = exp(tka + eta.ka) * depot - exp(tcl + eta.cl)/exp(tv + eta.v) * center
##   cp = center/exp(tv + eta.v)
##   #cp ~ add(add.sd)
## }), theta=c("tka", "tcl", "tv", "add.sd"),
## eta=c("eta.ka", "eta.cl", "eta.v"))


## ## This should be expanded to eta.ka mu-referenced variables
## rxMuRef(RxODE({
##   ka <- tka * exp(eta.ka)
##   cl <- tcl * exp(eta.cl)
##   v <- tv * exp(eta.v)
##   d/dt(depot) = -ka * depot
##   d/dt(center) = ka * depot - cl/v * center
##   cp = center/v
##   ## cp ~ add(add.sd)
## }), theta=c("tka", "tcl", "tv", "add.sd"),
## eta=c("eta.ka", "eta.cl", "eta.v"))
