rxodeTest({

  test_that("bounded functions needs numeric bounds", {

    testBounded <- function(type="expit") {

      expect_error(rxMuRef(paste0("a=", type, "(tka + eta.ka, a, b)"), theta=c("tka", "tcl", "tv", "add.sd"),
                           eta=c("eta.ka", "eta.cl", "eta.v")))

      expect_error(rxMuRef(paste0("a=", type, "(tka + eta.ka, 1, b)"), theta=c("tka", "tcl", "tv", "add.sd"),
                           eta=c("eta.ka", "eta.cl", "eta.v")))

      expect_error(rxMuRef(paste0("a=", type, "(tka + eta.ka, 1, b)"), theta=c("tka", "tcl", "tv", "add.sd"),
                           eta=c("eta.ka", "eta.cl", "eta.v")))

      expect_error(rxMuRef(paste0("a=", type, "(tka + eta.ka, 1, 2)"), theta=c("tka", "tcl", "tv", "add.sd"),
                           eta=c("eta.ka", "eta.cl", "eta.v")), NA)

      expect_error(rxMuRef(paste0("a=", type, "(tka + eta.ka, 2, 1)"), theta=c("tka", "tcl", "tv", "add.sd"),
                           eta=c("eta.ka", "eta.cl", "eta.v")))

      expect_error(rxMuRef(paste0("a=", type, "(tka + eta.ka, 0.5)"), theta=c("tka", "tcl", "tv", "add.sd"),
                           eta=c("eta.ka", "eta.cl", "eta.v")), NA)

      expect_error(rxMuRef(paste0("a=", type, "(tka + eta.ka, a)"), theta=c("tka", "tcl", "tv", "add.sd"),
                           eta=c("eta.ka", "eta.cl", "eta.v")))

      expect_error(rxMuRef(paste0("a=", type, "(tka + eta.ka, 4)"), theta=c("tka", "tcl", "tv", "add.sd"),
                           eta=c("eta.ka", "eta.cl", "eta.v")))
    }

    testBounded("logit")
    testBounded("expit")
    testBounded("probit")
    testBounded("probitInv")

  })

  test_that("bad mu referencing examples (throw error)", {

    expect_error(rxMuRef("a=theta1+theta2+theta3*wt+eta1", theta=c("theta1", "theta2", "theta3"),
                         eta="eta1"))

    expect_error(rxMuRef("a=theta1+theta2*wt+theta3*wt+eta1", theta=c("theta1", "theta2", "theta3"),
                         eta="eta1"))
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


## ## Test a duplicated eta; It shouldn't be counted as mu-referenced
## rxMuRef(RxODE({
##   EmaxA <- exp(t.EmaxA + eta.emax)
##   EmaxB <- exp(t.EmaxB + eta.emax)
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
