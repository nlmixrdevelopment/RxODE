rxodeTest({

  lmat <- lotri({
    ## You may label each parameter with a comment
    tka <- 0.45 # Log Ka
    tcl <- log(c(0, 2.7, 100)) # Log Cl
    ## This works with interactive models
    ## You may also label the preceding line with label("label text")
    tv <- 3.45; label("log V")
    tvp <- 3.45; label("log V")
    cl.wt <- 0.1
    v.wt <- 0.1
    cl.sex <- 0.1
    v.sex <- 0.1
    cl.age <- 0.1
    v.age <- 0.1
    vp.wt <- 1
    vp.sex <- 1
    vp.age <- 1
    ## the label("Label name") works with all models
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
    pow.sd <- 0.2
    pow  <- 1
    lambda  <- 0.5
  })

  df <- as.data.frame(lmat)

  test_that("error when errors have too many arguments", {
    expect_error(.errProcessExpression(quote({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + log(wt / 70) * cl.wt + sex * cl.sex + age * cl.age + 3)
      v  <- exp(tv + eta.v + wt * v.wt + sex * v.sex + age * v.age + 2)
      vp <- exp(tvp + wt * vp.wt + sex * vp.sex + age * vp.age)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl/v * center
      cp = center/v
      cp ~ add(add.sd, pow) + pow(pow.sd, pow) + boxCox(lambda) | cond
    }), df))
  })

  test_that("error when adding distributions that do not support the additive notation", {
    expect_error(.errProcessExpression(quote({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + log(wt / 70) * cl.wt + sex * cl.sex + age * cl.age + 3)
      v  <- exp(tv + eta.v + wt * v.wt + sex * v.sex + age * v.age + 2)
      vp <- exp(tvp + wt * vp.wt + sex * vp.sex + age * vp.age)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl/v * center
      cp = center/v
      cp ~ dbinom(add.sd, pow) + pow(pow.sd, pow) + boxCox(lambda) | cond
    }), df))
  })

  test_that("error for specifying distributions that have multiple numbers of arguments", {
    expect_error(.errProcessExpression(quote({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + log(wt / 70) * cl.wt + sex * cl.sex + age * cl.age + 3)
      v  <- exp(tv + eta.v + wt * v.wt + sex * v.sex + age * v.age + 2)
      vp <- exp(tvp + wt * vp.wt + sex * vp.sex + age * vp.age)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl/v * center
      cp = center/v
      cp ~ dbinom(add.sd, pow)
      center ~ pow(pow.sd, pow) + boxCox(lambda) | cond
    }), df))
  })

  test_that("error when adding algebraic expressions to known distributional abbreviations", {
    expect_error(.errProcessExpression(quote({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + log(wt / 70) * cl.wt + sex * cl.sex + age * cl.age + 3)
      v  <- exp(tv + eta.v + wt * v.wt + sex * v.sex + age * v.age + 2)
      vp <- exp(tvp + wt * vp.wt + sex * vp.sex + age * vp.age)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl/v * center
      cp = center/v
      cp ~ add(add.sd) + pow(pow.sd, pow) + boxCox(lambda) + tan(vp)| cond
    }), df))
  })

  test_that("The distribution names will transform to the preferred distributions", {
    expect_equal(rxPreferredDistributionName("dnorm"), "add")
    expect_equal(rxPreferredDistributionName("add"), "add")
    expect_equal(rxPreferredDistributionName("logitNorm"), "logitNorm")
  })


  test_that("non-numeric bounds for logitNorm and probitNorm", {

    expect_error(.errProcessExpression(quote({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + log(wt / 70) * cl.wt + sex * cl.sex + age * cl.age + 3)
      v  <- exp(tv + eta.v + wt * v.wt + sex * v.sex + age * v.age + 2)
      vp <- exp(tvp + wt * vp.wt + sex * vp.sex + age * vp.age)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl/v * center
      cp = center/v
      cp ~ logitNorm(add.sd, lower, upper) + pow(pow.sd, pow) + boxCox(lambda) | cond
    }), df))


    expect_error(.errProcessExpression(quote({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + log(wt / 70) * cl.wt + sex * cl.sex + age * cl.age + 3)
      v  <- exp(tv + eta.v + wt * v.wt + sex * v.sex + age * v.age + 2)
      vp <- exp(tvp + wt * vp.wt + sex * vp.sex + age * vp.age)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl/v * center
      cp = center/v
      cp ~ logitNorm(add.sd) + pow(pow.sd, pow) + boxCox(lambda) | cond
    }), df), NA)

    expect_error(.errProcessExpression(quote({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + log(wt / 70) * cl.wt + sex * cl.sex + age * cl.age + 3)
      v  <- exp(tv + eta.v + wt * v.wt + sex * v.sex + age * v.age + 2)
      vp <- exp(tvp + wt * vp.wt + sex * vp.sex + age * vp.age)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl/v * center
      cp = center/v
      cp ~ logitNorm(add.sd, 0) + pow(pow.sd, pow) + boxCox(lambda) | cond
    }), df), NA)

    expect_error(.errProcessExpression(quote({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + log(wt / 70) * cl.wt + sex * cl.sex + age * cl.age + 3)
      v  <- exp(tv + eta.v + wt * v.wt + sex * v.sex + age * v.age + 2)
      vp <- exp(tvp + wt * vp.wt + sex * vp.sex + age * vp.age)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl/v * center
      cp = center/v
      cp ~ logitNorm(add.sd, 0, 5) + pow(pow.sd, pow) + boxCox(lambda) | cond
    }), df), NA)


    .errProcessExpression(quote({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + log(wt / 70) * cl.wt + sex * cl.sex + age * cl.age + 3)
      v  <- exp(tv + eta.v + wt * v.wt + sex * v.sex + age * v.age + 2)
      vp <- exp(tvp + wt * vp.wt + sex * vp.sex + age * vp.age)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl/v * center
      cp = center/v
      cp ~ logitNorm(add.sd, -5, 5) + pow(pow.sd, pow) + boxCox(lambda) | cond
    }), df)


    .errProcessExpression(quote({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + log(wt / 70) * cl.wt + sex * cl.sex + age * cl.age + 3)
      v  <- exp(tv + eta.v + wt * v.wt + sex * v.sex + age * v.age + 2)
      vp <- exp(tvp + wt * vp.wt + sex * vp.sex + age * vp.age)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl/v * center
      cp = center/v
      cp ~ logitNorm(add.sd, -Inf, Inf) + pow(pow.sd, pow) + boxCox(lambda) | cond
    }), df)

    expect_error(.errProcessExpression(quote({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + log(wt / 70) * cl.wt + sex * cl.sex + age * cl.age + 3)
      v  <- exp(tv + eta.v + wt * v.wt + sex * v.sex + age * v.age + 2)
      vp <- exp(tvp + wt * vp.wt + sex * vp.sex + age * vp.age)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl/v * center
      cp = center/v
      cp ~ probitNorm(add.sd, lower, upper) + pow(pow.sd, pow) + boxCox(lambda) | cond
    }), df))

  })

  test_that("rxTransformCombine", {

    testCombine0 <- function(transforms) {
      .cmb <- rxDistributionCombine("matt")
      for (i in transforms) {
        .cmb <- rxDistributionCombine(.cmb, i)
      }
      .cmb
    }

    testCombine <- function(transforms) {
      .cmb <- testCombine0(transforms)
      .cmb2 <- testCombine0(rev(transforms))
      expect_equal(.cmb, .cmb2)
      .cmb
    }

    ## add + prop
    expect_equal(testCombine(c("add", "prop"))$errType,
                 testCombine(c("add", "propT"))$errType)

    expect_equal(testCombine(c("add", "prop"))$errType,
                 testCombine(c("add", "propF"))$errType)

    expect_equal(testCombine(c("logitNorm", "prop"))$errType,
                 testCombine(c("add", "prop"))$errType)

    expect_equal(testCombine(c("logitNorm", "propT"))$errType,
                 testCombine(c("add", "prop"))$errType)

    expect_equal(testCombine(c("logitNorm", "propF"))$errType,
                 testCombine(c("add", "prop"))$errType)

    expect_equal(testCombine(c("probitNorm", "prop"))$errType,
                 testCombine(c("add", "prop"))$errType)

    expect_equal(testCombine(c("probitNorm", "propT"))$errType,
                 testCombine(c("add", "prop"))$errType)

    expect_equal(testCombine(c("probitNorm", "propF"))$errType,
                 testCombine(c("add", "prop"))$errType)

    ## add + pow
    expect_equal(testCombine(c("add", "pow"))$errType,
                 testCombine(c("add", "powF"))$errType)

    expect_equal(testCombine(c("lnorm", "pow"))$errType,
                 testCombine(c("add", "pow"))$errType)

    expect_equal(testCombine(c("logitNorm", "pow"))$errType,
                 testCombine(c("add", "pow"))$errType)

    expect_equal(testCombine(c("probitNorm", "pow"))$errType,
                 testCombine(c("add", "pow"))$errType)

    # Test Error type F
    expect_equal(testCombine(c("add", "powF"))$errTypeF,
                 testCombine(c("add", "propF"))$errTypeF)

    expect_equal(testCombine(c("add", "powT"))$errTypeF,
                 testCombine(c("add", "propT"))$errTypeF)

    expect_equal(testCombine(c("add", "pow"))$errTypeF,
                 testCombine(c("add", "prop"))$errTypeF)

    expect_equal(testCombine("powF")$errTypeF,
                 testCombine("propF")$errTypeF)

    expect_equal(testCombine("powT")$errTypeF,
                 testCombine("propT")$errTypeF)

    expect_equal(testCombine("pow")$errTypeF,
                 testCombine("prop")$errTypeF)

    # Test Transformation Type
    expect_equal(testCombine(c("add", "prop", "boxCox"))$transform,
                 testCombine(c("add", "propT", "boxCox"))$transform)

    expect_equal(testCombine(c("add", "prop", "boxCox"))$transform,
                 testCombine(c("add", "propF", "boxCox"))$transform)

    expect_equal(testCombine(c("add", "pow", "boxCox"))$transform,
                 testCombine(c("add", "propF", "boxCox"))$transform)

    expect_equal(testCombine(c("pow", "boxCox"))$transform,
                 testCombine(c("add", "boxCox"))$transform)

    expect_equal(testCombine(c("prop", "boxCox"))$transform,
                 testCombine(c("add", "boxCox"))$transform)


    expect_equal(testCombine(c("lnorm", "prop", "boxCox")),
                 testCombine(c("lnorm", "propT", "boxCox")))

    expect_equal(testCombine(c("lnorm", "prop", "boxCox")),
                 testCombine(c("lnorm", "propF", "boxCox")))

    expect_equal(testCombine(c("lnorm", "pow", "boxCox")),
                 testCombine(c("lnorm", "propF", "boxCox")))

    expect_equal(testCombine(c("logitNorm", "prop"))$transform,
                 testCombine(c("logitNorm", "propT"))$transform)

    expect_equal(testCombine(c("logitNorm", "add"))$transform,
                 testCombine(c("logitNorm", "propT"))$transform)

    expect_equal(testCombine(c("logitNorm", "prop", "boxCox"))$transform,
                 testCombine(c("logitNorm", "propT", "boxCox"))$transform)

    expect_equal(testCombine(c("logitNorm", "prop", "boxCox"))$transform,
                 testCombine(c("logitNorm", "propF", "boxCox"))$transform)

    expect_equal(testCombine(c("logitNorm", "pow", "boxCox"))$transform,
                 testCombine(c("logitNorm", "propF", "boxCox"))$transform)

    expect_equal(testCombine(c("probitNorm", "prop", "boxCox"))$transform,
                 testCombine(c("probitNorm", "propT", "boxCox"))$transform)

    expect_equal(testCombine(c("probitNorm", "prop", "boxCox"))$transform,
                 testCombine(c("probitNorm", "propF", "boxCox"))$transform)

    expect_equal(testCombine(c("probitNorm", "pow", "boxCox"))$transform,
                 testCombine(c("probitNorm", "propF", "boxCox"))$transform)

    # Yeo Johnson
    expect_equal(testCombine(c("add", "prop", "yeoJohnson"))$transform,
                 testCombine(c("add", "propT", "yeoJohnson"))$transform)

    expect_equal(testCombine(c("add", "prop", "yeoJohnson"))$transform,
                 testCombine(c("add", "propF", "yeoJohnson"))$transform)

    expect_equal(testCombine(c("add", "pow", "yeoJohnson"))$transform,
                 testCombine(c("add", "propF", "yeoJohnson"))$transform)

    expect_equal(testCombine(c("pow", "yeoJohnson"))$transform,
                 testCombine(c("add", "yeoJohnson"))$transform)

    expect_equal(testCombine(c("prop", "yeoJohnson"))$transform,
                 testCombine(c("add", "yeoJohnson"))$transform)

    expect_equal(testCombine(c("lnorm", "prop", "yeoJohnson")),
                 testCombine(c("lnorm", "propT", "yeoJohnson")))

    expect_equal(testCombine(c("lnorm", "prop", "yeoJohnson")),
                 testCombine(c("lnorm", "propF", "yeoJohnson")))

    expect_equal(testCombine(c("lnorm", "pow", "yeoJohnson")),
                 testCombine(c("lnorm", "propF", "yeoJohnson")))

    expect_equal(testCombine(c("logitNorm", "prop", "yeoJohnson"))$transform,
                 testCombine(c("logitNorm", "propT", "yeoJohnson"))$transform)

    expect_equal(testCombine(c("logitNorm", "prop", "yeoJohnson"))$transform,
                 testCombine(c("logitNorm", "propF", "yeoJohnson"))$transform)

    expect_equal(testCombine(c("logitNorm", "pow", "yeoJohnson"))$transform,
                 testCombine(c("logitNorm", "propF", "yeoJohnson"))$transform)

    expect_equal(testCombine(c("probitNorm", "prop", "yeoJohnson"))$transform,
                 testCombine(c("probitNorm", "propT", "yeoJohnson"))$transform)

    expect_equal(testCombine(c("probitNorm", "prop", "yeoJohnson"))$transform,
                 testCombine(c("probitNorm", "propF", "yeoJohnson"))$transform)

    expect_equal(testCombine(c("probitNorm", "pow", "yeoJohnson"))$transform,
                 testCombine(c("probitNorm", "propF", "yeoJohnson"))$transform)

    # Box Cox and Yeo Johnson are not compatible
    expect_equal(testCombine(c("boxCox", "prop", "yeoJohnson")),
                 testCombine(c("boxCox", "propT", "yeoJohnson")))

    # Prop and pow are not compatible
    expect_equal(testCombine(c("pow", "prop")),
                 testCombine(c("pow", "prop")))

    # Make sure demotion works
    expect_equal(rxDemoteAddErr(testCombine(c("add", "prop"))),
                 testCombine("prop"))

    expect_equal(rxDemoteAddErr(testCombine(c("add", "pow"))),
                 testCombine("pow"))



  })


  .errProcessExpression(quote({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl + log(wt / 70) * cl.wt + sex * cl.sex + age * cl.age + 3)
    v  <- exp(tv + eta.v + wt * v.wt + sex * v.sex + age * v.age + 2)
    vp <- exp(tvp + wt * vp.wt + sex * vp.sex + age * vp.age)
    d/dt(depot) = -ka * depot
    d/dt(center) = ka * depot - cl/v * center
    cp = center/v
    cp ~ add(add.sd) + pow(pow.sd, pow) + boxCox(lambda) | cond
  }), df)

 }, test="lvl2")
