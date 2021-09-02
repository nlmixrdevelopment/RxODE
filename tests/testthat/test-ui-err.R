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
