rxPermissive({
  context("param order")

  mod <- RxODE({
    a <- 6
    b <- 0.6
    cmt(blood) # cmt = 1 now
    d / dt(intestine) <- -a * intestine
    d / dt(blood) <- a * intestine - b * blood
  })

  expect_equal(rxModelVars(mod)$param, c("a", "b"))

  mod2 <- RxODE({
    param(b, a)
    a <- 6
    b <- 0.6
    cmt(blood) # cmt = 1 now
    d / dt(intestine) <- -a * intestine
    d / dt(blood) <- a * intestine - b * blood
  })


})
