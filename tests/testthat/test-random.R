context("Random generator tests")

f <- RxODE({
    x1 <- rnorm()
    x2 <- rnorm(a)
    x3 <- rnorm(b, c)
})
