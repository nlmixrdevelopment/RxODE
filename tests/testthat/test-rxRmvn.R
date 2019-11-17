context("rxMvn")
## Tests adapted from
## https://github.com/cran/mvtnorm/blob/master/tests/rmvnorm.R
test_that("mvtnormal simulations make sense", {
    set.seed(1)
    m <- 1:3

    s <- diag(1:3)
    s[2,1] <- 1
    s[3,1] <- 2
    s[3,2] <- 3
    s <- s+t(s)

    set.seed(1)

    x <- rxRmvn(20000, m, s, ncores=1)

    expect_equal(m, colMeans(x), tolerance=0.01)
    expect_equal(s, var(x), tolerance=0.1)

    set.seed(1)
    x <- rxRmvn(20000, m, s, ncores=2)

    expect_equal(m, colMeans(x), tolerance=0.01)
    expect_equal(s, var(x), tolerance=0.1)

})

test_that("seeds behave as expected", {
    m <- 1:3

    s <- diag(1:3)
    s[2,1] <- 1
    s[3,1] <- 2
    s[3,2] <- 3
    s <- s+t(s)

    set.seed(1)
    x1 <- rxRmvn(10, m, s, ncores=1)

    set.seed(2)
    x2 <- rxRmvn(10, m, s, ncores=1)
    expect_false(isTRUE(all.equal(x1, x2)))

    set.seed(1)
    x2 <- rxRmvn(10, m, s, ncores=1)
    expect_equal(x1, x2)

    set.seed(1)
    x1 <- rxRmvn(10, m, s, ncores=1)
    set.seed(2)
    x2 <- rxRmvn(10, m, s, ncores=2)
    expect_false(isTRUE(all.equal(x1, x2)))

})


