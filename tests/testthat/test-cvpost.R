context("Test simulation cvPost functions")

test_that("single cvPost draw makes sense",{
    draw1 <- cvPost(3, matrix(c(1,.3,.3,1),2,2))
    expect_true(inherits(draw1, "matrix"))
    expect_equal(dim(draw1), c(2L,2L))
})

test_that("cvPost of 3 items make sense.",{
    set.seed(42)
    draw3 <- cvPost(3, matrix(c(1,.3,.3,1),2,2), n=3)
    set.seed(42)
    draw3c <- cvPost(3, matrix(c(1,.3,.3,1),2,2), n=3, returnChol=TRUE)
    expect_true(inherits(draw3, "list"))
    expect_true(inherits(draw3c, "list"))
    for (i in seq_along(draw3)){
        expect_equal(dim(draw3[[i]]),  c(2L, 2L))
        expect_equal(dim(draw3c[[i]]), c(2L, 2L))
        expect_equal(chol(draw3[[i]]), draw3c[[i]])
    }
    set.seed(42)
    draw3c <- cvPost(3, chol(matrix(c(1,.3,.3,1),2,2)), n=3, omegaIsChol=TRUE)

    set.seed(42)
    draw3 <- cvPost(3, matrix(c(1,.3,.3,1),2,2), n=3)
    for (i in seq_along(draw3)){
        expect_equal(draw3[[1]], draw3c[[1]])
    }
})
context("rinvchisq")
test_that("rinvchisq produces proper output",{
    expect_equal(length(rinvchisq(3, 4, 1)), 3) ## Scale = 1, degrees of freedom = 4
    expect_equal(length(rinvchisq(3, 4, 1)), 3) ## Scale = 1, degrees of freedom = 4
})
