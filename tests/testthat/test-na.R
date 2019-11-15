rxPermissive({
    context("NA tests")

    test_that("can assign NA", {
        mod1 <- RxODE("a=NA;\nb=2;\nc=a+b");
        ret <- rxSolve(mod1, et(0), returnType="data.frame")
        expect_true(all(is.na(ret$a)))
        expect_true(all(is.na(ret$b)))
    })


    test_that("logic NA", {
        mod1 <- RxODE("a=NA;\nif(is.na(a)){a=3}b=2;\nc=a+b+exp(3)");

        ret <- rxSolve(mod1, et(0), returnType="data.frame")

        expect_equal(ret$a, as.double(3))
        expect_equal(ret$c, 3 + 2 + exp(3))

        mod1 <- RxODE("a=NA;\nif(!is.na(a)){a=4} else { a=3 }; b=2;\nc=a+b+exp(3)");
    })


})
