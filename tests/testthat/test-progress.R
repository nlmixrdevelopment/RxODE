context("Test progress bar")

test_that("progress_bar", {

    f <- function(abort=FALSE){
        on.exit({rxProgressAbort()});
        rxProgress(100)
        for (i in 1:100) {
            rxTick()
            if (abort && i == 50) stop("here");
        }
        rxProgressStop();
        TRUE
    }

    expect_true(f())
    expect_error(f(TRUE))
})
