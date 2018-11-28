rxPermissive({
    context("Absolute Working Directory + Model Name (Issue #56)")
    test_that("Issue #56", {
        skip_on_os("windows")
        if (!dir.exists("/tmp/")){
            skip()
        }
        ode <- " d/dt(test) = 0; "
        m1 <- RxODE(model = ode, modName = "m1", wd = "/tmp")
        expect_true(dir.exists("/tmp/m1.d"))
        rxDelete(m1);
        unlink("/tmp/m1.d", recursive=TRUE)
    })
}, silent=TRUE, cran=TRUE)
