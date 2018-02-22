context("Test modelvars");
rxPermissive({

    rigid.txt <- "
y1(0)    = 1
y2(0)    = 0
y3(0)    = 0.9
a1       = -2
a2       = 1.25
a3       = -0.5
d/dt(y1) = a1*y2*y3
d/dt(y2) = a2*y1*y3
d/dt(y3) = a3*y1*y2
";
    rigid <- RxODE(rigid.txt)

    et <- eventTable();
    et$add.sampling(seq(0,20,by=0.01))

    out <- solve(rigid,et)


    test_that("modelvars", {
        expect_equal(rxModelVars(rigid), rxModelVars(rigid$cmpMgr$rxDll()))
        expect_equal(rxModelVars(rigid), rxModelVars(out))
    })
}, silent=TRUE);
