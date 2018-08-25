rxPermissive({

    context("Test Issue #50")
    test_that("Issue #50", {

        sink("rxodedemo_ode.txt")
        cat("C2 = centr/V2;
C3 = peri/V3;
d/dt(depot) =-KAdepot;
d/dt(centr) = KAdepot - CLC2 - QC2 + QC3;
d/dt(peri) = QC2 - QC3;
d/dt(eff) = Kin - Kout(1-C2/(EC50+C2))*eff;")
    sink();

    expect_true(file.exists("rxodedemo_ode.txt"))
    m1 <- try(RxODE(filename = 'rxodedemo_ode.txt',modName = 'm1'))
    expect_false(is(class(m1), "RxODE"))
    expect_true(file.exists("rxodedemo_ode.txt"))

    sink("rxodedemo_ode.txt")
    cat("C2 = centr/V2;
C3 = peri/V3;
d/dt(depot) =-KA*depot;
d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
d/dt(peri)  =                    Q*C2 - Q*C3;
d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;")
    sink();

    expect_true(file.exists("rxodedemo_ode.txt"))
    m1 <- RxODE(filename = 'rxodedemo_ode.txt',modName = 'm1')
    expect_true(is(m1, "RxODE"))
    expect_true(file.exists("rxodedemo_ode.txt"))

    rxDelete(m1)
    unlink("rxodedemo_ode.txt")
    if (dir.exists("m1.d")) unlink("m1.d", recursive=TRUE)
})
}, silent=TRUE, on.validate=TRUE)
