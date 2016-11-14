context("Test that model specification can take string, file name or rxode expressions")
rxPermissive({

    ode <- RxODE(model = 'd/dt(y) = r * y * (1 - y/K);');

    ode2 <- RxODE({
        d/dt(y) = r * y * (1 - y/K);
    });

    test_that("string and expression returns the same result", {
        expect_equal(rxModelVars(ode)$md5["parsed_md5"], rxModelVars(ode2)$md5["parsed_md5"])
    })

    sink("temp.rx")
    cat('d/dt(y) = r * y * (1 - y/K);\n');
    sink()
    ode3 <- RxODE('temp.rx');
    unlink("temp.rx");

    test_that("file and string returns the same result", {
        expect_equal(rxModelVars(ode)$md5["parsed_md5"], rxModelVars(ode3)$md5["parsed_md5"])
    })


    ode <- RxODE('
         C2 = centr/V2;
         C3 = peri/V3;
         d/dt(depot) =-KA*depot;
         d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
         d/dt(peri)  =                    Q*C2 - Q*C3;
         d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;');

    ode2 <- RxODE({
        C2 = centr/V2;
        C3 = peri/V3;
        d/dt(depot) =-KA*depot;
        d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
        d/dt(peri)  =                    Q*C2 - Q*C3;
        d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
    })

    sink("temp.rx")
    cat('
         C2 = centr/V2;
         C3 = peri/V3;
         d/dt(depot) =-KA*depot;
         d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
         d/dt(peri)  =                    Q*C2 - Q*C3;
         d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;');
    sink();
    ode3 <- RxODE('temp.rx');
    unlink("temp.rx");

    test_that("string and expression returns the same result", {
        expect_equal(rxModelVars(ode)$md5["parsed_md5"], rxModelVars(ode2)$md5["parsed_md5"])
    });
    test_that("file and string returns the same result", {
        expect_equal(rxModelVars(ode)$md5["parsed_md5"], rxModelVars(ode3)$md5["parsed_md5"])
    });
}, silent=TRUE);
