rxPermissive({
    ode <- RxODE({
                  C2 = centr/V2;
                  C3 = peri/V3;
                  d/dt(depot) =-KA*depot;
                  d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
                  d/dt(peri)  =                    Q*C2 - Q*C3;
                  d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
    });
    test_that("concentrations and volumes are correct", {
        expect_equal(rxModelVars(ode)$conc, c(centr="C2", peri="C3"))
        expect_equal(rxModelVars(ode)$vol, c(centr="V2", peri="V3"))
    });
})
