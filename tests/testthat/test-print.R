context("Test Printing and summary functions");

ode <- "
   C2 = centr/V2;
   C3 = peri/V3;
   d/dt(depot) =-KA*depot;
   d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
   d/dt(peri)  =                    Q*C2 - Q*C3;
   d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
"

test_that("Print for m1 works correctly",{
    m1 <- RxODE(model = ode)
    sink("test");
    print(m1)
    sink();
    p1 <- readLines("test");
    unlink("test");
    expect_equal(p1,sprintf("RxODE model named \"%s\" (ready to run)",basename(getwd())))
    rxDelete(m1);
    sink("test");
    print(m1);
    sink();
    p1 <- readLines("test");
    unlink("test");
    ## print(p1);
    expect_equal(p1,sprintf("RxODE model named \"%s\" (invalid object, needs to be re-created)",basename(getwd())))
});


