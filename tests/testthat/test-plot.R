context("Test Plot")
library(RxODE);
library(digest);
rxPermissive({
    ode <- "
   C2 = centr/V2;
   C3 = peri/V3;
   d/dt(depot) =-KA*depot;
   d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
   d/dt(peri)  =                    Q*C2 - Q*C3;
   d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
"
    m2 <- RxODE(model = ode)

    tmp <- RxODE:::igraph(m2$cmpMgr$rxDll())

    sink("test");
    print(tmp)
    sink()
    tst <- readLines("test");
    unlink("test")

    test_that("igraph works",{
        expect_true(any(igraph::edge.attributes(tmp)$label == "Indirect\nEffect (II)"));
    });



    ode <- "
   C2 = centr/V2;
   C3 = peri/V3;
   d/dt(depot) =-KA*depot;
   d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
   d/dt(peri)  =                    Q*C2 - Q*C3;
   d/dt(eff)  = Kin*(1-C2/(EC50+C2)) - Kout*eff;
"
    m1 <- RxODE(model = ode)

    tmp <- RxODE:::igraph(m1$cmpMgr$rxDll());

    test_that("igraph works",{
        expect_true(any(igraph::edge.attributes(tmp)$label == "Indirect\nEffect (I)"));
    })

    ode <- "
   C2 = centr/V2;
   C3 = peri/V3;
   d/dt(depot) =-KA*depot;
   d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
   d/dt(peri)  =                    Q*C2 - Q*C3;
   d/dt(eff)  = Kin*(1+C2/(EC50+C2)) - Kout*eff;
"

    m3 <- RxODE(model = ode)

    tmp <- RxODE:::igraph(m3$cmpMgr$rxDll());

    test_that("igraph works",{
        expect_true(any(igraph::edge.attributes(tmp)$label == "Indirect\nEffect (III)"));
    })


    ode <- "
   C2 = centr/V2;
   C3 = peri/V3;
   d/dt(depot) =-KA*depot;
   d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
   d/dt(peri)  =                    Q*C2 - Q*C3;
   d/dt(eff)  = Kin - Kout*eff*(1+C2/(EC50+C2));
"

    m4 <- RxODE(model = ode)

    tmp <- RxODE:::igraph(m4$cmpMgr$rxDll());

    test_that("igraph works",{
        expect_true(any(igraph::edge.attributes(tmp)$label == "Indirect\nEffect (IV)"));
    })

}, silent=TRUE, on.validate=TRUE);
