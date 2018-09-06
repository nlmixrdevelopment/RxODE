rxPermissive({
    context("lhs for different compartments make sense")

    et <- eventTable() %>% add.dosing(dose=3, nbr.doses=6, dosing.interval=8) %>%
        add.sampling(seq(0, 48, length.out=200))

    df <- et$get.sampling()

    df <- rbind(data.frame(et$get.dosing(), CMT=1),
                data.frame(df, CMT = 1),
                data.frame(df, CMT = 2),
                data.frame(df, CMT = 3))

    df <- df[order(df$time, df$CMT), ]

    ode.2c.ka <- RxODE({
        C2 = centr/V;
        C3 = peri/V2;
        d/dt(depot) =-KA*depot;
        d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
        d/dt(peri)  =                    Q*C2 - Q*C3;
        C4 = CMT
        if (CMT == 1){
            prd <- depot
        }
        if (CMT == 2){
            prd <- centr
        }
        if (CMT == 3){
            prd <- peri
        }
    })

    test_that("multi-compartment solves", {
        s <- rxSolve(ode.2c.ka, params=c(V=40, CL=18, V2=297, Q=10, KA= 0.3), df)

        expect_equal(sort(unique(s$C4)), c(1, 2, 3))
        expect_equal(s[s$C4 == 1, "prd"], s[s$C4 == 1, "depot"])
        expect_equal(s[s$C4 == 2, "prd"], s[s$C4 == 2, "centr"])
        expect_equal(s[s$C4 == 3, "prd"], s[s$C4 == 3, "peri"])

    })
})
