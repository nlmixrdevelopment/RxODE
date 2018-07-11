library(RxODE)
library(testthat)
rxPermissive({
    context("Focei Inner");

    d <- as.data.frame(list(time = c(72, 95.99, 96, 119.99, 120, 143.99, 144,
                                     144.25, 144.5, 144.75, 145, 145.5, 146, 146.5, 147, 148, 150,
                                     152, 156, 160, 164, 167.99, 168, 191.99, 192, 215.99, 216, 216.25,
                                     216.5, 216.75, 217, 217.5, 218, 218.5, 219, 220, 222, 224, 228,
                                     232, 236, 240, 252, 264, 276, 288),
                            evid = c(101, 0, 101, 0,
                                     101, 0, 101, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 101,
                                     0, 101, 0, 101, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0),
                            amt = c(60000, NA, 60000, NA, 60000, NA, 60000,
                                    NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 60000,
                                    NA, 60000, NA, 60000, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                    NA, NA, NA, NA, NA, NA, NA, NA, NA),
                            dv = c(NA, 263.6, NA, 164.7,
                                   NA, 287.3, NA, 1248.7, 1211.5, 1017.7, 1690.1, 1029.8, 890.7,
                                   598.4, 1009.3, 1159.8, 742.2, 724.6, 728.2, 509.7, 243.1, 259.9,
                                   NA, 242.2, NA, 281.4, NA, 1500.1, 1281.4, 1200.2, 1378.8, 1373.2,
                                   582.9, 960.2, 720.3, 852.6, 950.3, 654.7, 402.5, 456, 346.5,
                                   268.2, 134.2, 42.6, 25.9, 14.6)))
    d$id <- 1;


    inner <- RxODE("d/dt(centr)~-centr*exp(ETA[1]+THETA[1])*exp(-ETA[2]-THETA[2]);\nd/dt(rx__sens_centr_BY_ETA_1___)~-centr*exp(ETA[1]+THETA[1])*exp(-ETA[2]-THETA[2])-rx__sens_centr_BY_ETA_1___*exp(ETA[1]+THETA[1])*exp(-ETA[2]-THETA[2]);\nd/dt(rx__sens_centr_BY_ETA_2___)~centr*exp(ETA[1]+THETA[1])*exp(-ETA[2]-THETA[2])-rx__sens_centr_BY_ETA_2___*exp(ETA[1]+THETA[1])*exp(-ETA[2]-THETA[2]);\nrx_pred_=centr*exp(-ETA[2]-THETA[2]);\nrx__sens_rx_pred__BY_ETA_1___=rx__sens_centr_BY_ETA_1___*exp(-ETA[2]-THETA[2]);\nrx__sens_rx_pred__BY_ETA_2___=-centr*exp(-ETA[2]-THETA[2])+rx__sens_centr_BY_ETA_2___*exp(-ETA[2]-THETA[2]);\nrx_r_=Rx_pow_di(THETA[3],2)*Rx_pow_di(centr,2)*exp(-2*ETA[2]-2*THETA[2]);\nrx__sens_rx_r__BY_ETA_1___=2*Rx_pow_di(THETA[3],2)*centr*rx__sens_centr_BY_ETA_1___*exp(-2*ETA[2]-2*THETA[2]);\nrx__sens_rx_r__BY_ETA_2___=-2*Rx_pow_di(THETA[3],2)*Rx_pow_di(centr,2)*exp(-2*ETA[2]-2*THETA[2])+2*Rx_pow_di(THETA[3],2)*centr*rx__sens_centr_BY_ETA_2___*exp(-2*ETA[2]-2*THETA[2]);\n");

    theta <- structure(c(1.6, 4.5, sqrt(0.1)), .Names=sprintf("THETA[%d]", 1:3));
    eta <- structure(c(-0.147736086922763, -0.294637022436797), .Names=sprintf("ETA[%d]", 1:2))


    rxInv <- rxSymInvCholCreate(matrix(c(0.1, 0, 0, 0.1), nrow=2))

    etaMat <- t(matrix(eta))

    ## rxSolve setupOnly needs to be called first.
    print(RxODE:::.foceiSetup(inner, d, theta, rxInv = rxInv, etaMat = etaMat,
                              odeOpt=foceiControl(maxOuterIterations=0, maxInnerIterations=0)))

    expect_equal(RxODE:::foceiLik(rep(1, 5)), -209.467627968204)

    expect_equal(RxODE:::foceiInnerLp(rep(0, 2)), c(31.4177868473072, 56.8293635831368))
    expect_equal(RxODE:::foceiInnerLp(eta), c(-0.000367366009438985, -2.46513261989989e-05))

    rxSolveFree()
}, on.validate="NLMIXR_VALIDATION", silent=TRUE)
