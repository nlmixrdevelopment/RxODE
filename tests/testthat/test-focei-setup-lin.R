 context("Linear Model Focei Setup checks");
rxPermissive({

    pred <- function () {return(Central)}

    pk <- function ()
    {
        lCl = THETA[1]
        lVc = THETA[2]
        prop.err = THETA[3]
        eta.Vc = ETA[1]
        eta.Cl = ETA[2]
        Vc <- exp(lVc + eta.Vc)
        Cl <- exp(lCl + eta.Cl)
    }

    err <- function ()
    {
        return(prop(prop.err))
    }

    mod <- RxODE({
        rx_ka ~ 0
        rx_tlag ~ 0
        rx_v ~ Vc
        rx_k ~ Cl/Vc
        rx_alpha ~ rx_k
        rx_A ~ 1.0 / rx_v
        rx_beta ~ 0
        rx_B ~ 0
        rx_gamma ~ 0
        rx_C ~ 0
        Central=solveLinB(rx__PTR__, t, 0, 0, 0, rx_A, rx_alpha, rx_B, rx_beta, rx_C, rx_gamma, rx_ka, rx_tlag);
    })

    ## 1 compartment model
    pk1 <- rxSymPySetupPred(mod, predfn=pred, pkpars=pk, err=err)

    pk2 <- rxSymPySetupPred(mod, predfn=pred, pkpars=pk, err=err, grad=TRUE)

    test_that("One compartment model", {
        expect_equal(class(pk1), "rxFocei")
        expect_equal(class(pk2), "rxFocei")
    })

    pk <- function()
    {
        lCl = THETA[1]
        lVc = THETA[2]
        lKA = THETA[3]
        prop.err = THETA[4]
        eta.Cl = ETA[1]
        eta.Vc = ETA[2]
        eta.KA = ETA[3]
        Cl <- exp(lCl + eta.Cl)
        Vc <- exp(lVc + eta.Vc)
        KA <- exp(lKA + eta.KA)
    }

    mod <- RxODE({
        rx_ka ~ KA
        rx_tlag ~ 0
        rx_v ~ Vc
        rx_k ~ Cl/Vc
        rx_alpha ~ rx_k
        rx_A ~ rx_ka / (rx_ka - rx_alpha) / rx_v
        rx_beta ~ 0
        rx_B ~ 0
        rx_gamma ~ 0
        rx_C ~ 0
        Central=solveLinB(rx__PTR__, t, 0, 0, 0, rx_A, rx_alpha, rx_B, rx_beta, rx_C, rx_gamma, rx_ka, rx_tlag);
    })

    ## 1 compartment oral
    pk1 <- rxSymPySetupPred(mod, predfn=pred, pkpars=pk, err=err)
    pk2 <- rxSymPySetupPred(mod, predfn=pred, pkpars=pk, err=err, grad=TRUE)

    test_that("One compartment oral model", {
        expect_equal(class(pk1), "rxFocei")
        expect_equal(class(pk2), "rxFocei")
    })

    mod <- RxODE({
        rx_ka ~ 0
        rx_tlag ~ 0
        rx_v ~ Vc
        rx_k ~ Cl/Vc
        rx_k12 ~ Q/Vc
        rx_k21 ~ Q/Vp
        rx_beta  ~ 0.5 * (rx_k12 + rx_k21 + rx_k - sqrt((rx_k12 + rx_k21 + rx_k) * (rx_k12 + rx_k21 + rx_k) - 4.0 * rx_k21 * rx_k))
        rx_alpha ~ rx_k21 * rx_k / rx_beta
        rx_A ~ (rx_alpha - rx_k21) / (rx_alpha - rx_beta) / rx_v
        rx_B ~ (rx_beta - rx_k21) / (rx_beta - rx_alpha) / rx_v;
        rx_gamma ~ 0
        rx_C ~ 0
        Central=solveLinB(rx__PTR__, t, 0, 0, 0, rx_A, rx_alpha, rx_B, rx_beta, rx_C, rx_gamma, rx_ka, rx_tlag);
    })

    pk <- function () {
        lCl = THETA[1]
        lVc = THETA[2]
        lQ = THETA[3]
        lVp = THETA[4]
        prop.err = THETA[5]
        eta.Vc = ETA[1]
        eta.Cl = ETA[2]
        eta.Vp = ETA[3]
        eta.Q = ETA[4]
        Vc <- exp(lVc + eta.Vc)
        Cl <- exp(lCl + eta.Cl)
        Vp <- exp(lVp + eta.Vp)
        Q <- exp(lQ + eta.Q)
    }

    pk1 <- rxSymPySetupPred(mod, predfn=pred, pkpars=pk, err=err)

    pk2 <- rxSymPySetupPred(mod, predfn=pred, pkpars=pk, err=err, grad=TRUE)

    test_that("Two compartment model", {
        expect_equal(class(pk1), "rxFocei")
        expect_equal(class(pk2), "rxFocei")
    })

    pk <- function ()
    {
        lCl = THETA[1]
        lVc = THETA[2]
        lQ = THETA[3]
        lVp = THETA[4]
        lKA = THETA[5]
        prop.err = THETA[6]
        eta.Vc = ETA[1]
        eta.Cl = ETA[2]
        eta.Vp = ETA[3]
        eta.Q = ETA[4]
        eta.KA = ETA[5]
        Vc <- exp(lVc + eta.Vc)
        Cl <- exp(lCl + eta.Cl)
        Vp <- exp(lVp + eta.Vp)
        Q <- exp(lQ + eta.Q)
        KA <- exp(lKA + eta.KA)
    }

    mod <- RxODE({
        rx_ka ~ KA
        rx_tlag ~ 0
        rx_v ~ Vc
        rx_k ~ Cl/Vc
        rx_k12 ~ Q/Vc
        rx_k21 ~ Q/Vp
        rx_beta  ~ 0.5 * (rx_k12 + rx_k21 + rx_k - sqrt((rx_k12 + rx_k21 + rx_k) * (rx_k12 + rx_k21 + rx_k) - 4.0 * rx_k21 * rx_k))
        rx_alpha ~ rx_k21 * rx_k / rx_beta
        rx_A ~ rx_ka / (rx_ka - rx_alpha) * (rx_alpha - rx_k21) / (rx_alpha - rx_beta) / rx_v
        rx_B ~ rx_ka / (rx_ka - rx_beta) * (rx_beta - rx_k21) / (rx_beta - rx_alpha) / rx_v;
        rx_gamma ~ 0
        rx_C ~ 0
        Central=solveLinB(rx__PTR__, t, 0, 0, 0, rx_A, rx_alpha, rx_B, rx_beta, rx_C, rx_gamma, rx_ka, rx_tlag);
    })

    pk1 <- rxSymPySetupPred(mod, predfn=pred, pkpars=pk, err=err)

    pk2 <- rxSymPySetupPred(mod, predfn=pred, pkpars=pk, err=err, grad=TRUE)

    test_that("Two compartment oral model", {
        expect_equal(class(pk1), "rxFocei")
        expect_equal(class(pk2), "rxFocei")
    })

    pk <- function () {
        lCl = THETA[1]
        lVc = THETA[2]
        lQ = THETA[3]
        lVp = THETA[4]
        lV3 = THETA[5]
        lQ2 = THETA[6]
        prop.err = THETA[7]
        eta.Vc = ETA[1]
        eta.Cl = ETA[2]
        eta.Vp = ETA[3]
        eta.Q = ETA[4]
        eta.V3 = ETA[5]
        eta.Q2 = ETA[6]
        Vc <- exp(lVc + eta.Vc)
        Cl <- exp(lCl + eta.Cl)
        Vp <- exp(lVp + eta.Vp)
        Q <- exp(lQ + eta.Q)
        V3 <- exp(lV3 + eta.V3)
        Q2 <- exp(lQ2 + eta.Q2)
    }

    mod <- RxODE({
        rx_tlag ~ 0
        rx_ka ~ 0
        rx_v ~ Vc
        rx_k ~ Cl/Vc
        rx_k12 ~ Q/Vc
        rx_k21 ~ Q/Vp
        rx_k13 ~ Q2/Vc
        rx_k31 ~ Q2/V2
        rx_a0 ~ rx_k * rx_k21 * rx_k31
        rx_a1 ~ rx_k * rx_k31 + rx_k21 * rx_k31 + rx_k21 * rx_k13 + rx_k * rx_k21 + rx_k31 * rx_k12
        rx_a2 ~ rx_k + rx_k12 + rx_k13 + rx_k21 + rx_k3
        rx_p ~ rx_a1 - rx_a2 * rx_a2 / 3.0
        rx_q ~ 2.0 * rx_a2 * rx_a2 * rx_a2 / 27.0 - rx_a1 * rx_a2 /3.0 + rx_a
        rx_r1 ~ sqrt(-rx_p * rx_p * rx_p / 27.0)
        rx_r2 ~ 2 * rx_r1^(1.0/3.0)
        rx_theta ~ acos(-rx_q / (2.0 * rx_r1)) / 3
        rx_alpha ~ -(cos(rx_theta) * rx_r2 - rx_a2 / 3.0)
        rx_beta ~ -(cos(rx_theta + 2.0 / 3.0 * pi) * rx_r2 - rx_a2 / 3.0)
        rx_gamma ~ -(cos(rx_theta + 4.0 / 3.0 * pi) * rx_r2 - rx_a2 / 3.0)
        rx_A ~ (rx_k21 - rx_alpha) * (rx_k31 - rx_alpha) / (rx_alpha - rx_beta) / (rx_alpha - rx_gamma) / rx_v
        rx_B ~ (rx_k21 - rx_beta) * (rx_k31 - rx_beta) / (rx_beta - rx_alpha) / (rx_beta - rx_gamma) / rx_v
        rx_C ~ (rx_k21 - rx_gamma) * (rx_k31 - rx_gamma) / (rx_gamma - rx_alpha) / (rx_gamma - rx_beta) / rx_v
        Central=solveLinB(rx__PTR__, t, 0, 0, 0, rx_A, rx_alpha, rx_B, rx_beta, rx_C, rx_gamma, rx_ka, rx_tlag);
    })

    pk1 <- rxSymPySetupPred(mod, predfn=pred, pkpars=pk, err=err)

    ## pk2 <- rxSymPySetupPred(mod, predfn=pred, pkpars=pk, err=err, grad=TRUE)  # 231 cmt model

    test_that("Three compartment model", {
        expect_equal(class(pk1), "rxFocei")
        ## expect_equal(class(pk2), "rxFocei")
    })

    pk <- function () {
        lCl = THETA[1]
        lVc = THETA[2]
        lQ = THETA[3]
        lVp = THETA[4]
        lV3 = THETA[5]
        lQ2 = THETA[6]
        lKa = THETA[7]
        prop.err = THETA[8]
        eta.Vc = ETA[1]
        eta.Cl = ETA[2]
        eta.Vp = ETA[3]
        eta.Q = ETA[4]
        eta.V3 = ETA[5]
        eta.Q2 = ETA[6]
        eta.Ka = ETA[7]
        Vc <- exp(lVc + eta.Vc)
        Cl <- exp(lCl + eta.Cl)
        Vp <- exp(lVp + eta.Vp)
        Q <- exp(lQ + eta.Q)
        V3 <- exp(lV3 + eta.V3)
        Q2 <- exp(lQ2 + eta.Q2)
        Ka <- exp(lKa + eta.Ka)
    }

    mod <- RxODE({
        rx_tlag ~ 0
        rx_ka ~ Ka
        rx_v ~ Vc
        rx_k ~ Cl/Vc
        rx_k12 ~ Q/Vc
        rx_k21 ~ Q/Vp
        rx_k13 ~ Q2/Vc
        rx_k31 ~ Q2/V2
        rx_a0 ~ rx_k * rx_k21 * rx_k31
        rx_a1 ~ rx_k * rx_k31 + rx_k21 * rx_k31 + rx_k21 * rx_k13 + rx_k * rx_k21 + rx_k31 * rx_k12
        rx_a2 ~ rx_k + rx_k12 + rx_k13 + rx_k21 + rx_k3
        rx_p ~ rx_a1 - rx_a2 * rx_a2 / 3.0
        rx_q ~ 2.0 * rx_a2 * rx_a2 * rx_a2 / 27.0 - rx_a1 * rx_a2 /3.0 + rx_a
        rx_r1 ~ sqrt(-rx_p * rx_p * rx_p / 27.0)
        rx_r2 ~ 2 * rx_r1^(1.0/3.0)
        rx_theta ~ acos(-rx_q / (2.0 * rx_r1)) / 3
        rx_alpha ~ -(cos(rx_theta) * rx_r2 - rx_a2 / 3.0)
        rx_beta ~ -(cos(rx_theta + 2.0 / 3.0 * pi) * rx_r2 - rx_a2 / 3.0)
        rx_gamma ~ -(cos(rx_theta + 4.0 / 3.0 * pi) * rx_r2 - rx_a2 / 3.0)
        rx_A ~ (rx_k21 - rx_alpha) * (rx_k31 - rx_alpha) / (rx_alpha - rx_beta) / (rx_alpha - rx_gamma) / rx_v
        rx_B ~ (rx_k21 - rx_beta) * (rx_k31 - rx_beta) / (rx_beta - rx_alpha) / (rx_beta - rx_gamma) / rx_v
        rx_C ~ (rx_k21 - rx_gamma) * (rx_k31 - rx_gamma) / (rx_gamma - rx_alpha) / (rx_gamma - rx_beta) / rx_v
        Central=solveLinB(rx__PTR__, t, 0, 0, 0, rx_A, rx_alpha, rx_B, rx_beta, rx_C, rx_gamma, rx_ka, rx_tlag);
    })

    pk1 <- rxSymPySetupPred(mod, predfn=pred, pkpars=pk, err=err)

    ## pk2 <- rxSymPySetupPred(mod, predfn=pred, pkpars=pk, err=err, grad=TRUE);

    test_that("Three compartment oral model", {
        expect_equal(class(pk1), "rxFocei")
        ## expect_equal(class(pk2), "rxFocei")
    })

},  silent=TRUE, on.validate=TRUE)
