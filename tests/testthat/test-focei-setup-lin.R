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

},  silent=TRUE, on.validate=TRUE)
