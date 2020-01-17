## Prior FOCEi
rxPermissive({

    context("focei cross-check with old sympy RxODE");
    ## Test Wang 2007 focei earlier
    test_that("wang 2007 focei", {

        ## SymPy RxODE cross-check
        rx <- RxODE({
            rx_expr_001 ~ -THETA[1]
            rx_expr_002 ~ rx_expr_001 * t
            rx_expr_003 ~ exp(ETA[1])
            rx_expr_004 ~ -2 * THETA[1]
            rx_expr_005 ~ rx_expr_004 * t
            rx_expr_006 ~ Rx_pow_di(THETA[2], 2)
            rx_expr_007 ~ rx_expr_002 * rx_expr_003
            rx_expr_008 ~ rx_expr_005 * rx_expr_003
            rx_expr_009 ~ exp(rx_expr_007)
            rx_expr_010 ~ exp(rx_expr_008)
            rx_yj_ ~ 2
            rx_lambda_ ~ 1
            rx_pred_ = 10 * rx_expr_009
            rx__sens_rx_pred__BY_ETA_1___ = -10 * THETA[1] * t * rx_expr_003 *
                rx_expr_009
            rx_r_ = 100 * rx_expr_006 * rx_expr_010
            rx__sens_rx_r__BY_ETA_1___ = -200 * THETA[1] * rx_expr_006 *
                t * rx_expr_003 * rx_expr_010
        })


        .pars <- data.frame("ID"=1:10, "THETA[1]" = 0.5, "THETA[2]" = 0.316227766016838,
                            "ETA[1]" = c(0.0715450734486017, 0.00570452906524344,  0.0249471848462488, 0.0319502378032186, -0.00478736811351454,  0.0039789224128819, -0.0118035084999257, 0.00588012198220198,  -0.00313671907426165, -0.000666414340548085),
                            check.names = FALSE)

        .dat <- data.frame(ID = c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L,  6L, 6L, 7L, 7L, 8L, 8L, 9L, 9L, 10L, 10L),
                           Time = c(0L, 1L, 0L,  1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L,  1L),
                           DV = c(10.68, 3.6837, 10.402, 6.454, 9.8814, 5.8565, 9.3408,  5.6209, 10.082, 6.7583, 9.8938, 6.5049, 9.8908, 6.9557, 10.234,  6.4488, 9.9882, 6.7112, 9.6736, 6.6402))

        solve1 <- rxSolve(rx, .dat, .pars, returnType="data.frame")


        mypar1 = function ()
        {
            ke = theta[1] * exp(eta[1]);
        }

        mypar2 = function ()
        {
            k = theta[1] * exp(eta[1]);
            v = 1
        }

        mod <- RxODE({
            ipre = 10 * exp(-ke * t)
        })

        pred = function() ipre

        err2 <- function(f){
            return(theta[2]); ## SD
        }

        m2d <- rxSymPySetupPred(mod, pred, mypar1,
                                function(){return(prop(.1))})

        solve2 <- rxSolve(m2d$inner, .dat, .pars, returnType="data.frame")

        expect_equal(solve1, solve2)

        m2d <- rxSymPySetupPred(mod, pred, mypar1,
                                function(){return(prop(.1))},
                                optExpression = FALSE)

        solve2 <- rxSolve(m2d$inner, .dat, .pars, returnType="data.frame")

        expect_equal(solve1, solve2)
    })

    context("foce cross-check with old sympy RxODE");

})
