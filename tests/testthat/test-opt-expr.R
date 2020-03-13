rxPermissive({
    test_that("simple expression optimization", {
    context("Expression optimization tests")

    exp1 <- "rx_yj_~2;\nrx_lambda_~1;\nrx_pred_=10*exp(-THETA[1]*t*exp(ETA[1]));\nrx__sens_rx_pred__BY_ETA_1___=-10*THETA[1]*t*exp(ETA[1])*exp(-THETA[1]*t*exp(ETA[1]));\nrx_r_=100*Rx_pow_di(THETA[2],2)*exp(-2*THETA[1]*t);\ndvid(3,4)\n"

    expect_equal(rxOptExpr(exp1), "rx_yj_~2\nrx_lambda_~1\nrx_expr_0~exp(ETA[1])\nrx_expr_1~exp(-THETA[1]*t*rx_expr_0)\nrx_pred_=10*rx_expr_1\nrx__sens_rx_pred__BY_ETA_1___=-10*THETA[1]*t*rx_expr_0*rx_expr_1\nrx_r_=100*Rx_pow_di(THETA[2], 2)*exp(-2*THETA[1]*t)\ndvid(3, 4)")
    })
}, test="parsing")

