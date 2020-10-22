context("Expression optimization tests")
test_that("simple expression optimization", {
    exp1 <- "rx_yj_~2;\nrx_lambda_~1;\nrx_pred_=10*exp(-THETA[1]*t*exp(ETA[1]));\nrx__sens_rx_pred__BY_ETA_1___=-10*THETA[1]*t*exp(ETA[1])*exp(-THETA[1]*t*exp(ETA[1]));\nrx_r_=100*Rx_pow_di(THETA[2],2)*exp(-2*THETA[1]*t);\n"

    expect_equal(rxOptExpr(exp1), "rx_expr_001 ~ -THETA[1]\nrx_expr_002 ~ rx_expr_001*t\nrx_expr_003 ~ exp(ETA[1])\nrx_expr_004 ~ rx_expr_002*rx_expr_003\nrx_expr_005 ~ exp(rx_expr_004)\nrx_yj_ ~ 2\nrx_lambda_ ~ 1\nrx_pred_ = 10*rx_expr_005\nrx__sens_rx_pred__BY_ETA_1___ = -10*THETA[1]*t*rx_expr_003*rx_expr_005\nrx_r_ = 100*Rx_pow_di(THETA[2], 2)*exp(-2*THETA[1]*t)")
})

