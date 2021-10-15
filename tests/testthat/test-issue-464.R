rxodeTest({
  test_that("Calculating the number of compartments is off #464", {

    rx <- RxODE({
      cmt(SC)
      cmt(CENTRAL)
      rx_expr_0 ~ exp(THETA[2])
      d/dt(SC) = -rx_expr_0 * SC
      rx_expr_1 ~ ETA[1] + THETA[3]
      rx_expr_2 ~ ETA[2] + THETA[4]
      rx_expr_8 ~ rx_expr_1 - (rx_expr_2)
      rx_expr_9 ~ exp(rx_expr_8)
      rx_expr_10 ~ rx_expr_9 * CENTRAL
      d/dt(CENTRAL) = rx_expr_0 * SC - rx_expr_10
      f(SC) = exp(THETA[1])
      d/dt(rx__sens_SC_BY_ETA_1___) = -rx_expr_0 * rx__sens_SC_BY_ETA_1___
      d/dt(rx__sens_CENTRAL_BY_ETA_1___) = rx_expr_0 * rx__sens_SC_BY_ETA_1___ -
        rx_expr_10 - rx_expr_9 * rx__sens_CENTRAL_BY_ETA_1___
      d/dt(rx__sens_SC_BY_ETA_2___) = -rx_expr_0 * rx__sens_SC_BY_ETA_2___
      d/dt(rx__sens_CENTRAL_BY_ETA_2___) = rx_expr_0 * rx__sens_SC_BY_ETA_2___ +
        rx_expr_10 - rx_expr_9 * rx__sens_CENTRAL_BY_ETA_2___
      rx_yj_ ~ 2
      rx_lambda_ ~ 1
      rx_hi_ ~ 1
      rx_low_ ~ 0
      rx_expr_4 ~ exp(-(rx_expr_2))
      rx_expr_5 ~ 1000 * rx_expr_4
      rx_pred_ = rx_expr_5 * CENTRAL
      rx__sens_rx_pred__BY_ETA_1___ = rx_expr_5 * rx__sens_CENTRAL_BY_ETA_1___
      rx__sens_rx_pred__BY_ETA_2___ = -1000 * rx_expr_4 * CENTRAL +
        rx_expr_5 * rx__sens_CENTRAL_BY_ETA_2___
      rx_expr_3 ~ Rx_pow_di(THETA[5], 2)
      rx_expr_7 ~ rx_expr_4 * CENTRAL
      rx_r_ = 1e+06 * Rx_pow_di((rx_expr_7), 2) * rx_expr_3 + Rx_pow_di(THETA[6],
                                                                        2)
      rx_expr_6 ~ 2e+06 * rx_expr_4
      rx__sens_rx_r__BY_ETA_1___ = rx_expr_6 * (rx_expr_7) * rx__sens_CENTRAL_BY_ETA_1___ *
        rx_expr_3
      rx__sens_rx_r__BY_ETA_2___ = -2e+06 * rx_expr_4 * (rx_expr_7) *
        CENTRAL * rx_expr_3 + rx_expr_6 * (rx_expr_7) * rx__sens_CENTRAL_BY_ETA_2___ *
        rx_expr_3
      cmt(cp)
    })

    d <- qs::qread("issue-464.qs")

    theta <- c(-0.6931472, -6.9077553, -3.5065579, 0.4054651, 0.0100000, 1.0000000)

    theta <- setNames(theta, paste0("THETA[", seq_along(theta), "]"))

    eta <- c(0, 0)

    eta <- setNames(eta, paste0("ETA[", seq_along(eta), "]"))

    p <- c(theta, eta)

    expect_error(rxSolve(rx, p, d), NA)

  })
})
