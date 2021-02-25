rxodeTest(
  {
    context("phi/pnorm/qnorm")

    expect_equal(phi(1:3), pnorm(1:3))
    expect_equal(phi(as.double(1:3)), pnorm(as.double(1:3)))

    o <- RxODE({
      o <- phi(a)
    })

    expect_equal(
      rxSolve(o, data.frame(a = 1:3), et(0))$o,
      pnorm(as.double(1:3))
    )

    o <- RxODE({
      o <- pnorm(a)
    })

    expect_equal(
      rxSolve(o, data.frame(a = 1:3), et(0))$o,
      pnorm(as.double(1:3))
    )

    o <- RxODE({
      o <- pnorm(a, 0.5)
    })

    expect_equal(
      rxSolve(o, data.frame(a = 1:3), et(0))$o,
      pnorm(as.double(1:3), 0.5)
    )

    o <- RxODE({
      o <- pnorm(a, 0.5, 2)
    })

    expect_equal(
      rxSolve(o, data.frame(a = 1:3), et(0))$o,
      pnorm(as.double(1:3), 0.5, 2)
    )

    expect_error(RxODE({
      o <- pnorm()
    }))

    expect_error(RxODE({
      o <- pnorm(a, b, c, d)
    }))

    o <- RxODE({
      o <- qnorm(a)
    })

    expect_equal(
      rxSolve(o, data.frame(a = 1:3), et(0))$o,
      suppressWarnings(qnorm(as.double(1:3)))
    )

    o <- RxODE({
      o <- qnorm(a, 0.5)
    })

    expect_equal(
      rxSolve(o, data.frame(a = 1:3), et(0))$o,
      suppressWarnings(qnorm(as.double(1:3), 0.5))
    )

    o <- RxODE({
      o <- qnorm(a, 0.5, 2)
    })

    expect_equal(
      rxSolve(o, data.frame(a = 1:3), et(0))$o,
      suppressWarnings(qnorm(as.double(1:3), 0.5, 2))
    )

    expect_error(RxODE({
      o <- qnorm()
    }))

    expect_error(RxODE({
      o <- qnorm(a, b, c, d)
    }))

    m <- RxODE({
      o <- pnorm(a)
    })

    if (requireNamespace("units", quietly = TRUE)) {

      expect_error(rxS(m), NA)

      m <- RxODE({
        o <- pnorm(a, b)
      })

      expect_error(rxS(m), NA)

      m <- RxODE({
        o <- pnorm(a, b, c)
      })

      expect_error(rxS(m), NA)
    }
  },
  test = "lvl2"
)
