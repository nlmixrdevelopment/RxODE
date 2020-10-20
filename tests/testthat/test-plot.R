rxPermissive(
  {
    library(units)
    context("plot tests")

    test_that("plot tests", {
      skip_if(utils::packageVersion("ggplot2") < "3.3.0")

      ## Model from RxODE tutorial
      m1 <- RxODE({
        KA <- 2.94E-01
        CL <- 1.86E+01
        V2 <- 4.02E+01
        Q <- 1.05E+01
        V3 <- 2.97E+02
        Kin <- 1
        Kout <- 1
        EC50 <- 200
        ## Added modeled bioavaiblity, duration and rate
        fdepot <- 1
        durDepot <- 8
        rateDepot <- 1250
        C2 <- centr / V2
        C3 <- peri / V3
        d/dt(depot) <- -KA * depot
        f(depot) <- fdepot
        dur(depot) <- durDepot
        rate(depot) <- rateDepot
        d/dt(centr) <- KA * depot - CL * C2 - Q * C2 + Q * C3
        d/dt(peri) <- Q * C2 - Q * C3
        d/dt(eff) <- Kin - Kout * (1 - C2 / (EC50 + C2)) * eff
        eff(0) <- 1
      })

      ev <- et(timeUnits = "hr") %>%
        et(amt = 10000, ii = 12, until = 24) %>%
        et(seq(0, 24, length.out = 100))

      evR <- et(timeUnits="hr") %>% et(evid=3)

      evR <- dplyr::bind_rows(as.data.frame(ev), as.data.frame(evR), as.data.frame(ev))

      evR4 <- rbind(data.frame(id=1, evR), data.frame(id=2, evR), data.frame(id=3, evR), data.frame(id=4, evR))

      s <- rxSolve(m1, ev)

      sR <- rxSolve(m1, evR)

      ev <- et(timeUnits = "hr") %>%
        et(amt = 10000, until = set_units(3, days), ii = 12) %>% # loading doses
        et(seq(0, 48, length.out = 200)) %>%
        et(id = 1:4)

      set.seed(32)
      s2 <- rxSolve(m1, ev, params = data.frame(KA = 0.294 * exp(rnorm(4)), 18.6 * exp(rnorm(4))))

      set.seed(32)
      s2R <- rxSolve(m1, evR4, params = data.frame(KA = 0.294 * exp(rnorm(4)), 18.6 * exp(rnorm(4))))

      set.seed(32)
      s20 <- rxSolve(m1, ev, params = data.frame(KA = 0.294 * exp(rnorm(20)), 18.6 * exp(rnorm(20))))

      set.seed(32)
      s20R <- rxSolve(m1, evR4, params = data.frame(KA = 0.294 * exp(rnorm(20)), 18.6 * exp(rnorm(20))))

      m2 <- RxODE({
        KA <- 2.94E-01
        CL <- 1.86E+01 * exp(eta.Cl)
        V2 <- 4.02E+01
        Q <- 1.05E+01
        V3 <- 2.97E+02
        Kin <- 1
        Kout <- 1
        EC50 <- 200
        ## Added modeled bioavaiblity, duration and rate
        fdepot <- 1
        durDepot <- 8
        rateDepot <- 1250
        C2 <- centr / V2
        C3 <- peri / V3
        d / dt(depot) <- -KA * depot
        f(depot) <- fdepot
        dur(depot) <- durDepot
        rate(depot) <- rateDepot
        d / dt(centr) <- KA * depot - CL * C2 - Q * C2 + Q * C3
        d / dt(peri) <- Q * C2 - Q * C3
        d / dt(eff) <- Kin - Kout * (1 - C2 / (EC50 + C2)) * eff
        eff(0) <- 1
      })

      omega <- lotri(eta.Cl ~ 0.4^2)

      ev <- et(amount.units = "mg", time.units = "hours") %>%
        et(amt = 10000, cmt = "centr", until = 48, ii = 8) %>%
        et(0, 48, length.out = 100)

      evR <- et(amount.units = "mg", time.units = "hours") %>% et(evid=3)

      evR <- dplyr::bind_rows(as.data.frame(ev), as.data.frame(evR), as.data.frame(ev))

      sim <- rxSolve(m2, ev, omega = omega, nSub = 100)

      simR <- rxSolve(m2, evR, omega = omega, nSub = 100)

      expect_error(plot(sim, log = "f"))

      ev2 <- et() %>%
        et(amt = 10000, cmt = "centr", until = 48, ii = 8) %>%
        et(0, 48, length.out = 100)

      ev2R <- et(evid=3)

      ev2R <- dplyr::bind_rows(as.data.frame(ev2), as.data.frame(ev2R), as.data.frame(ev2))

      sim3 <- rxSolve(m2, ev2, omega = omega, nSub = 3)
      sim3R <- rxSolve(m2, ev2R, omega = omega, nSub = 3)

      vdiffr::expect_doppelganger("sim.id-unitless", plot(sim3, C2))

      options(RxODE.theme = FALSE)
      vdiffr::expect_doppelganger("sim.id-unitless-notheme", plot(sim3, C2))
      options(RxODE.theme = TRUE)

      ci1.C2 <- confint(sim, c("C2"))

      ci1.C2.eff <- confint(sim, c("C2", "eff"))

      sim2 <- rxSolve(m2, ev, omega = omega, nSub = 2500)
      sim2R <- rxSolve(m2, evR, omega = omega, nSub = 2500)

      ci2.C2 <- confint(sim2, c("C2"))

      ci2.C2.eff <- confint(sim2, c("C2", "eff"))

      f <- function(xgxr = FALSE, repel = FALSE) {
        if (xgxr) {
          .xgxtxt <- "xgxr-"
          options(RxODE.xgxr = TRUE)
        } else {
          .xgxtxt <- ""
          options(RxODE.xgxr = FALSE)
        }

        vdiffr::expect_doppelganger(paste0("plot-", .xgxtxt, "C2"), s %>% plot(C2))
        vdiffr::expect_doppelganger(paste0("plot-", .xgxtxt, "C2-log-x"), s %>% plot(C2, log = "x"))
        vdiffr::expect_doppelganger(paste0("plot-", .xgxtxt, "C2-log-y"), s %>% plot(C2, log = "y"))
        vdiffr::expect_doppelganger(paste0("plot-", .xgxtxt, "C2-log-xy"), s %>% plot(C2, log = "xy"))
        vdiffr::expect_doppelganger(paste0("plot-", .xgxtxt, "C2-log-yx"), s %>% plot(C2, log = "yx"))
        vdiffr::expect_doppelganger(paste0("plot-", .xgxtxt, "all"), s %>% plot())
        vdiffr::expect_doppelganger(paste0("plot-", .xgxtxt, "all-log-x"), s %>% plot(log = "x"))
        vdiffr::expect_doppelganger(paste0("plot-", .xgxtxt, "all-log-y"), s %>% plot(log = "y"))
        vdiffr::expect_doppelganger(paste0("plot-", .xgxtxt, "all-log-xy"), s %>% plot(log = "xy"))
        vdiffr::expect_doppelganger(paste0("plot-", .xgxtxt, "all-log-yx"), s %>% plot(log = "yx"))

        vdiffr::expect_doppelganger(paste0("plot-", .xgxtxt, "C2-r"), sR %>% plot(C2))
        vdiffr::expect_doppelganger(paste0("plot-", .xgxtxt, "C2-log-x-r"), sR %>% plot(C2, log = "x"))
        vdiffr::expect_doppelganger(paste0("plot-", .xgxtxt, "C2-log-y-r"), sR %>% plot(C2, log = "y"))
        vdiffr::expect_doppelganger(paste0("plot-", .xgxtxt, "C2-log-xy-r"), sR %>% plot(C2, log = "xy"))
        vdiffr::expect_doppelganger(paste0("plot-", .xgxtxt, "C2-log-yx-r"), sR %>% plot(C2, log = "yx"))
        vdiffr::expect_doppelganger(paste0("plot-", .xgxtxt, "all-r"), sR %>% plot())
        vdiffr::expect_doppelganger(paste0("plot-", .xgxtxt, "all-log-x-r"), sR %>% plot(log = "x"))
        vdiffr::expect_doppelganger(paste0("plot-", .xgxtxt, "all-log-y-r"), sR %>% plot(log = "y"))
        vdiffr::expect_doppelganger(paste0("plot-", .xgxtxt, "all-log-xy-r"), sR %>% plot(log = "xy"))
        vdiffr::expect_doppelganger(paste0("plot-", .xgxtxt, "all-log-yx-r"), sR %>% plot(log = "yx"))

        vdiffr::expect_doppelganger(paste0("plot-ci1c2", .xgxtxt), ci1.C2 %>% plot())
        vdiffr::expect_doppelganger(paste0("plot-ci1c2", .xgxtxt, "log-x"), ci1.C2 %>% plot(log = "x"))
        vdiffr::expect_doppelganger(paste0("plot-ci1c2", .xgxtxt, "log-y"), ci1.C2 %>% plot(log = "y"))
        vdiffr::expect_doppelganger(paste0("plot-ci1c2", .xgxtxt, "log-xy"), ci1.C2 %>% plot(log = "xy"))
        vdiffr::expect_doppelganger(paste0("plot-ci1c2", .xgxtxt, "log-yx"), ci1.C2 %>% plot(log = "yx"))

        vdiffr::expect_doppelganger(paste0("plot-ci1c2-eff", .xgxtxt), ci1.C2.eff %>% plot())
        vdiffr::expect_doppelganger(paste0("plot-ci1c2-eff", .xgxtxt, "log-x"), ci1.C2.eff %>% plot(log = "x"))
        vdiffr::expect_doppelganger(paste0("plot-ci1c2-eff", .xgxtxt, "log-y"), ci1.C2.eff %>% plot(log = "y"))
        vdiffr::expect_doppelganger(paste0("plot-ci1c2-eff", .xgxtxt, "log-xy"), ci1.C2.eff %>% plot(log = "xy"))
        vdiffr::expect_doppelganger(paste0("plot-ci1c2-eff", .xgxtxt, "log-yx"), ci1.C2.eff %>% plot(log = "yx"))

        vdiffr::expect_doppelganger(paste0("plot-ci2c2", .xgxtxt), ci2.C2 %>% plot())
        vdiffr::expect_doppelganger(paste0("plot-ci2c2", .xgxtxt, "log-x"), ci2.C2 %>% plot(log = "x"))
        vdiffr::expect_doppelganger(paste0("plot-ci2c2", .xgxtxt, "log-y"), ci2.C2 %>% plot(log = "y"))
        vdiffr::expect_doppelganger(paste0("plot-ci2c2", .xgxtxt, "log-xy"), ci2.C2 %>% plot(log = "xy"))
        vdiffr::expect_doppelganger(paste0("plot-ci2c2", .xgxtxt, "log-yx"), ci2.C2 %>% plot(log = "yx"))

        vdiffr::expect_doppelganger(paste0("plot-ci2c2-eff", .xgxtxt), ci2.C2.eff %>% plot())
        vdiffr::expect_doppelganger(paste0("plot-ci2c2-eff", .xgxtxt, "log-x"), ci2.C2.eff %>% plot(log = "x"))
        vdiffr::expect_doppelganger(paste0("plot-ci2c2-eff", .xgxtxt, "log-y"), ci2.C2.eff %>% plot(log = "y"))
        vdiffr::expect_doppelganger(paste0("plot-ci2c2-eff", .xgxtxt, "log-xy"), ci2.C2.eff %>% plot(log = "xy"))
        vdiffr::expect_doppelganger(paste0("plot-ci2c2-eff", .xgxtxt, "log-yx"), ci2.C2.eff %>% plot(log = "yx"))

        ##

        ## large
        vdiffr::expect_doppelganger(paste0("plot-sp-", .xgxtxt, "C2"), s20 %>% plot(C2))
        vdiffr::expect_doppelganger(paste0("plot-sp-", .xgxtxt, "C2-log-x"), s20 %>% plot(C2, log = "x"))
        vdiffr::expect_doppelganger(paste0("plot-sp-", .xgxtxt, "C2-log-y"), s20 %>% plot(C2, log = "y"))
        vdiffr::expect_doppelganger(paste0("plot-sp-", .xgxtxt, "C2-log-xy"), s20 %>% plot(C2, log = "xy"))
        vdiffr::expect_doppelganger(paste0("plot-sp-", .xgxtxt, "C2-log-yx"), s20 %>% plot(C2, log = "yx"))
        vdiffr::expect_doppelganger(paste0("plot-sp-", .xgxtxt, "all"), s20 %>% plot())
        vdiffr::expect_doppelganger(paste0("plot-sp-", .xgxtxt, "all-log-x"), s20 %>% plot(log = "x"))
        vdiffr::expect_doppelganger(paste0("plot-sp-", .xgxtxt, "all-log-y"), s20 %>% plot(log = "y"))
        vdiffr::expect_doppelganger(paste0("plot-sp-", .xgxtxt, "all-log-xy"), s20 %>% plot(log = "xy"))
        vdiffr::expect_doppelganger(paste0("plot-sp-", .xgxtxt, "all-log-yx"), s20 %>% plot(log = "yx"))

        vdiffr::expect_doppelganger(paste0("plot-sp-", .xgxtxt, "C2-r"), s20R %>% plot(C2))
        vdiffr::expect_doppelganger(paste0("plot-sp-", .xgxtxt, "C2-log-x-r"), s20R %>% plot(C2, log = "x"))
        vdiffr::expect_doppelganger(paste0("plot-sp-", .xgxtxt, "C2-log-y-r"), s20R %>% plot(C2, log = "y"))
        vdiffr::expect_doppelganger(paste0("plot-sp-", .xgxtxt, "C2-log-xy-r"), s20R %>% plot(C2, log = "xy"))
        vdiffr::expect_doppelganger(paste0("plot-sp-", .xgxtxt, "C2-log-yx-r"), s20R %>% plot(C2, log = "yx"))
        vdiffr::expect_doppelganger(paste0("plot-sp-", .xgxtxt, "all-r"), s20R %>% plot())
        vdiffr::expect_doppelganger(paste0("plot-sp-", .xgxtxt, "all-log-x-r"), s20R %>% plot(log = "x"))
        vdiffr::expect_doppelganger(paste0("plot-sp-", .xgxtxt, "all-log-y-r"), s20R %>% plot(log = "y"))
        vdiffr::expect_doppelganger(paste0("plot-sp-", .xgxtxt, "all-log-xy-r"), s20R %>% plot(log = "xy"))
        vdiffr::expect_doppelganger(paste0("plot-sp-", .xgxtxt, "all-log-yx-r"), s20R %>% plot(log = "yx"))

        for (repel in c(TRUE, FALSE)) {
          if (repel) {
            .repel <- "repel-"
            options(RxODE.ggrepel = TRUE)
          } else {
            .repel <- ""
            options(RxODE.ggrepel = FALSE)
          }

          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "C2"), s2 %>% plot(C2))
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "C2-log-x"), s2 %>% plot(C2, log = "x"))
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "C2-log-y"), s2 %>% plot(C2, log = "y"))
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "C2-log-xy"), s2 %>% plot(C2, log = "xy"))
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "C2-log-yx"), s2 %>% plot(C2, log = "yx"))
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "all"), s2 %>% plot())
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "all-log-x"), s2 %>% plot(log = "x"))
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "all-log-y"), s2 %>% plot(log = "y"))
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "all-log-xy"), s2 %>% plot(log = "xy"))
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "all-log-yx"), s2 %>% plot(log = "yx"))

          ## Issue #284
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "284-log-yx"), s2 %>% plot(C2, eff, log = "yx"))
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "284-log-y"), s2 %>% plot(C2, eff, log = "y"))
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "284-log-x"), s2 %>% plot(C2, eff, log = "x"))

          # reset
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "C2-r"), s2R %>% plot(C2))
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "C2-log-x-r"), s2R %>% plot(C2, log = "x"))
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "C2-log-y-r"), s2R %>% plot(C2, log = "y"))
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "C2-log-xy-r"), s2R %>% plot(C2, log = "xy"))
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "C2-log-yx-r"), s2R %>% plot(C2, log = "yx"))
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "all-r"), s2R %>% plot())
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "all-log-x-r"), s2R %>% plot(log = "x"))
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "all-log-y-r"), s2R %>% plot(log = "y"))
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "all-log-xy-r"), s2R %>% plot(log = "xy"))
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "all-log-yx-r"), s2R %>% plot(log = "yx"))

          ## Issue #284
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "284-log-yx-r"), s2R %>% plot(C2, eff, log = "yx"))
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "284-log-y-r"), s2R %>% plot(C2, eff, log = "y"))
          vdiffr::expect_doppelganger(paste0("plot-multi-", .repel, .xgxtxt, "284-log-x-r"), s2R %>% plot(C2, eff, log = "x"))
        }
      }
      f(TRUE)
      f(FALSE)
    })

    options(RxODE.xgxr = TRUE)
    options(RxODE.ggrepel = TRUE)
    options(RxODE.theme = TRUE)
  },
  silent = TRUE,
  test = "plot"
)
