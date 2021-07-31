rxodeTest({
  test_that("Parallel solve", {

    TV_CLr=6.54 # L/h, (CLr/F)
    TV_CLnr=2.39 # L/h, (CLnr/F)
    TV_Vc=95.1 # L, (V/F)
    TV_alag=0.145 # h,
    TV_D=0.512 # h,
    TV_Q=2.1 #L/h, (Q/F)
    TV_Vp=23.3 #L, (Vp/F)

    OM_D_normal = log((128/100)^2+1)
    D_trans = 0.0819
    OM_CLr = log((36.2/100)^2+1)
    OM_CLnr = log((43.6/100)^2+1)
    OM_Vc = log((14.4/100)^2+1)
    OM_Q = log((15.1/100)^2+1)
    OM_Vp = log((37.6/100)^2+1)
    OM_CLr_CLnr = 0.101
    OM_CLr_Vc = 0.0066

    rwishart <- function(df, p = nrow(SqrtSigma), Sigma, SqrtSigma = diag(p)) {
      if (!missing(Sigma)) {
        tmp <- svd(Sigma)
        SqrtSigma <- sqrt(tmp$d) * t(tmp$u)
      }
      if ((Ident <- missing(SqrtSigma)) && missing(p))
        stop("either p, Sigma or SqrtSigma must be specified")
      Z <- matrix(0, p, p)
      diag(Z) <- sqrt(rchisq(p, df:(df - p + 1)))
      if (p > 1) {
        pseq <- 1:(p - 1)
        Z[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p *
                                                                      (p - 1)/2)
      }
      if (Ident)
        crossprod(Z)
      else crossprod(Z %*% SqrtSigma)
    }

    sample.etas <- function(df, R, n) {
      R.inv = solve(df*R)
      R0 = chol(.5*(R.inv+t(R.inv)))  ## R0 is an upper-diag matrix
      R1 = solve(rwishart(df, SqrtSigma=R0))
      eta = MASS::mvrnorm(n = n,
                          mu= rep(0, nrow(R)),
                          Sigma= .5*(R1+t(R1)))
      eta
    }

    set.seed(100)
    nsubj <- 200

    n.pts.pk <- 300  # number of subjects in PK model fitting

    eta_CLr_CLnr_Vc <- as.data.frame(sample.etas(df = n.pts.pk,
                                                R = lotri({eta_CLr + eta_CLnr + eta_Vc ~
                                                             c(OM_CLr,
                                                               OM_CLr_CLnr, OM_CLnr,
                                                               OM_CLr_Vc, 0, OM_Vc)
                                                }),
                                                n = nsubj))

    eta_D_trans <- data.frame(eta_D_normal = rnorm(mean=0,sd=sqrt(OM_D_normal),n=nsubj)) %>%
      dplyr::mutate(eta_D = ((exp(eta_D_normal))^D_trans-1)/D_trans)

    par.pk <- data.frame(sim.id = seq(nsubj),
                        D = TV_D * exp(eta_D_trans$eta_D),
                        CLr = TV_CLr * exp(eta_CLr_CLnr_Vc$eta_CLr),
                        CLnr = TV_CLnr * exp(eta_CLr_CLnr_Vc$eta_CLnr),
                        Vc = TV_Vc * exp(eta_CLr_CLnr_Vc$eta_Vc),
                        Vp = TV_Vp * exp(rnorm(nsubj, 0, sqrt(OM_Vp))),
                        Q = TV_Q * exp(rnorm(nsubj, 0, sqrt(OM_Q))),
                        alag = TV_alag)

    # zero-order absorption with lag time 2-compartment
    mod <- RxODE({
      CL <- CLr+CLnr
      C2 <- central/Vc*1000
      all<- central+periph+output

      d/dt(central) <- - CL/Vc*central - Q/Vc*central + Q/Vp*periph
      d/dt(periph)  <- Q/Vc*central - Q/Vp*periph
      d/dt(output)  <- CL/Vc*central
      alag(central) <- alag
      dur(central) <- D
    })

    ev <-  et(amt=2,cmt="central",rate=-2,ii=24,addl=4) %>%
      et(seq(0,120,0.1))

    bar2x <- rxSolve(mod, ev, params=par.pk, cores=2L, returnType="data.frame")

    bar1x <- rxSolve(mod, ev, params=par.pk, cores=1L, returnType="data.frame")

    expect_equal(bar1x, bar2x)

  })
}, test="lvl2")
