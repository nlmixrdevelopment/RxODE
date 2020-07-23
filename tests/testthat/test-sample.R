rxPermissive({
  context("resample tests")
  test_that("resample tests", {

    m1 <- RxODE({
      CL ~ (1-0.2*SEX)*(0.807+0.00514*(CRCL-91.2))*exp(eta.cl)
      V1 ~ 4.8*exp(eta.v1)
      Q ~ (3.46+0.0593*(WT-75.1))*exp(eta.q);
      V2 ~ 1.93*(3.13+0.0458*(WT-75.1))*exp(eta.v2)
      cp = max(linCmt() + err.sd, 0.01)
      if (cp == 0.01) cp = NA
      mWT = WT
      mCRCL = CRCL
    })

    library(dplyr)
    nsub=30
    # Simulate Weight based on age and gender
    AGE<-round(runif(nsub,min=18,max=70))
    SEX<-round(runif(nsub,min=0,max=1))
    HTm<-round(rnorm(nsub,176.3,0.17*sqrt(4482)),digits=1)
    HTf<-round(rnorm(nsub,162.2,0.16*sqrt(4857)),digits=1)
    WTm<-round(exp(3.28+1.92*log(HTm/100))*exp(rnorm(nsub,0,0.14)),digits=1)
    WTf<-round(exp(3.49+1.45*log(HTf/100))*exp(rnorm(nsub,0,0.17)),digits=1)
    WT<-ifelse(SEX==1,WTf,WTm)
    CRCL<-round(runif(nsub,30,140))
    ## id is in lower case to match the event table
    cov.df <- tibble(id=seq_along(AGE), AGE=AGE, SEX=SEX, WT=WT, CRCL=CRCL)

    s<-c(0,0.25,0.5,0.75,1,1.5,seq(2,24,by=1))
    ## Add 10% diff
    s <- lapply(s, function(x){d <- x * 0.1; c(x - d, x + d)})

    e <- et(time.units="hr") %>%
      ## Specify the id and weight based dosing from covariate data.frame
      ## This requires RxODE XXX
      et(id=cov.df$id, amt=6*cov.df$WT, rate=6 * cov.df$WT) %>%
      ## Sampling is added for each ID
      et(s) %>%
      as.data.frame %>%
      ## Merge the event table with the covarite information
      merge(cov.df, by="id") %>%
      as_tibble

    f1 <- rxSolve(m1, e,
                 ## Lotri uses lower-triangular matrix rep. for named matrix
                 omega=lotri(eta.cl ~ .306,
                             eta.q ~0.0652,
                             eta.v1 ~.567,
                             eta.v2 ~ .191),
                 sigma=lotri(err.sd ~ 0.5), addCov = TRUE)

    f2 <- rxSolve(m1, e,
                 ## Lotri uses lower-triangular matrix rep. for named matrix
                 omega=lotri(eta.cl ~ .306,
                             eta.q ~0.0652,
                             eta.v1 ~.567,
                             eta.v2 ~ .191),
                 sigma=lotri(err.sd ~ 0.5), addCov = TRUE,
                 resample=c("SEX", "WT", "CRCL"))

    r1 <- f1[!duplicated(f1$id), c("id", "SEX", "WT", "CRCL")]
    r2 <- f2[!duplicated(f2$id), c("id", "SEX", "WT", "CRCL")]

    expect_false(isTRUE(all.equal(r1, r2)))

    # Make these time-varying covariates

    e$WT <- e$WT + rnorm(length(e$WT), sd=1)
    e$CRCL <- e$CRCL + rnorm(length(e$CRCL), sd=1)

    f1 <- rxSolve(m1, e,
                 ## Lotri uses lower-triangular matrix rep. for named matrix
                 omega=lotri(eta.cl ~ .306,
                             eta.q ~0.0652,
                             eta.v1 ~.567,
                             eta.v2 ~ .191),
                 sigma=lotri(err.sd ~ 0.5), addCov = TRUE)

    f2 <- rxSolve(m1, e,
                 ## Lotri uses lower-triangular matrix rep. for named matrix
                 omega=lotri(eta.cl ~ .306,
                             eta.q ~0.0652,
                             eta.v1 ~.567,
                             eta.v2 ~ .191),
                 sigma=lotri(err.sd ~ 0.5), addCov = TRUE,
                 resample=c("SEX", "WT", "CRCL"))

    f3 <- rxSolve(m1, e,
                  omega=lotri(eta.cl ~ .306,
                             eta.q ~0.0652,
                             eta.v1 ~.567,
                             eta.v2 ~ .191),
                 sigma=lotri(err.sd ~ 0.5), keep = c("SEX", "WT", "CRCL"),
                 resample=c("SEX", "WT", "CRCL"))

    r1 <- f1[!duplicated(f1$id), c("id", "SEX", "WT", "CRCL")]
    r2 <- f2[!duplicated(f2$id), c("id", "SEX", "WT", "CRCL")]

    expect_false(isTRUE(all.equal(r1$WT, r2$WT)))

    ## Now test keep case

    r1 <- f1[!duplicated(f1$id), c("id", "SEX", "WT", "CRCL")]
    r3 <- f3[!duplicated(f3$id), c("id", "SEX", "WT", "CRCL")]

    all.equal(r1$WT, r3$WT)


  })

},
silent = TRUE,
test = "lvl2")
