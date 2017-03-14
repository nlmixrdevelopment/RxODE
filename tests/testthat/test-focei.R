context("Focei checks");
rxPermissive({

    mat.indices <- function(nETA){
        idx = do.call("rbind",
                      lapply(1:nETA, function(k) cbind(k:nETA, k)))
        H = matrix(1:(nETA^2), nETA, nETA)
        Hlo.idx = row(H)>=col(H)
        lo.idx = H[row(H)>col(H)]
        hi.idx = t(H)[row(H)>col(H)]

        list(idx=idx,                       # (r, c) of lo-half
             Hlo.idx=Hlo.idx,       # index of lo-half
             lo.idx=lo.idx,         # index of strict lo-half
             hi.idx=hi.idx)         # index of strict hi-half
    }

    m1 <- RxODE({
        d/dt(centr) = - CL / V*centr;
    })

    m1a <- RxODE(m1, calcSens=TRUE);

    test_that("m1a created successfully.", {
        expect_equal(class(m1a), "RxODE");
    })

    m2 <- RxODE({
        d/dt(depot) = -KA*depot;
        d/dt(centr) = KA*depot - CL / V*centr;
    })

    et <- eventTable() %>% add.dosing(dose=30000) %>%
        add.sampling(c(0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 12, 16, 20, 24, 36, 48, 60, 71.99))

    pk <- function(){
        KA = exp(THETA[1])
        CL = exp(THETA[2] + ETA[1])
        V = exp(THETA[3] + ETA[2])
    }

    pred <- function(){
        return(cntr);
    }

    test_that("Error when pred dosen't depend on state varaibles", {
        expect_error(rxSymPySetupPred(m2, pred, pk));
    })

    pred <- function(){
        return(centr);
    }
    ## err ~ prop(.1) + add(.2)

    err <- function(f){
        return(f ^ 2* theta[4] ^ 2); ## Theta 4 is residual sd for proportional error.
    }

    err2 <- function(f){
        return(theta[4]); ## SD
    }

    m2a1 <- rxSymPySetupPred(m2, pred, pk)

    m2a2 <- rxSymPySetupPred(m2, pred)

    m2a <- rxSymPySetupPred(m2, pred, pk, err)

    err2 <- function(f){
        return(prop(0.1));
    }

    m2b <- rxSymPySetupPred(m2, pred, pk, err2)

    err3 <- function(f){
        return(add(0.05) + prop(0.1));
    }

    m2c <- rxSymPySetupPred(m2, pred, pk, err3)

    pred2 <- function(){
        if (cmt == 2){
            return(cntr);
        } else {
            return(depot)
        }
    }

    test_that("A warning should occur if some of the prediction doesn't depend on the pred.", {
        expect_warning(rxSymPySetupPred(m2, pred2, pk, err3))
    })

    pred2 <- function(){
        if (cmt == 2){
            return(centr);
        } else {
            return(depot)
        }
    }

    m2d <- rxSymPySetupPred(m2, pred2, pk, err3)


    ## FIXME should work without a return...
    ## err4 <- function(){
    ##     if (cmt == 2){
    ##         add(0.3);
    ##     } else {
    ##         prop(0.3) + add(0.2);
    ##     }
    ## }

    err4 <- function(){
        if (cmt == 2){
            return(add(0.3));
        } else {
            return(prop(0.3) + add(0.2));
        }
    }

    m2e <- rxSymPySetupPred(m2, pred2, pk, err4)

    err5 <- function() add(0.3);

    m2f <- rxSymPySetupPred(m2, pred2, pk, err5)

    test_that("1, 2 and 3 parameter Pred Setup works", {
        expect_equal(class(m2a1), "rxFocei")
        expect_equal(class(m2a2), "rxFocei")
        expect_equal(class(m2a), "rxFocei")
        expect_equal(class(m2b), "rxFocei")
        expect_equal(class(m2c), "rxFocei")
        expect_equal(class(m2d), "rxFocei")
        expect_equal(class(m2e), "rxFocei")
        expect_equal(class(m2f), "rxFocei")
        expect_true(length(rxInit(m2f$inner)) == 0)
    })


    dv <- c(1126.1, 869.9, 883.6, 1244, 995.2, 946.4, 589.2, 754.4, 1060.8, 433, 609.5, 630.1, 342.1, 322.9, 196.3, 89.8, 42.8, 16.6, 8.4)

    omega <- matrix(c(0.1, 0, 0, 0.1), nrow=2)

    symo <- rxSymInvCreate(omega);

    symenv <- rxSymInv(symo, c(sqrt(0.1), sqrt(0.1)))

    ## test_that("theta/eta solve works", {
    ##     expect_equal(suppressWarnings(digest(m2a %>% solve(et, theta=c(2, 1.6, 4.5, 0.1), eta=c(0.01, -0.01)) %>% as.data.frame %>% round(3))),
    ##                  "13baa6cabb14367e1328149d3c14ebf7");
    ## })

    ## m2a <- rxSymPySetupPred(m2, pred, pk, err)


    tmp1 <- m2a$inner %>% solve(et, theta=c(2, 1.6, 4.5,0.01), eta=c(0.01, -0.01))

    tmp2 <- m2a %>% rxFoceiEta(et, theta=c(2, 1.6, 4.5,0.01), eta=c(0.01, -0.01),dv=dv, inv.env=symenv)

    tmp2.nm <- m2a %>% rxFoceiEta(et, theta=c(2, 1.6, 4.5,0.01), eta=c(0.01, -0.01),dv=dv, inv.env=symenv, nonmem=TRUE)

    tmp3 <- m2a %>% rxFoceiLik(et, theta=c(2, 1.6, 4.5,0.01), eta=c(0.01, -0.01),dv=dv, inv.env=tmp2)

    tmp4 <- m2a %>% rxFoceiLp(et, theta=c(2, 1.6, 4.5,0.01), eta=c(0.01, -0.01),dv=dv, inv.env=symenv)

    tmp5 <- m2a %>%rxFoceiInner(et, theta=c(2, 1.6, 4.5,0.01), eta=c(10, 10),dv=dv, inv.env=symenv, invisible=1)

    test_that("rxFoceiEta makes sense", {
        expect_equal(tmp1$rx_pred_, tmp2$f); ## F
        err <- matrix(dv - tmp1$rx_pred_, ncol=1)
        expect_equal(err, tmp2$err) ## Err
        R <- matrix(tmp1$rx_r_, ncol=1)
        expect_equal(R, tmp2$R) ## R (Varinace)
        m <- as.matrix(tmp1[,c("_sens_rx_pred__ETA_1_", "_sens_rx_pred__ETA_2_")])
        dimnames(m) <- list(NULL, NULL)
        m2 <- tmp2$dErr
        dimnames(m2) <- list(NULL, NULL)
        expect_equal(m, m2) ## Check dErr
        m <- as.matrix(tmp1[, c("_sens_rx_r__ETA_1_", "_sens_rx_r__ETA_2_")]);
        dimnames(m) <- list(NULL, NULL);
        m2 <- tmp2$dR
        dimnames(m2) <- list(NULL, NULL)
        expect_equal(m, m2) ## Check dR
        c <- list(matrix(tmp1[["_sens_rx_r__ETA_1_"]] / tmp1[["rx_r_"]],ncol=1),matrix(tmp1[["_sens_rx_r__ETA_2_"]] / tmp1[["rx_r_"]],ncol=1))
        expect_equal(c, tmp2$c) ## c
        B <- matrix(2 / tmp1[["rx_r_"]],ncol=1)
        expect_equal(B, tmp2$B)
        a <- list(matrix(tmp1[["_sens_rx_pred__ETA_1_"]],ncol=1) - err/R*matrix(tmp1[["_sens_rx_r__ETA_1_"]]),
                  matrix(tmp1[["_sens_rx_pred__ETA_2_"]],ncol=1)- err/R*matrix(tmp1[["_sens_rx_r__ETA_2_"]]));
        expect_equal(a, tmp2$a);

        a <- list(matrix(tmp2.nm$dErr[, 1], ncol=1),
                  matrix(tmp2.nm$dErr[, 2], ncol=1))
        expect_equal(tmp2.nm$a, a)

        ## Does not include the matrix multilpaction part (done with RcppArmadillo)
        lp <- matrix(c(NA, NA), ncol=1)
        c <- matrix(tmp1[["_sens_rx_r__ETA_1_"]] / tmp1[["rx_r_"]],ncol=1)
        fp <- as.matrix(tmp1[,c("_sens_rx_pred__ETA_1_" )])
        lp[1, 1] <- .5*apply(err*fp*B + .5*err^2*B*c - c, 2, sum)
        c <- matrix(tmp1[["_sens_rx_r__ETA_2_"]] / tmp1[["rx_r_"]],ncol=1)
        fp <- as.matrix(tmp1[,c("_sens_rx_pred__ETA_2_" )])
        lp[2, 1] <- .5*apply(err*fp*B + .5*err^2*B*c - c, 2, sum)
        expect_equal(lp, tmp2$lp)

        llik <- -0.5 * sum(err ^ 2 / R + log(R));

        expect_equal(llik, tmp2$llik);

        eta <- c(0.01, -0.01)
        omegaInv <- symenv$omegaInv
        llik <- -0.5 * sum(err ^ 2 / R + log(R)) - 0.5 * t(matrix(eta,ncol=1)) %*% omegaInv %*% matrix(eta,ncol=1);
        llik <- -llik;

        expect_equal(as.vector(llik), tmp3);

        lp2 <- lp - omegaInv %*% matrix(eta, ncol=1)

        expect_equal(-lp2, tmp4);

        ## Now Test the Hessian calculation
        Hidx <- mat.indices(2);
        H <- matrix(1:(2^2), 2, 2)
        k = Hidx$idx[,1];
        l = Hidx$idx[,2];
        a = cbind(tmp2$a[[1]], tmp2$a[[2]])
        B = as.vector(tmp2$B);
        c = cbind(tmp2$c[[1]], tmp2$c[[2]])
        NONMEM <- 0;
        H[Hidx$Hlo.idx] = apply(a[, k]*B*a[, l] + c(-1, 1)[NONMEM+1] *c[,k]*c[,l], 2, sum)
        H[Hidx$hi.idx] = H[Hidx$lo.idx]
        H = -.5*H - omegaInv;

        rxHessian(tmp2);
        expect_equal(tmp2$H, H);

        H.neg.5 = tryCatch({
            chol(-H)
        }, error = function(e) {
            cat("Warning: Hessian not positive definite\n")
            print(-H)
            .m <- -H
            .md = matrix(0, nETA, nETA)
            diag(.md) = abs(diag(.m)) * 1.1 + .001
            chol(.md)
        })
                                        #H.neg.5 = chol(-H)
        log.det.H.neg.5 = sum(log(diag(H.neg.5)))

        RxODE_focei_eta_lik(tmp2$eta, tmp2)
        log.det.OMGAinv.5 <- symenv$log.det.OMGAinv.5
        llik = -tmp2$llik2 + log.det.OMGAinv.5             #note no -1/2 with OMGAinv.5
        llik.lapl =llik - log.det.H.neg.5               #note no 1/2
        llik.lapl = as.vector(llik.lapl);
        attr(llik.lapl, "fitted") <- tmp2$f;
        attr(llik.lapl, "posthoc") <- tmp2$eta;

        lik2 <- RxODE_focei_finalize_llik(tmp2);

        expect_equal(lik2, llik.lapl);

        ## Now test NONMEM approximation.
        NONMEM <- 1;
        Hidx <- mat.indices(2);
        H <- matrix(1:(2^2), 2, 2)
        a = cbind(tmp2.nm$a[[1]], tmp2.nm$a[[2]])
        B = as.vector(tmp2.nm$B);
        c = cbind(tmp2.nm$c[[1]], tmp2.nm$c[[2]])
        H[Hidx$Hlo.idx] = apply(a[, k]*B*a[, l] + c(-1, 1)[NONMEM+1] *c[,k]*c[,l], 2, sum)
        H[Hidx$hi.idx] = H[Hidx$lo.idx]
        H = -.5*H - omegaInv;

        rxHessian(tmp2.nm);
        expect_equal(tmp2.nm$H, H);

        H.neg.5 = tryCatch({
            chol(-H)
        }, error = function(e) {
            cat("Warning: Hessian not positive definite\n")
            print(-H)
            .m <- -H
            .md = matrix(0, nETA, nETA)
            diag(.md) = abs(diag(.m)) * 1.1 + .001
            chol(.md)
        })
        ##H.neg.5 = chol(-H)
        log.det.H.neg.5 = sum(log(diag(H.neg.5)))

        RxODE_focei_eta_lik(tmp2.nm$eta, tmp2.nm)

        llik = -tmp2.nm$llik2 + log.det.OMGAinv.5             #note no -1/2 with OMGAinv.5
        llik.lapl =llik - log.det.H.neg.5               #note no 1/2
        llik.lapl = as.vector(llik.lapl);
        attr(llik.lapl, "fitted") <- tmp2.nm$f;
        attr(llik.lapl, "posthoc") <- tmp2.nm$eta;

        lik2 <- RxODE_focei_finalize_llik(tmp2.nm);

        expect_equal(lik2, llik.lapl);
    })

    m2ag <- rxSymPySetupPred(m2, pred, pk, err, grad=TRUE)

    tmp1 <- m2ag$outer %>% solve(et, theta=c(2, 1.6, 4.5,0.01), eta=c(0.01, -0.01))

    tmp2 <- m2ag %>% rxFoceiTheta(et, theta=c(2, 1.6, 4.5,0.01), eta=c(0.01, -0.01),dv=dv, inv.env=symenv)

    tmp2.nm <- m2ag %>% rxFoceiTheta(et, theta=c(2, 1.6, 4.5,0.01), eta=c(0.01, -0.01),dv=dv, inv.env=symenv, nonmem=TRUE)

    tmp5 <- m2ag %>% rxFoceiInner(et, theta=c(2, 1.6, 4.5,0.01), eta=c(0.01, -0.01),dv=dv, inv.env=symenv, invisible=TRUE)

    test_that("rxFoceiTheta makes sense", {
        expect_equal(tmp1$rx_pred_, tmp2$f); ## F
        err <- matrix(dv - tmp1$rx_pred_, ncol=1)
        expect_equal(err, tmp2$err) ## Err
        R <- matrix(tmp1$rx_r_, ncol=1)
        expect_equal(R, tmp2$R) ## R (Varinace)
        m <- as.matrix(tmp1[,c("_sens_rx_pred__ETA_1_", "_sens_rx_pred__ETA_2_")])
        dimnames(m) <- list(NULL, NULL)
        m2 <- tmp2$dErr
        dimnames(m2) <- list(NULL, NULL)
        expect_equal(m, m2) ## Check dErr
        m <- as.matrix(tmp1[, c("_sens_rx_r__ETA_1_", "_sens_rx_r__ETA_2_")]);
        dimnames(m) <- list(NULL, NULL);
        m2 <- tmp2$dR
        dimnames(m2) <- list(NULL, NULL)
        expect_equal(m, m2) ## Check dR
        c <- list(matrix(tmp1[["_sens_rx_r__ETA_1_"]] / tmp1[["rx_r_"]],ncol=1),matrix(tmp1[["_sens_rx_r__ETA_2_"]] / tmp1[["rx_r_"]],ncol=1))
        expect_equal(c, tmp2$c) ## c
        B <- matrix(2 / tmp1[["rx_r_"]],ncol=1)
        expect_equal(B, tmp2$B)

        a <- list(matrix(tmp1[["_sens_rx_pred__ETA_1_"]],ncol=1) - err/R*matrix(tmp1[["_sens_rx_r__ETA_1_"]]),
                  matrix(tmp1[["_sens_rx_pred__ETA_2_"]],ncol=1)- err/R*matrix(tmp1[["_sens_rx_r__ETA_2_"]]));

        expect_equal(a, tmp2$a);

        a <- list(matrix(tmp2.nm$dErr[, 1], ncol=1),
                  matrix(tmp2.nm$dErr[, 2], ncol=1))
        expect_equal(tmp2.nm$a, a)

        ## Does not include the matrix multilpaction part (done with RcppArmadillo)
        lp <- matrix(c(NA, NA), ncol=1)
        c <- matrix(tmp1[["_sens_rx_r__ETA_1_"]] / tmp1[["rx_r_"]],ncol=1)
        fp <- as.matrix(tmp1[,c("_sens_rx_pred__ETA_1_" )])
        lp[1, 1] <- .5*apply(err*fp*B + .5*err^2*B*c - c, 2, sum)
        c <- matrix(tmp1[["_sens_rx_r__ETA_2_"]] / tmp1[["rx_r_"]],ncol=1)
        fp <- as.matrix(tmp1[,c("_sens_rx_pred__ETA_2_" )])
        lp[2, 1] <- .5*apply(err*fp*B + .5*err^2*B*c - c, 2, sum)
        expect_equal(lp, tmp2$lp)

        llik <- -0.5 * sum(err ^ 2 / R + log(R));

        expect_equal(llik, tmp2$llik);

        ##
        expect_equal(tmp2$dErr.dTheta,
                     matrix(c(tmp1[["_sens_rx_pred__THETA_1_"]],tmp1[["_sens_rx_pred__THETA_2_"]],tmp1[["_sens_rx_pred__THETA_3_"]],tmp1[["_sens_rx_pred__THETA_4_"]],
                              rep(0, dim(tmp2$dErr.dTheta)[1] * 2)),ncol=6))

        expect_equal(tmp2$dR.dTheta,
                     matrix(c(tmp1[["_sens_rx_r__THETA_1_"]],tmp1[["_sens_rx_r__THETA_2_"]],tmp1[["_sens_rx_r__THETA_3_"]],tmp1[["_sens_rx_r__THETA_4_"]],
                              rep(0, dim(tmp2$dErr.dTheta)[1] * 2)),ncol=6))

        expect_equal(tmp2$dR2,
                     list(matrix(c(tmp1[["_sens_rx_r__BY_ETA_1__ETA_1_"]],tmp1[["_sens_rx_r__BY_ETA_1__ETA_2_"]]),ncol=2),
                          matrix(c(tmp1[["_sens_rx_r__BY_ETA_1__ETA_2_"]],tmp1[["_sens_rx_r__BY_ETA_2__ETA_2_"]]),ncol=2)
                          ))

        expect_equal(tmp2$dErr2,
                     list(matrix(c(tmp1[["_sens_rx_pred__BY_ETA_1__ETA_1_"]],tmp1[["_sens_rx_pred__BY_ETA_1__ETA_2_"]]),ncol=2),
                          matrix(c(tmp1[["_sens_rx_pred__BY_ETA_1__ETA_2_"]],tmp1[["_sens_rx_pred__BY_ETA_2__ETA_2_"]]),ncol=2)))

    })

}, silent=TRUE)
