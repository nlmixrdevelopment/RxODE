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

    ## Fixme?
    ## m2a2 <- rxSymPySetupPred(m2, pred)

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
        ## expect_equal(class(m2a2), "rxFocei")
        expect_equal(class(m2a), "rxFocei")
        expect_equal(class(m2b), "rxFocei")
        expect_equal(class(m2c), "rxFocei")
        expect_equal(class(m2d), "rxFocei")
        expect_equal(class(m2e), "rxFocei")
        expect_equal(class(m2f), "rxFocei")
        expect_true(length(rxInit(m2f$inner)) == 0)
    })

    context("Inner Problem Tests")

    ev <- eventTable() %>%
        add.sampling(c(95.99, 119.99, 143.99, 144.25, 144.5, 144.75,
                       145, 145.5, 146, 146.5, 147, 148, 150, 152, 156, 160, 164, 167.99,
                       191.99, 215.99, 216.25, 216.5, 216.75, 217, 217.5, 218, 218.5, 219,
                       220, 222, 224, 228, 232, 236, 240, 252, 264, 276, 288)) %>%
        add.dosing(dose=60000, start.time=72, nbr.doses=7, dosing.interval=24)

    dv <- c(263.6, 164.7, 287.3, 1248.7, 1211.5, 1017.7, 1690.1, 1029.8,
            890.7, 598.4, 1009.3, 1159.8, 742.2, 724.6, 728.2, 509.7, 243.1,
            259.9, 242.2, 281.4, 1500.1, 1281.4, 1200.2, 1378.8, 1373.2,
            582.9, 960.2, 720.3, 852.6, 950.3, 654.7, 402.5, 456, 346.5,
            268.2, 134.2, 42.6, 25.9, 14.6)

    m1 <- RxODE({
        C2 = centr/V;
        d/dt(centr) = - CL*C2
    })


    pred = function() C2


    mypar1 = function ()
    {
        CL = exp(THETA[1] + ETA[1])
        V = exp(THETA[2] + ETA[2])
    }

    m2a <- rxSymPySetupPred(m1, pred, mypar1, function(){err ~ prop(0.1)})

    omega <- matrix(c(0.1, 0, 0, 0.1), nrow=2)

    symo <- rxSymInvCreate(omega);

    symenv <- rxSymInv(symo, c(sqrt(0.1), sqrt(0.1)))

    THETA <- c(1.6, 4.5, sqrt(0.1));
    ETA <- c(-0.147736086922763, -0.294637022436797)

    tmp1 <- m2a$inner %>% rxSolve(ev, theta=THETA, eta=ETA)

    tmp2 <- m2a %>% rxFoceiEta(ev, theta=THETA, eta=ETA,dv=dv, inv.env=symenv)

    tmp2.nm <- m2a %>% rxFoceiEta(ev, theta=THETA, eta=ETA,dv=dv, inv.env=symenv, nonmem=TRUE)

    tmp3 <- m2a %>% rxFoceiLik(ev, theta=THETA, eta=ETA,dv=dv, inv.env=tmp2)

    tmp4 <- m2a %>% rxFoceiLp(ev, theta=THETA, eta=ETA,dv=dv, inv.env=symenv)

    tmp5 <- m2a %>%rxFoceiInner(ev, theta=THETA, eta=c(10, 10),dv=dv, inv.env=symenv, invisible=1)

    tmp5 <- m2a %>%rxFoceiInner(ev, theta=THETA, eta=c(10, 10),dv=dv, inv.env=symenv, invisible=1, pred.minus.dv=FALSE)

    tmp5a <- m2a %>%rxFoceiInner(ev, theta=THETA, eta=c(0, 0),dv=dv, inv.env=symenv, invisible=1, nonmem=TRUE)

    test_that("rxFoceiEta makes sense", {

        expect_equal(tmp1$rx_pred_, tmp2$f); ## F
        err <- matrix(tmp1$rx_pred_ - dv , ncol=1)
        expect_equal(err, tmp2$err) ## Err
        R <- matrix(tmp1$rx_r_, ncol=1) ## check
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
        expect_equal(c, tmp2$c) ## check
        B <- matrix(2 / tmp1[["rx_r_"]],ncol=1)
        expect_equal(B, tmp2$B) ## check

        ## This is different...  According to the Almquist paper 2015
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
        lp[1, 1] <- .5*apply(-err*fp*B + .5*err^2*B*c - c, 2, sum)
        c <- matrix(tmp1[["_sens_rx_r__ETA_2_"]] / tmp1[["rx_r_"]],ncol=1)
        fp <- as.matrix(tmp1[,c("_sens_rx_pred__ETA_2_" )])
        lp[2, 1] <- .5*apply(-err*fp*B + .5*err^2*B*c - c, 2, sum)
        expect_equal(lp, tmp2$lp)

        llik <- -0.5 * sum(err ^ 2 / R + log(R));

        expect_equal(llik, tmp2$llik);

        eta <- ETA
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
        dimnames(tmp2$H) <- list(NULL, NULL);
        dimnames(H) <- list(NULL, NULL);
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
        ##H.neg.5 = chol(-H)
        log.det.H.neg.5 = sum(log(diag(H.neg.5)))

        RxODE_focei_eta_lik(tmp2$eta, tmp2)
        log.det.OMGAinv.5 <- symenv$log.det.OMGAinv.5
        llik = -tmp2$llik2 + log.det.OMGAinv.5             #note no -1/2 with OMGAinv.5
        llik.lapl =llik - log.det.H.neg.5               #note no 1/2
        llik.lapl = as.vector(llik.lapl);
        attr(llik.lapl, "fitted") <- tmp2$f;
        attr(llik.lapl, "posthoc") <- tmp2$eta;
        attr(llik.lapl, "corrected") <- 0L

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
        dimnames(tmp2.nm$H) <- list(NULL, NULL);
        dimnames(H) <- list(NULL, NULL);
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
        attr(llik.lapl, "corrected") <- 0L

        lik2 <- RxODE_focei_finalize_llik(tmp2.nm);

        expect_equal(lik2, llik.lapl);

        expect_equal(as.vector(lik2), -209.467627968204);

        expect_equal(as.vector(tmp5a), -209.467627968204);

    })

    ## Make sure the eta f/grads are correct...

    m1 <- RxODE({
        C2 = centr/V;
        d/dt(centr) = - CL*C2
    })


    pred = function() C2


    mypar1 = function ()
    {
        CL = exp(THETA[1] + ETA[1])
        V = exp(THETA[2] + ETA[2])
    }

    context("Outer Problem Gradient Tests")

    m2ag <- rxSymPySetupPred(m1, pred, mypar1, function(){err ~ prop(0.1)}, grad=TRUE)

    tmp1 <- m2ag$outer %>% solve(ev, theta=THETA, eta=ETA)

    tmp2 <- m2ag %>% rxFoceiTheta(ev, theta=THETA, eta=ETA,dv=dv, inv.env=symenv)

    tmp2.nm <- m2ag %>% rxFoceiTheta(ev, theta=THETA, eta=ETA,dv=dv, inv.env=symenv, nonmem=TRUE)

    tmp5.g <- m2ag %>%rxFoceiInner(ev, theta=THETA, eta=ETA,dv=dv, inv.env=symenv, invisible=1)

    tmp5.g2 <- m2ag %>%rxFoceiInner(ev, theta=THETA, eta=ETA,dv=dv, inv.env=symenv, invisible=1,
                                    inits.vec=rep(0.5, 5))

    test_that("rxFoceiTheta makes sense", {

        expect_equal(tmp1$rx_pred_, tmp2$f); ## F
        err <- matrix(tmp1$rx_pred_ - dv, ncol=1)
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

        ## Does not include the matrix mult part (done with RcppArmadillo)
        lp <- matrix(c(NA, NA), ncol=1)
        c <- matrix(tmp1[["_sens_rx_r__ETA_1_"]] / tmp1[["rx_r_"]],ncol=1)
        expect_equal(c, tmp2$c[[1]])
        fp <- as.matrix(tmp1[,c("_sens_rx_pred__ETA_1_" )])
        lp[1, 1] <- .5*apply(-err*fp*B + .5*err^2*B*c - c, 2, sum)
        c <- matrix(tmp1[["_sens_rx_r__ETA_2_"]] / tmp1[["rx_r_"]],ncol=1)
        expect_equal(c, tmp2$c[[2]])
        fp <- as.matrix(tmp1[,c("_sens_rx_pred__ETA_2_" )])
        lp[2, 1] <- .5*apply(-err*fp*B + .5*err^2*B*c - c, 2, sum)

        expect_equal(lp, tmp2$lp)

        llik <- -0.5 * sum(err ^ 2 / R + log(R));

        expect_equal(llik, tmp2$llik);

        ##
        expect_equal(tmp2$dErr.dTheta,
                     matrix(c(tmp1[["_sens_rx_pred__THETA_1_"]],tmp1[["_sens_rx_pred__THETA_2_"]],tmp1[["_sens_rx_pred__THETA_3_"]],
                              rep(0, dim(tmp2$dErr.dTheta)[1] * 2)),ncol=5))

        expect_equal(tmp2$dR.dTheta,
                     matrix(c(tmp1[["_sens_rx_r__THETA_1_"]],tmp1[["_sens_rx_r__THETA_2_"]],tmp1[["_sens_rx_r__THETA_3_"]],
                              rep(0, dim(tmp2$dErr.dTheta)[1] * 2)),ncol=5))

        expect_equal(tmp2$dR2,
                     list(matrix(c(tmp1[["_sens_rx_r__BY_ETA_1__ETA_1_"]],tmp1[["_sens_rx_r__BY_ETA_1__ETA_2_"]]),ncol=2),
                          matrix(c(tmp1[["_sens_rx_r__BY_ETA_1__ETA_2_"]],tmp1[["_sens_rx_r__BY_ETA_2__ETA_2_"]]),ncol=2)
                          ))

        expect_equal(tmp2$dErr2,
                     list(matrix(c(tmp1[["_sens_rx_pred__BY_ETA_1__ETA_1_"]],tmp1[["_sens_rx_pred__BY_ETA_1__ETA_2_"]]),ncol=2),
                          matrix(c(tmp1[["_sens_rx_pred__BY_ETA_1__ETA_2_"]],tmp1[["_sens_rx_pred__BY_ETA_2__ETA_2_"]]),ncol=2)));

        ## Now Test dErr.dEta.dTheta
        expect_equal(tmp2$dErr.dEta.dTheta,
                     list(matrix(c(tmp1[["_sens_rx_pred__BY_ETA_1__THETA_1_"]],
                                   tmp1[["_sens_rx_pred__BY_ETA_2__THETA_1_"]]), ncol=2),
                          matrix(c(tmp1[["_sens_rx_pred__BY_ETA_1__THETA_2_"]],
                                   tmp1[["_sens_rx_pred__BY_ETA_2__THETA_2_"]]), ncol=2),
                          matrix(c(tmp1[["_sens_rx_pred__BY_ETA_1__THETA_3_"]],
                                   tmp1[["_sens_rx_pred__BY_ETA_2__THETA_3_"]]), ncol=2),
                          matrix(rep(0, 2 * length(tmp1[, 1])), ncol=2),
                          matrix(rep(0, 2 * length(tmp1[, 1])), ncol=2)));

        ## Now test dR.dEta.dTheta
        expect_equal(tmp2$dR.dEta.dTheta,
                     list(matrix(c(tmp1[["_sens_rx_r__BY_ETA_1__THETA_1_"]],
                                   tmp1[["_sens_rx_r__BY_ETA_2__THETA_1_"]]), ncol=2),
                          matrix(c(tmp1[["_sens_rx_r__BY_ETA_1__THETA_2_"]],
                                   tmp1[["_sens_rx_r__BY_ETA_2__THETA_2_"]]), ncol=2),
                          matrix(c(tmp1[["_sens_rx_r__BY_ETA_1__THETA_3_"]],
                                   tmp1[["_sens_rx_r__BY_ETA_2__THETA_3_"]]), ncol=2),
                          matrix(rep(0, 2 * length(tmp1[, 1])), ncol=2),
                          matrix(rep(0, 2 * length(tmp1[, 1])), ncol=2)))

        ## Now test H2.
        h2f <- function(k, l){
            ## Equation 13
            dErr.k <- tmp2$dErr[, k];
            dErr.l <- tmp2$dErr[, l];
            dR.k <- tmp2$dR[, k];
            dR.l <- tmp2$dR[, l];
            dErr.k.l <- tmp2$dErr2[[k]][, l];
            dR.k.l <- tmp2$dR2[[k]][, l];
            err <- tmp2$err;
            R <- tmp2$R;
            return(-0.5 * sum(2 * dErr.l * dErr.k / R - 2 * err * dR.l * dErr.k / (R * R) +
                              2 * err * dErr.k.l / R - err * err * dR.k.l / (R * R) +
                              2 * err * err * dR.k * dR.l / (R * R * R) -
                              2 * err * dR.k * dErr.l / (R * R) -
                              dR.k * dR.l / (R * R) + dR.k.l / R) - symenv$omegaInv[k, l]);
        }

        h2 <- matrix(c(h2f(1, 1), h2f(2, 1), h2f(1, 2), h2f(2, 2)), ncol=2)

        expect_equal(h2, tmp2$H2);

        ## Now test l.dEta.dTheta

        lEH <- function(k, m){
            ## Equation 47
            dErr.m <- tmp2$dErr.dTheta[, m];
            dErr.k <- tmp2$dErr[, k];
            dErr.k.m <- tmp2$dErr.dEta.dTheta[[m]][, k];
            ##
            dR.m <- tmp2$dR.dTheta[, m];
            dR.k <- tmp2$dR[, k];
            dR.k.m <- tmp2$dR.dEta.dTheta[[m]][, k];
            ##
            err <- tmp2$err;
            R <- tmp2$R;
            nomega <- length(tmp2$dOmega)
            ntheta <- tmp2$ntheta;
            if (m > ntheta){
                ome <- tmp2$omega.47[k, m - ntheta];
            } else {
                ome <- 0;
            }
            return(-0.5 * sum(2 * dErr.m * dErr.k / R -
                              2 * err * dR.m * dErr.k / (R * R) +
                              2 * err * dErr.k.m / R -
                              err * err * dR.k.m / (R * R) +
                              2 * err * err * dR.k * dR.m / (R * R * R) -
                              2 * err * dR.k * dErr.m / (R * R) +
                              dR.m * dR.k / (R * R) +
                              dR.k.m / R) - ome);
        }

        df <- expand.grid(k=c(1, 2), theta=1:5);

        tmp2a <- matrix(as.vector(apply(df, 1, function(x) {return(lEH(x[1], x[2]))})), nrow=2)

        expect_equal(tmp2a, tmp2$l.dEta.dTheta)

        ## Now test  #46
        expect_equal(tmp2$dEta.dTheta, -solve(tmp2$H2) %*% tmp2$l.dEta.dTheta);

        ## Now test dErr.dTheta. (Equation #33)
        f <- function(m){
            tmp2$dErr.dTheta[, m] + tmp2$dErr %*% tmp2$dEta.dTheta[, m]
        }

        expect_equal(tmp2$dErr.dTheta., cbind(f(1), f(2), f(3), f(4), f(5)))

        f <- function(m){
            tmp2$dR.dTheta[, m] + tmp2$dR %*% tmp2$dEta.dTheta[, m]
        }

        expect_equal(tmp2$dR.dTheta., cbind(f(1), f(2), f(3), f(4), f(5)));

        ## Now #37
        f <- function(k, m){
            dErr.k.m <- tmp2$dErr.dEta.dTheta[[m]][, k];
            err1 <- tmp2$dErr2[[k]];
            dNdH <- tmp2$dEta.dTheta[, m];
            return(dErr.k.m + err1 %*% dNdH);
        }

        expect_equal(tmp2$dErr.dEta.dTheta., list(matrix(c(f(1, 1), f(1, 2), f(1, 3), f(1, 4), f(1, 5)), ncol=5),
                                                  matrix(c(f(2, 1), f(2, 2), f(2, 3), f(2, 4), f(2, 5)), ncol=5)))

        ## Now #37 for dR
        f <- function(k, m){
            dR.k.m <- tmp2$dR.dEta.dTheta[[m]][, k];
            err1 <- tmp2$dR2[[k]];
            dNdH <- tmp2$dEta.dTheta[, m];
            return(dR.k.m + err1 %*% dNdH);
        }

        expect_equal(tmp2$dR.dEta.dTheta., list(matrix(c(f(1, 1), f(1, 2), f(1, 3), f(1, 4), f(1, 5)), ncol=5),
                                                matrix(c(f(2, 1), f(2, 2), f(2, 3), f(2, 4), f(2, 5)), ncol=5)))

        ## Now dc/dTheta
        f <- function(k, m){
            dR.m <- tmp2$dR.dTheta.[, m];
            dR.k <- tmp2$dR[, k];
            dR.k.m <- tmp2$dR.dEta.dTheta.[[k]][, m];
            R <- tmp2$R;
            return(-dR.m * dR.k / (R * R) + dR.k.m / R);
        }

        expect_equal(tmp2$dc.dTheta, list(matrix(c(f(1, 1), f(1, 2), f(1, 3), f(1, 4), f(1, 5)), ncol=5),
                                          matrix(c(f(2, 1), f(2, 2), f(2, 3), f(2, 4), f(2, 5)), ncol=5)))

        ## Now dB/dTheta
        f <- function(m){
            dR.m <- tmp2$dR.dTheta.[, m];
            R <- tmp2$R;
            return(-2 * dR.m / (R * R));
        }

        expect_equal(tmp2$dB.dTheta, matrix(c(f(1), f(2), f(3), f(4), f(5)), ncol=5))

        ## Now da/dTheta
        f <- function(k, m){
            dErr.m.k <- tmp2$dErr.dEta.dTheta.[[k]][, m];
            dErr.m <- tmp2$dErr.dTheta.[, m];
            dR.k <- tmp2$dR[, k];
            dR.m <- tmp2$dR.dTheta.[, m];
            dR.m.k <- tmp2$dR.dEta.dTheta.[[k]][, m];
            err <- tmp2$err;
            R <- tmp2$R;
            return(dErr.m.k - dErr.m * dR.k / R +
                   err * dR.m * dR.k / (R * R) -
                   err * dR.m.k / R);
        }

        expect_equal(tmp2$da.dTheta,
                     list(matrix(c(f(1, 1), f(1, 2), f(1, 3), f(1, 4), f(1, 5)), ncol=5),
                          matrix(c(f(2, 1), f(2, 2), f(2, 3), f(2, 4), f(2, 5)), ncol=5)));

        f <- function(k, m){
            dErr.m.k <- tmp2$dErr.dEta.dTheta.[[k]][, m];
            return(dErr.m.k);
        }

        expect_equal(tmp2.nm$da.dTheta,
                     list(matrix(c(f(1, 1), f(1, 2), f(1, 3), f(1, 4), f(1, 5)), ncol=5),
                          matrix(c(f(2, 1), f(2, 2), f(2, 3), f(2, 4), f(2, 5)), ncol=5)));

        ## Test dH/dTheta

        f <- function(m, k, l){
            da.l.m <- tmp2$da.dTheta[[l]][, m];
            B <- tmp2$B
            a.k <- tmp2$a[[k]];
            a.l <- tmp2$a[[l]];
            B.m <- tmp2$dB.dTheta[, m];
            da.k.m <- tmp2$da.dTheta[[k]][, m];
            dc.l.m <- tmp2$dc.dTheta[[l]][, m];
            dc.k.m <- tmp2$dc.dTheta[[k]][, m];
            c.k <- tmp2$c[[k]];
            c.l <- tmp2$c[[l]];
            if (m > 3){
                ome <- symenv$dOmegaInv[[m - 3]][k, l];
            } else {
                ome <- 0;
            }
            return(-0.5 * sum(da.l.m * B * a.k +
                              a.l * B.m * a.k +
                              a.l * B * da.k.m -
                              dc.l.m * c.k -
                              c.l * dc.k.m) - ome);
        }

        df <- expand.grid(m=1:5, k=1:2, l=1:2);
        df <- df[order(df$m, df$k, df$l), ];

        v <- apply(df, 1, function(x) {return(f(x[1], x[2], x[3]))})

        expect_equal(tmp2$dH.dTheta,
                     list(matrix(v[1:4], 2),
                          matrix(v[5:8], 2),
                          matrix(v[9:12], 2),
                          matrix(v[13:16], 2),
                          matrix(v[17:20], 2)))

        f <- function(m, k, l){
            da.l.m <- tmp2.nm$da.dTheta[[l]][, m];
            B <- tmp2.nm$B
            a.k <- tmp2.nm$a[[k]];
            a.l <- tmp2.nm$a[[l]];
            B.m <- tmp2.nm$dB.dTheta[, m];
            da.k.m <- tmp2.nm$da.dTheta[[k]][, m];
            dc.l.m <- tmp2.nm$dc.dTheta[[l]][, m];
            dc.k.m <- tmp2.nm$dc.dTheta[[k]][, m];
            c.k <- tmp2.nm$c[[k]];
            c.l <- tmp2.nm$c[[l]];
            if (m > 3){
                ome <- symenv$dOmegaInv[[m - 3]][k, l];
            } else {
                ome <- 0;
            }
            return(-0.5 * sum(da.l.m * B * a.k +
                              a.l * B.m * a.k +
                              a.l * B * da.k.m +
                              dc.l.m * c.k +
                              c.l * dc.k.m) - ome);
        }

        v <- apply(df, 1, function(x) {return(f(x[1], x[2], x[3]))})

        expect_equal(tmp2.nm$dH.dTheta,
                     list(matrix(v[1:4], 2),
                          matrix(v[5:8], 2),
                          matrix(v[9:12], 2),
                          matrix(v[13:16], 2),
                          matrix(v[17:20], 2)));

        ## Test Inverses
        expect_equal(solve(tmp2$H), tmp2$Hinv);
        expect_equal(solve(tmp2.nm$H), tmp2.nm$Hinv);

        f <- function(m){
            dErr.m <- tmp2$dErr.dTheta[, m];
            dR.m <- tmp2$dR.dTheta[, m];
            err <- tmp2$err;
            R <- tmp2$R;
            eta <- tmp2$eta.mat
            if (m > 3){
                ome <- tmp2$omega.28[m - 3]
            } else {
                ome <- 0;
            }
            Hinv <- tmp2$Hinv;
            dH <- tmp2$dH.dTheta[[m]];
            return(-0.5 * sum(2 * err * dErr.m / R -
                              err * err * dR.m / (R * R) +
                              dR.m / R)  + ome - 0.5 * sum(diag(Hinv %*% dH)));
        }

        expect_equal(tmp2$l.dTheta, c(f(1), f(2), f(3), f(4), f(5)))

        f <- function(m){
            dErr.m <- tmp2.nm$dErr.dTheta[, m];
            dR.m <- tmp2.nm$dR.dTheta[, m];
            err <- tmp2.nm$err;
            R <- tmp2.nm$R;
            eta <- tmp2.nm$eta.mat
            if (m > 3){
                ome <- tmp2.nm$omega.28[m - 3]
            } else {
                ome <- 0;
            }
            Hinv <- tmp2.nm$Hinv;
            dH <- tmp2.nm$dH.dTheta[[m]];
            return(-0.5 * sum(2 * err * dErr.m / R -
                              err * err * dR.m / (R * R) +
                              dR.m / R)  + ome - 0.5 * sum(diag(Hinv %*% dH)));
        }

        expect_equal(tmp2.nm$l.dTheta, c(f(1), f(2), f(3), f(4), f(5)))

        ## Now test omega.28
        f <- function(m){
            eta <- tmp2$eta.mat;
            omegaInv <- symenv$omegaInv;
            dOmega <- symenv$dOmega[[m]];
            return(0.5 * t(eta) %*% omegaInv %*% dOmega %*% omegaInv %*% eta - 0.5 * sum(diag(omegaInv %*% dOmega)))
        }

        expect_equal(tmp2$omega.28, c(f(1), f(2)))

        ## Now test omega.47

        f <- function(m, k){
            eta <- tmp2$eta.mat;
            eta1 <- rep(0, length(as.vector(eta)));
            eta1[k] <- 1;
            eta1 <- matrix(eta1, ncol=1);
            omegaInv <- symenv$omegaInv;
            dOmega <- symenv$dOmega[[m]];
            return(t(eta) %*% omegaInv %*% dOmega %*% omegaInv %*% eta1);
        }

        df <- expand.grid(m=1:2, k=1:2);
        df <- df[order(df$m), ];

        v <- matrix(apply(df, 1, function(x) {return(f(x[1], x[2]))}), nrow=2);

        expect_equal(v, tmp2$omega.47);

        ## Test scaling
        expect_equal(attr(tmp5.g,"grad") * 2, attr(tmp5.g2,"grad"));
        expect_equal(attr(tmp5.g,"dEta.dTheta"), attr(tmp5.g2,"dEta.dTheta"));

    })

    context("Test Initial conditions -> sensitivity initial conditions")

    fini <- RxODE({
        C2 = centr/V2;
        C3 = peri/V3;
        d/dt(depot) =-KA*depot;
        d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
        d/dt(peri)  =                    Q*C2 - Q*C3;
        d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
        eff(0) = theta[5]  + eta[4] + sqrt(theta[5] * eta[4]); ## doesn't work here.
    })

    pk <- function(){
        KA = exp(THETA[1])
        CL = exp(THETA[2] + ETA[1])
        V2 = exp(THETA[3] + ETA[2])
        V3 = exp(THETA[4] + ETA[3])
    }

    pred <- function(){
        if (cmt == 1){
            return(centr);
        } else {
            return(eff);
        }
    }

    err <- function(f){
        return(f ^ 2* theta[6] ^ 2); ## Theta 4 is residual sd for proportional error.
    }

    finip <- rxSymPySetupPred(fini, pred, pk, err, grad=TRUE);

    inner <- strsplit(rxNorm(finip$inner), "\n")[[1]];
    inner <- gsub(rex::rex(start, anything, "=", capture(anything), ";", end), "\\1",
                  inner[regexpr(rex::rex("rx__sens_eff_BY_ETA_4___(0)="), inner) != -1]);

    test_that("Inner sensitivites initial conditions are calculated appropriately.", {
        expect_equal(length(inner), 1)
        expect_true(regexpr(rex::rex("THETA[5]"), inner) != -1);
        expect_true(regexpr(rex::rex("ETA[4]"), inner) != -1);
    })

    out <- strsplit(rxNorm(finip$out), "\n")[[1]];
    out <- gsub(rex::rex(start, anything, "=", capture(anything), ";", end), "\\1",
                  out[regexpr(rex::rex(or("rx__sens_eff_BY_ETA_4___(0)=",
                                       "rx__sens_eff_BY_ETA_4__BY_ETA_4___(0)",
                                       "rx__sens_eff_BY_THETA_5___(0)")), out) != -1]);

    test_that("Out sensitivites initial conditions are calculated appropriately.", {
        expect_equal(length(out), 3)
        expect_true(all(regexpr(rex::rex("THETA[5]"), out) != -1));
        expect_true(all(regexpr(rex::rex("ETA[4]"), out) != -1));
        expect_true(all(!duplicated(out)));
    })


    fini <- RxODE({
        C2 = centr/V2;
        C3 = peri/V3;
        d/dt(depot) =-KA*depot;
        d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
        d/dt(peri)  =                    Q*C2 - Q*C3;
        d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
    })

    pk <- function(){
        KA = exp(THETA[1])
        CL = exp(THETA[2] + ETA[1])
        V2 = exp(THETA[3] + ETA[2])
        V3 = exp(THETA[4] + ETA[3])
        eff(0) = theta[5]  + eta[4] + sqrt(theta[5] * eta[4]); ## doesn't work here.
    }

    pred <- function() eff

    err <- function(f){
        return(f ^ 2* theta[6] ^ 2); ## Theta 4 is residual sd for proportional error.
    }

    finip <- rxSymPySetupPred(fini, pred, pk, err, grad=TRUE);

    inner <- strsplit(rxNorm(finip$inner), "\n")[[1]];
    inner <- gsub(rex::rex(start, anything, "=", capture(anything), ";", end), "\\1",
                  inner[regexpr(rex::rex("rx__sens_eff_BY_ETA_4___(0)="), inner) != -1]);

    test_that("Inner sensitivites initial conditions are calculated appropriately #2", {
        expect_equal(length(inner), 1)
        expect_true(regexpr(rex::rex("THETA[5]"), inner) != -1);
        expect_true(regexpr(rex::rex("ETA[4]"), inner) != -1);
    })

    out <- strsplit(rxNorm(finip$out), "\n")[[1]];
    out <- gsub(rex::rex(start, anything, "=", capture(anything), ";", end), "\\1",
                out[regexpr(rex::rex(or("rx__sens_eff_BY_ETA_4___(0)=",
                                        "rx__sens_eff_BY_ETA_4__BY_ETA_4___(0)",
                                        "rx__sens_eff_BY_THETA_5___(0)")), out) != -1]);

    test_that("Out sensitivites initial conditions are calculated appropriately #2", {
        expect_equal(length(out), 3)
        expect_true(all(regexpr(rex::rex("THETA[5]"), out) != -1));
        expect_true(all(regexpr(rex::rex("ETA[4]"), out) != -1));
        expect_true(all(!duplicated(out)));
    })

    ## Now use initCondition statement.

    pk <- function(){
        KA = exp(THETA[1])
        CL = exp(THETA[2] + ETA[1])
        V2 = exp(THETA[3] + ETA[2])
        V3 = exp(THETA[4] + ETA[3])
        initCondition = c(0, 0, 0, theta[5]  + eta[4] + sqrt(theta[5] * eta[4]));
    }

    finip <- rxSymPySetupPred(fini, pred, pk, err, grad=TRUE);

    inner <- strsplit(rxNorm(finip$inner), "\n")[[1]];
    inner <- gsub(rex::rex(start, anything, "=", capture(anything), ";", end), "\\1",
                  inner[regexpr(rex::rex("rx__sens_eff_BY_ETA_4___(0)="), inner) != -1]);

    test_that("Inner sensitivites initial conditions are calculated appropriately #3", {
        expect_equal(length(inner), 1)
        expect_true(regexpr(rex::rex("THETA[5]"), inner) != -1);
        expect_true(regexpr(rex::rex("ETA[4]"), inner) != -1);
    })

    out <- strsplit(rxNorm(finip$out), "\n")[[1]];
    out <- gsub(rex::rex(start, anything, "=", capture(anything), ";", end), "\\1",
                out[regexpr(rex::rex(or("rx__sens_eff_BY_ETA_4___(0)=",
                                        "rx__sens_eff_BY_ETA_4__BY_ETA_4___(0)",
                                        "rx__sens_eff_BY_THETA_5___(0)")), out) != -1]);

    test_that("Out sensitivites initial conditions are calculated appropriately #3", {
        expect_equal(length(out), 3)
        expect_true(all(regexpr(rex::rex("THETA[5]"), out) != -1));
        expect_true(all(regexpr(rex::rex("ETA[4]"), out) != -1));
        expect_true(all(!duplicated(out)));
    })

    mod <- RxODE({
        C2 = centr/V2;
        C3 = peri/V3;
        d/dt(centr) =  - CL*C2 - Q*C2 + Q*C3;
        d/dt(peri)  = Q*C2 - Q*C3;
    })



    par <- function ()
    {
        CL <- THETA[1] * exp(ETA[1])
        V2 <- THETA[2]
        Q <- THETA[3]  * exp(ETA[2]);
        V3 <- THETA[4]
    }

    pred <- function() C2


    focei.mod1 <- rxSymPySetupPred(mod, pred, par, err=function(){return(prop(0.1) + add(0.1))});

    focei.mod2 <- rxSymPySetupPred(mod, pred, par, err=function(){return(prop(0.1) + add(0.1))}, grad=TRUE);

    test_that("Can use a LHS quantity to caluclate PRED.", {
        expect_equal(class(focei.mod1), "rxFocei");
        expect_equal(class(focei.mod2), "rxFocei");
    })

    context("Test actual gradients")
    ## Test the gradient for a single subject
    ev <- eventTable() %>%
        add.sampling(c(0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 12, 16, 20, 24,
                       36, 48, 60, 71.99, 95.99, 119.99, 143.99, 144.25, 144.5, 144.75,
                       145, 145.5, 146, 146.5, 147, 148, 150, 152, 156, 160, 164, 167.99,
                       191.99, 215.99, 216.25, 216.5, 216.75, 217, 217.5, 218, 218.5, 219,
                       220, 222, 224, 228, 232, 236, 240, 252, 264, 276, 288)) %>%
        add.dosing(dose=1000, start.time=0, nbr.doses=1) %>%
        add.dosing(dose=10000, start.time=72, nbr.doses=7, dosing.interval=24)

    DV <- c(0, 137.6, 142.1, 113, 134.8, 143.4, 106.9, 64.2, 118.5, 73.9, 90.8, 83.1,
            61, 75.4, 46.3, 54.5, 29.8, 29.3, 12.5, 13.2, 0, 64.9, 0, 81.5, 0, 108.3, 0, 243.3,
            185.5, 190.5, 169.9, 180.7, 205, 163.2, 169.5, 224.2, 126.3, 129.1, 107.7, 146.2, 96,
            113.2, 0, 79.7, 0, 85.3, 0, 226.9, 172, 164.2, 148.7, 200.9, 151, 160.6, 142.4, 154.8,
            126.9, 134.8, 129.1, 110.4, 121.6, 78.9, 37, 37.2, 28.3, 27.7)

    THETA <- c(1.59170802337068, 4.45624870814774, 0.330113291259874, 0.151067662733064, 0.132428482079086)

    mypar1 = function ()
    {
        CL = exp(THETA[1] + ETA[1])
        V = exp(THETA[2] + ETA[2])
    }

    m1 <- RxODE({
        C2 = centr/V;
        d/dt(centr) = - CL*C2
    })

    pred = function() C2

    err = function(){err ~ prop(.1)}

    m1g <- rxSymPySetupPred(m1, pred, mypar1, err, grad=TRUE, logify=FALSE)

    ETA <- c(0, 0);

    symo <- rxSymInvCreate(structure(c(0.1, 0, 0, 0.11), .Dim = c(2L, 2L)),
                           diag.xform="sqrt")

    symenv <- rxSymInv(symo, THETA[4:5])

    ret <- m1g %>% rxFoceiInner(ev, theta=THETA[-(4:5)], eta=ETA,
                                dv=DV, inv.env=symenv, NONMEM=1, invisible=1,
                                rtol.outer=1e-12, atol.outer= 1e-13)

    m1g <- rxSymPySetupPred(m1, pred, mypar1, err, grad=TRUE, logify=TRUE)

    ## ret2 <- m1g %>% rxFoceiTheta(et, theta=THETA, eta=ETA,dv=DV, inv.env=symenv)

    ret2 <- m1g %>% rxFoceiInner(ev, theta=THETA[-(4:5)], eta=ETA,
                                dv=DV, inv.env=symenv, NONMEM=1, invisible=1,
                                rtol.outer=1e-12, atol.outer= 1e-13)

    m1g2 <- rxSymPySetupPred(m1, pred, mypar1, err, grad=TRUE, logify=TRUE, pred.minus.dv=FALSE)

    ## ret2 <- m1g %>% rxFoceiTheta(et, theta=THETA, eta=ETA,dv=DV, inv.env=symenv)

    ret2a <- m1g2 %>% rxFoceiInner(ev, theta=THETA[-(4:5)], eta=ETA,
                                 dv=DV, inv.env=symenv, NONMEM=1, invisible=1,
                                 rtol.outer=1e-12, atol.outer= 1e-13, pred.minus.dv=FALSE)

    library(numDeriv)

    f <- function(th=THETA[-(4:5)]){
        return(suppressWarnings(-2 * {as.vector(m1g %>% rxFoceiInner(ev, theta=th, eta=as.vector(attr(ret2,"posthoc")),
                                                                     dv=DV, inv.env=symenv, NONMEM=1, invisible=1,
                                                                     estimate=FALSE))}))
    }

    gr.richard <- grad(f, THETA[-(4:5)])

    gr.simple <- grad(f, THETA[-(4:5)], method="simple")

    ## While it is close, it isn't exactly the same.
    gr.calc <- -2 * attr(ret2, "grad")


    f <- function(th=THETA[-(4:5)]){
        return(suppressWarnings({as.vector(m1g %>% rxFoceiInner(ev, theta=th, eta=as.vector(attr(ret2a,"posthoc")),
                                                                dv=DV, inv.env=symenv, NONMEM=1, invisible=1))}))
    }

    gr.richard <- grad(f, THETA[-(4:5)])

    gr.simple <- grad(f, THETA[-(4:5)], method="simple")

    ## While it is close, it isn't exactly the same.
    gr.calc <- attr(ret2a, "grad")

    ## m1g$outer <- RxODE(rxLogifyModel(m1g$outer))

    ## ret2 <- m1g %>% rxFoceiInner(ev, theta=THETA[-(4:5)], eta=ETA,
    ##                             dv=DV, inv.env=symenv, NONMEM=1, invisible=1,
    ##                             rtol.outer=1e-12, atol.outer= 1e-13)
    ## The numerical values may not be right from NumDeriv either
    ## gr2.calc <- attr(ret2, "grad")

    context("Michelis Menton test")
    ## Michelis Menton test
    mm <- RxODE({
        C2 = centr / V;
        d/dt(centr)  = -(VM*C2)/(KM+C2);
    })

    mypar3 <- function(lVM, lKM, lV )
    {
        VM = exp(theta[1] + eta[1])
        KM = exp(theta[2] + eta[2])
        V  = exp(theta[3] + eta[3])
    }

    pred <- function() C2

    ## focei.mm.mod2 <- rxSymPySetupPred(mm, pred, par, err=function(){err ~ prop(0.1)}, grad=TRUE, run.internal=TRUE);

    ## Now test solved focei capability.
    sol.1c.ka <- RxODE({
        ## double solvedC(double t, int parameterization, int cmt, unsigned int col, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double p8);
        sf=solvedC(t, 1, 1, 1, V, CL, 0, 0, 0, 0, KA, 0);
    })

    prd <- function(){
        return(sf);
    }

    err <- function(f){
        return(f ^ 2* theta[4] ^ 2); ## Theta 4 is residual sd for proportional error.
    }

    mypar1 = function ()
    {
        CL = exp(THETA[1] + ETA[1])
        V = exp(THETA[2] + ETA[2])
    }

    rxSymPySetupPred(sol.1c.ka, prd, mypar1, err, run.internal=TRUE)

}, silent=TRUE)
