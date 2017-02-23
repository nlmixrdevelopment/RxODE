context("Test RxODE THETA/ETA support")
library(digest)
rxPermissive({

    rigid <- RxODE({
        y1(0)    = 1
        y2(0)    = 0
        y3(0)    = 0.9
        a1       = theta[1]
        a2       = theta[2]
        a3       = theta[3] + eta[1]
        d/dt(y1) = a1*y2*y3
        d/dt(y2) = a2*y1*y3
        d/dt(y3) = a3*y1*y2
    })

    et <- eventTable();
    et$add.sampling(seq(0,20,by=0.01))

    out <- solve(rigid,et, theta=c(-2, 1.25, -0.5), eta=c(0))

    test_that("Test rigid body example",{
        expect_equal(digest(round(as.data.frame(out),3)),
                     "14ac439006ceb32455ba00c7817f9163")
    })

    out <- solve(rigid,et, theta=c(-2, 1.25, -0.5), eta=c(1))

    test_that("Test rigid body example",{
        expect_equal(digest(round(as.data.frame(out),3)),
                     "c7cffaa650a47e2b28e4cba99c603dde")
    })

    m1 <- RxODE({
        KA = exp(THETA[1])
        CL = exp(THETA[2] + ETA[1])
        V = exp(THETA[3] + ETA[2])
        d/dt(depot) = -KA*depot;
        d/dt(centr) = KA*depot - CL / V*centr;
    });


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
        return(f * theta[4]); ## Theta 4 is residual sd for proportional error.
    }

    m2a <- rxSymPySetupPred(m2, pred, pk, err)

    dv <- c(1126.1, 869.9, 883.6, 1244, 995.2, 946.4, 589.2, 754.4, 1060.8, 433, 609.5, 630.1, 342.1, 322.9, 196.3, 89.8, 42.8, 16.6, 8.4)

    omega <- matrix(c(0.1, 0, 0, 0.1), nrow=2)

    omegaInv <- solve(omega)

    ## test_that("theta/eta solve works", {
    ##     expect_equal(suppressWarnings(digest(m2a %>% solve(et, theta=c(2, 1.6, 4.5, 0.1), eta=c(0.01, -0.01)) %>% as.data.frame %>% round(3))),
    ##                  "13baa6cabb14367e1328149d3c14ebf7");
    ## })

    ## m2a <- rxSymPySetupPred(m2, pred, pk, err)

    tmp1 <- m2a %>% solve(et, theta=c(2, 1.6, 4.5,0.01), eta=c(0.01, -0.01))
    tmp2 <- m2a %>% rxFoceiEta(et, theta=c(2, 1.6, 4.5,0.01), eta=c(0.01, -0.01),dv=dv, omegaInv=omegaInv)

    tmp3 <- m2a %>% rxFoceiLik(et, theta=c(2, 1.6, 4.5,0.01), eta=c(0.01, -0.01),dv=dv, omegaInv=omegaInv)
    tmp4 <- m2a %>% rxFoceiLp(et, theta=c(2, 1.6, 4.5,0.01), eta=c(0.01, -0.01),dv=dv, omegaInv=omegaInv)

    tmp5 <- m2a %>%rxFoceiInner(et, theta=c(2, 1.6, 4.5,0.01), eta=c(10, 10),dv=dv, omegaInv=omegaInv, invisible=1)

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
        llik <- -0.5 * sum(err ^ 2 / R + log(R)) - 0.5 * t(matrix(eta,ncol=1)) %*% omegaInv %*% matrix(eta,ncol=1);
        llik <- -llik;

        expect_equal(llik, tmp3);

        lp2 <- lp - omegaInv %*% matrix(eta, ncol=1)

        expect_equal(-lp2, tmp4);

    })


}, silent=TRUE)
