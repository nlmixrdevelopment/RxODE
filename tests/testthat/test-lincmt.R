context("Test 1 compartment parameterizations (Cl, V)")
e <- new.env(parent = emptyenv());
test_that("1 cmt CL/V", {
    for (i  in 1:25){
        e$CL <- runif(1, 0.001, 1000)
        e$V <- runif(1, 0.001, 1000)
        k <- e$CL / e$V;
        V <- e$V;
        A <- 1.0 / V;
        alpha <- k;
        getMacroConstants(e);
        expect_equal(e$A, A)
        expect_equal(e$alpha, alpha);
        g <- matrix(c(alpha, A),nrow=1)
        expect_equal(e$g, g);
        expect_equal(e$oral, 0L);
        expect_equal(e$parameterization, 1L);
        expect_equal(e$ncmt, 1L);
        getLinDerivs(e);
        expect_true(all(ls(e)[regexpr("^d",ls(e))!= -1] %in% c("dCL", "dV")))
    }
})
context("\tOral");

e <- new.env(parent = emptyenv());
test_that("oral (Ka) 1 cmt CL/V", {
    for (i  in 1:25){
        e$CL <- runif(1, 0.001, 1000)
        e$V <- runif(1, 0.001, 1000)
        e$KA <- runif(1, 0.001, 1000)
        k <- e$CL / e$V;
        V <- e$V;
        A <- 1.0 / V;
        alpha <- k;
        ka <- e$KA;
        getMacroConstants(e);
        expect_equal(e$A, A)
        expect_equal(e$alpha, alpha);
        g <- matrix(c(alpha, ka / (ka - alpha) *A),nrow=1)
        expect_equal(e$g, g);
        expect_equal(e$oral, 1L);
        expect_equal(e$parameterization, 1L);
        expect_equal(e$ncmt, 1L);
        getLinDerivs(e);
        expect_true(all(ls(e)[regexpr("^d",ls(e))!= -1] %in% c("dCL", "dV", "dKA")))
    }
})

context("Test 1 compartment parameterizations (K, V)")
test_that("1 cmt K/V", {
    e <- new.env(parent = emptyenv());
    for (i  in 1:25){
        e$K <- runif(1, 0.001, 1000)
        e$V <- runif(1, 0.001, 1000)
        k <- e$K
        V <- e$V;
        A <- 1.0 / V;
        alpha <- k;
        getMacroConstants(e);
        expect_equal(e$A, A)
        expect_equal(e$alpha, alpha);
        g <- matrix(c(alpha, A),nrow=1)
        expect_equal(e$g, g);
        expect_equal(e$oral, 0L);
        expect_equal(e$parameterization, 2L);
        expect_equal(e$ncmt, 1L);
        getLinDerivs(e);
        expect_true(all(ls(e)[regexpr("^d",ls(e))!= -1] %in% c("dK", "dV")))
    }
})
context("\tOral");
e <- new.env(parent = emptyenv());
test_that("oral (Ka) 1 cmt K/V", {
    for (i  in 1:25){
        e$K <- runif(1, 0.001, 1000)
        e$V <- runif(1, 0.001, 1000)
        e$KA <- runif(1, 0.001, 1000)
        k <- e$K;
        V <- e$V;
        A <- 1.0 / V;
        alpha <- k;
        ka <- e$KA;
        getMacroConstants(e);
        expect_equal(e$A, A)
        expect_equal(e$alpha, alpha);
        g <- matrix(c(alpha, ka / (ka - alpha) *A),nrow=1)
        expect_equal(e$g, g);
        expect_equal(e$oral, 1L);
        expect_equal(e$parameterization, 2L);
        expect_equal(e$ncmt, 1L);
        getLinDerivs(e);
        expect_true(all(ls(e)[regexpr("^d",ls(e))!= -1] %in% c("dK", "dV", "dKA")))
    }
})


context("Test 2 compartment parameterizations (Cl, V, Q, V2)")

test_that("2 cmt CL V Q, V2", {
    e <- new.env(parent = emptyenv());
    for (i  in 1:25){
        e$CL <- runif(1, 0.001, 1000)
        e$V <- runif(1, 0.001, 1000)
        e$Q <- runif(1, 0.001, 1000)
        e$V2 <- runif(1, 0.001, 1000)
        k <- e$CL / e$V;
        V <- e$V;
        k12 <- e$Q / e$V;
        k21 <- e$Q / e$V2;
        beta <- 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k) * (k12 + k21 + k) - 4.0 * k21 * k));
        alpha <- k21 * k / beta;
        A <- (alpha - k21) / (alpha - beta) / V;
        B <- (beta - k21) / (beta - alpha) / V;
        getMacroConstants(e);
        expect_equal(e$A, A)
        expect_equal(e$alpha, alpha);
        expect_equal(e$B, B)
        expect_equal(e$beta, beta);
        g <- matrix(c(alpha, A,
                      beta, B),nrow=2, byrow=TRUE);
        expect_equal(e$g, g);
        expect_equal(e$oral, 0L);
        expect_equal(e$parameterization, 1L);
        expect_equal(e$ncmt, 2L);
        getLinDerivs(e);
        expect_true(all(ls(e)[regexpr("^d",ls(e))!= -1] %in% c("dCL", "dV", "dQ", "dV2")))
    }
})

context("\tOral");

test_that("oral 2 cmt CL V Q, V2", {
    e <- new.env(parent = emptyenv());
    for (i  in 1:25){
        e$CL <- runif(1, 0.001, 1000)
        e$V <- runif(1, 0.001, 1000);
        e$Q <- runif(1, 0.001, 1000);
        e$V2 <- runif(1, 0.001, 1000);
        e$KA <- runif(1, 0.001, 1000);
        ka <- e$KA;
        k <- e$CL / e$V;
        V <- e$V;
        k12 <- e$Q / e$V;
        k21 <- e$Q / e$V2;
        beta <- 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k) * (k12 + k21 + k) - 4.0 * k21 * k));
        alpha <- k21 * k / beta;
        A <- (alpha - k21) / (alpha - beta) / V;
        B <- (beta - k21) / (beta - alpha) / V;
        getMacroConstants(e);
        expect_equal(e$A, A)
        expect_equal(e$alpha, alpha);
        expect_equal(e$B, B)
        expect_equal(e$beta, beta);
        g <- matrix(c(alpha, ka / (ka - alpha) * A,
                      beta, ka / (ka - beta) * B),
                    nrow=2, byrow=TRUE);
        expect_equal(e$g, g);
        expect_equal(e$oral, 1L);
        expect_equal(e$parameterization, 1L);
        expect_equal(e$ncmt, 2L);
        getLinDerivs(e);
        expect_true(all(ls(e)[regexpr("^d",ls(e))!= -1] %in% c("dCL", "dV", "dQ", "dV2", "dKA")));
    }
});

context("Test 2 compartment parameterizations (k, V, k12, k21)")

test_that("2 cmt k, V, k12, k21", {
    e <- new.env(parent = emptyenv());
    for (i  in 1:25){
        e$K <- runif(1, 0.001, 1000)
        e$V <- runif(1, 0.001, 1000)
        e$K12 <- runif(1, 0.001, 1000)
        e$K21 <- runif(1, 0.001, 1000)
        k <- e$K;
        V <- e$V;
        k12 <- e$K12
        k21 <- e$K21;
        beta <- 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k) * (k12 + k21 + k) - 4.0 * k21 * k));
        alpha <- k21 * k / beta;
        A <- (alpha - k21) / (alpha - beta) / V;
        B <- (beta - k21) / (beta - alpha) / V;
        getMacroConstants(e);
        expect_equal(e$A, A)
        expect_equal(e$alpha, alpha);
        expect_equal(e$B, B)
        expect_equal(e$beta, beta);
        g <- matrix(c(alpha, A,
                      beta, B),nrow=2, byrow=TRUE);
        expect_equal(e$g, g);
        expect_equal(e$oral, 0L);
        expect_equal(e$parameterization, 2L);
        expect_equal(e$ncmt, 2L);
        getLinDerivs(e);
        expect_true(all(ls(e)[regexpr("^d",ls(e))!= -1] %in% c("dK", "dV", "dK12", "dK21")));
    }
})
context("\tOral");
test_that("oral 2 cmt k, V, k12, k21", {
    e <- new.env(parent = emptyenv());
    for (i  in 1:25){
        e$K <- runif(1, 0.001, 1000);
        e$V <- runif(1, 0.001, 1000);
        e$K12 <- runif(1, 0.001, 1000);
        e$K21 <- runif(1, 0.001, 1000);
        e$KA <- runif(1, 0.0001, 1000);
        ka <- e$KA;
        k <- e$K;
        V <- e$V;
        k12 <- e$K12
        k21 <- e$K21;
        beta <- 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k) * (k12 + k21 + k) - 4.0 * k21 * k));
        alpha <- k21 * k / beta;
        A <- (alpha - k21) / (alpha - beta) / V;
        B <- (beta - k21) / (beta - alpha) / V;
        getMacroConstants(e);
        expect_equal(e$A, A)
        expect_equal(e$alpha, alpha);
        expect_equal(e$B, B)
        expect_equal(e$beta, beta);
        g <- matrix(c(alpha, ka / (ka - alpha) * A,
                      beta, ka / (ka - beta) * B),
                    nrow=2, byrow=TRUE);
        expect_equal(e$g, g);
        expect_equal(e$oral, 1L);
        expect_equal(e$parameterization, 2L);
        expect_equal(e$ncmt, 2L);
        getLinDerivs(e);
        expect_true(all(ls(e)[regexpr("^d",ls(e))!= -1] %in% c("dK", "dV", "dK12", "dK21", "dKA")));
    }
})


context("Test 2 compartment parameterizations (Cl, V, Q, VSS)")

test_that("2 cmt CL V Q, VSS", {
    e <- new.env(parent = emptyenv());
    for (i  in 1:25){
        e$CL <- runif(1, 0.001, 1000)
        e$V <- runif(1, 0.001, 1000)
        e$Q <- runif(1, 0.001, 1000)
        e$VSS <- e$V + runif(1, 0.001, 1000)
        V <- e$V;
        k <- e$CL / e$V;
        k12 <- e$Q / V;
        k21 <- e$Q / (e$VSS-V);
        beta <- 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k) * (k12 + k21 + k) - 4.0 * k21 * k));
        alpha <- k21 * k / beta;
        A <- (alpha - k21) / (alpha - beta) / V;
        B <- (beta - k21) / (beta - alpha) / V;
        getMacroConstants(e);
        expect_equal(e$A, A)
        expect_equal(e$alpha, alpha);
        expect_equal(e$B, B)
        expect_equal(e$beta, beta);
        g <- matrix(c(alpha, A,
                      beta, B),nrow=2, byrow=TRUE);
        expect_equal(e$g, g);
        expect_equal(e$oral, 0L);
        expect_equal(e$parameterization, 3L);
        expect_equal(e$ncmt, 2L);
        getLinDerivs(e);
        expect_true(all(ls(e)[regexpr("^d",ls(e))!= -1] %in% c("dCL", "dV", "dQ", "dVSS")));
    }
})

context("\tOral");
test_that("oral 2 cmt CL V Q, VSS", {
    e <- new.env(parent = emptyenv());
    for (i  in 1:25){
        e$CL <- runif(1, 0.001, 1000)
        e$V <- runif(1, 0.001, 1000)
        e$Q <- runif(1, 0.001, 1000)
        e$VSS <- e$V + runif(1, 0.001, 1000)
        e$KA <- runif(1, 0.001, 1000)
        ka <- e$KA;
        V <- e$V;
        k <- e$CL / e$V;
        k12 <- e$Q / V;
        k21 <- e$Q / (e$VSS-V);
        beta <- 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k) * (k12 + k21 + k) - 4.0 * k21 * k));
        alpha <- k21 * k / beta;
        A <- (alpha - k21) / (alpha - beta) / V;
        B <- (beta - k21) / (beta - alpha) / V;
        getMacroConstants(e);
        expect_equal(e$A, A)
        expect_equal(e$alpha, alpha);
        expect_equal(e$B, B)
        expect_equal(e$beta, beta);
        g <- matrix(c(alpha, ka / (ka - alpha) * A,
                      beta, ka / (ka - beta) * B),
                    nrow=2, byrow=TRUE);
        expect_equal(e$g, g);
        expect_equal(e$oral, 1L);
        expect_equal(e$parameterization, 3L);
        expect_equal(e$ncmt, 2L);
        getLinDerivs(e);
        expect_true(all(ls(e)[regexpr("^d",ls(e))!= -1] %in% c("dCL", "dV", "dQ", "dVSS", "dKA")));
    }
})


context("Test 2 compartment parameterizations (ALPHA, BETA, AOB, V)")

test_that("2 cmt ALPHA, BETA, AOB, V", {
    e <- new.env(parent = emptyenv());
    for (i  in 1:25){
        e$ALPHA <- runif(1, 0.001, 1000)
        e$BETA <- runif(1, 0.001, 1000)
        e$AOB <- runif(1, 0.001, 1000);
        e$V <- runif(1, 0.001, 1000);
        aob <- e$AOB;
        alpha <- e$ALPHA;
        beta <- e$BETA;
        k21 <- (aob*beta+alpha)/(aob+1);
        k <- (alpha*beta)/k21;
        k12 <- alpha+beta-k21-k;
        V <- e$V;
        beta <- 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k) * (k12 + k21 + k) - 4.0 * k21 * k));
        alpha <- k21 * k / beta;
        A <- (alpha - k21) / (alpha - beta) / V;
        B <- (beta - k21) / (beta - alpha) / V;
        getMacroConstants(e);
        expect_equal(e$A, A)
        expect_equal(e$alpha, alpha);
        expect_equal(e$B, B)
        expect_equal(e$beta, beta);
        g <- matrix(c(alpha, A,
                      beta, B),nrow=2, byrow=TRUE);
        expect_equal(e$g, g);
        expect_equal(e$oral, 0L);
        expect_equal(e$parameterization, 4L);
        expect_equal(e$ncmt, 2L);
        getLinDerivs(e);
        expect_true(all(ls(e)[regexpr("^d",ls(e))!= -1] %in% c("dALPHA", "dBETA", "dAOB", "dV")));
    }
})

context("\tOral");
test_that("2 cmt ALPHA, BETA, AOB, V, KA", {
    e <- new.env(parent = emptyenv());
    for (i  in 1:25){
        e$ALPHA <- runif(1, 0.001, 1000)
        e$BETA <- runif(1, 0.001, 1000)
        e$AOB <- runif(1, 0.001, 1000);
        e$V <- runif(1, 0.001, 1000);
        e$KA <- runif(1, 0.001, 1000)
        ka <- e$KA;
        aob <- e$AOB;
        alpha <- e$ALPHA;
        beta <- e$BETA;
        k21 <- (aob*beta+alpha)/(aob+1);
        k <- (alpha*beta)/k21;
        k12 <- alpha+beta-k21-k;
        V <- e$V;
        beta <- 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k) * (k12 + k21 + k) - 4.0 * k21 * k));
        alpha <- k21 * k / beta;
        A <- (alpha - k21) / (alpha - beta) / V;
        B <- (beta - k21) / (beta - alpha) / V;
        getMacroConstants(e);
        expect_equal(e$A, A)
        expect_equal(e$alpha, alpha);
        expect_equal(e$B, B)
        expect_equal(e$beta, beta);
        g <- matrix(c(alpha, ka / (ka - alpha) * A,
                      beta, ka / (ka - beta) * B),
                    nrow=2, byrow=TRUE);
        expect_equal(e$g, g);
        expect_equal(e$oral, 1L);
        expect_equal(e$parameterization, 4L);
        expect_equal(e$ncmt, 2L);
        getLinDerivs(e);
        expect_true(all(ls(e)[regexpr("^d",ls(e))!= -1] %in% c("dALPHA", "dBETA", "dAOB", "dV", "dKA")));
    }
})


context("Test 2 compartment parameterizations (ALPHA, BETA, K21, V)")

test_that("2 cmt ALPHA, BETA, AOB, V", {
    e <- new.env(parent = emptyenv());
    for (i  in 1:25){
        e$ALPHA <- runif(1, 0.001, 1000)
        e$BETA <- runif(1, 0.001, 1000)
        e$K21 <- runif(1, 0.001, 1000);
        e$V <- runif(1, 0.001, 1000);
        k21 = e$K21;
        alpha <- e$ALPHA;
        beta <- e$BETA;
        k <- (alpha*beta)/k21;
        k12 <- alpha+beta-k21-k;
        V <- e$V;
        beta <- 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k) * (k12 + k21 + k) - 4.0 * k21 * k));
        alpha <- k21 * k / beta;
        A <- (alpha - k21) / (alpha - beta) / V;
        B <- (beta - k21) / (beta - alpha) / V;
        getMacroConstants(e);
        expect_equal(e$A, A)
        expect_equal(e$alpha, alpha);
        expect_equal(e$B, B)
        expect_equal(e$beta, beta);
        g <- matrix(c(alpha, A,
                      beta, B),nrow=2, byrow=TRUE);
        expect_equal(e$g, g);
        expect_equal(e$oral, 0L);
        expect_equal(e$parameterization, 5L);
        expect_equal(e$ncmt, 2L);
        getLinDerivs(e);
        expect_true(all(ls(e)[regexpr("^d",ls(e))!= -1] %in% c("dALPHA", "dBETA", "dK21", "dV")));
    }
})

context("\tOral")
test_that("oral 2 cmt ALPHA, BETA, k21, V", {
    e <- new.env(parent = emptyenv());
    for (i  in 1:25){
        e$ALPHA <- runif(1, 0.001, 1000)
        e$BETA <- runif(1, 0.001, 1000)
        e$K21 <- runif(1, 0.001, 1000);
        e$V <- runif(1, 0.001, 1000);
        e$KA <- runif(1, 0.001, 1000);
        ka <- e$KA;
        k21 = e$K21;
        alpha <- e$ALPHA;
        beta <- e$BETA;
        k <- (alpha*beta)/k21;
        k12 <- alpha+beta-k21-k;
        V <- e$V;
        beta <- 0.5 * (k12 + k21 + k - sqrt((k12 + k21 + k) * (k12 + k21 + k) - 4.0 * k21 * k));
        alpha <- k21 * k / beta;
        A <- (alpha - k21) / (alpha - beta) / V;
        B <- (beta - k21) / (beta - alpha) / V;
        getMacroConstants(e);
        expect_equal(e$A, A)
        expect_equal(e$alpha, alpha);
        expect_equal(e$B, B)
        expect_equal(e$beta, beta);
        g <- matrix(c(alpha, ka / (ka - alpha) * A,
                      beta, ka / (ka - beta) * B),
                    nrow=2, byrow=TRUE);
        expect_equal(e$g, g);
        expect_equal(e$oral, 1L);
        expect_equal(e$parameterization, 5L);
        expect_equal(e$ncmt, 2L);
        getLinDerivs(e);
        expect_true(all(ls(e)[regexpr("^d",ls(e))!= -1] %in% c("dALPHA", "dBETA", "dK21", "dV", "dKA")));
    }
})

## parameterization: CL V Q V2 Q2 V3

context("3 Compartmentment: CL V Q V2 Q2 V3")
test_that("3 cmt CL V Q, V2", {
    e <- new.env(parent = emptyenv());
    for (i  in 1:25){
        e$CL <- runif(1, 0.001, 1000)
        e$V <- runif(1, 0.001, 1000);
        e$Q <- runif(1, 0.001, 1000);
        e$V2 <- runif(1, 0.001, 1000);
        e$Q2 <- runif(1, 0.001, 1000);
        e$V3 <- runif(1, 0.001, 1000);
        V <- e$V;
        k <- e$CL / V;
        k12 <- e$Q / V;
        k21 <- e$Q / e$V2;
        k13 <- e$Q2 / V;
        k31 <- e$Q2 / e$V3;
        a0 <- k * k21 * k31;
        a1 <- k * k31 + k21 * k31 + k21 * k13 + k * k21 + k31 * k12;
        a2 <- k + k12 + k13 + k21 + k31;
        p <- a1 - a2 * a2 / 3.0;
        q <- 2.0 * a2 * a2 * a2 / 27.0 - a1 * a2 /3.0 + a0;
        r1 <- sqrt(-p * p * p / 27.0);
        r2 <- 2 *r1 ^ (1.0 / 3.0);
        theta <- acos(-q / (2.0 * r1)) / 3.0;
        alpha <- -(cos(theta) * r2 - a2 / 3.0);
        bet <- -(cos(theta + 2.0 / 3.0 * pi) * r2 - a2 / 3.0);
        gam <- -(cos(theta + 4.0 / 3.0 * pi) * r2 - a2 / 3.0);
        A <- (k21 - alpha) * (k31 - alpha) / (alpha - bet) / (alpha - gam) / V;
        B <- (k21 - bet) * (k31 - bet) / (bet - alpha) / (bet - gam) / V;
        C <- (k21 - gam) * (k31 - gam) / (gam - alpha) / (gam - bet) / V;
        getMacroConstants(e);
        expect_equal(e$A, A)
        expect_equal(e$alpha, alpha);
        expect_equal(e$B, B)
        expect_equal(e$beta, bet);
        expect_equal(e$C, C)
        expect_equal(e$gamma, gam);
        g <- matrix(c(alpha,  A,
                      bet,  B,
                      gam, C),
                    nrow=3, byrow=TRUE);
        expect_equal(e$g, g);
        expect_equal(e$oral, 0L);
        expect_equal(e$parameterization, 1L);
        expect_equal(e$ncmt, 3L);
        getLinDerivs(e);
        expect_true(all(ls(e)[regexpr("^d",ls(e))!= -1] %in% c("dCL", "dV", "dQ", "dV2", "dQ2", "dV3")));
    }
})

context("\tOral:");

test_that("oral 3 cmt CL V Q, V2", {
    e <- new.env(parent = emptyenv());
    for (i  in 1:25){
        e$CL <- runif(1, 0.001, 1000)
        e$V <- runif(1, 0.001, 1000);
        e$Q <- runif(1, 0.001, 1000);
        e$V2 <- runif(1, 0.001, 1000);
        e$Q2 <- runif(1, 0.001, 1000);
        e$V3 <- runif(1, 0.001, 1000);
        e$KA <- runif(1, 0.001, 1000);
        ka <- e$KA;
        V <- e$V;
        k <- e$CL / V;
        k12 <- e$Q / V;
        k21 <- e$Q / e$V2;
        k13 <- e$Q2 / V;
        k31 <- e$Q2 / e$V3;
        a0 <- k * k21 * k31;
        a1 <- k * k31 + k21 * k31 + k21 * k13 + k * k21 + k31 * k12;
        a2 <- k + k12 + k13 + k21 + k31;
        p <- a1 - a2 * a2 / 3.0;
        q <- 2.0 * a2 * a2 * a2 / 27.0 - a1 * a2 /3.0 + a0;
        r1 <- sqrt(-p * p * p / 27.0);
        r2 <- 2 *r1 ^ (1.0 / 3.0);
        theta <- acos(-q / (2.0 * r1)) / 3.0;
        alpha <- -(cos(theta) * r2 - a2 / 3.0);
        bet <- -(cos(theta + 2.0 / 3.0 * pi) * r2 - a2 / 3.0);
        gam <- -(cos(theta + 4.0 / 3.0 * pi) * r2 - a2 / 3.0);
        A <- (k21 - alpha) * (k31 - alpha) / (alpha - bet) / (alpha - gam) / V;
        B <- (k21 - bet) * (k31 - bet) / (bet - alpha) / (bet - gam) / V;
        C <- (k21 - gam) * (k31 - gam) / (gam - alpha) / (gam - bet) / V;
        getMacroConstants(e);
        expect_equal(e$A, A)
        expect_equal(e$alpha, alpha);
        expect_equal(e$B, B)
        expect_equal(e$beta, bet);
        expect_equal(e$C, C)
        expect_equal(e$gamma, gam);
        g <- matrix(c(alpha,  A  * ka / (ka - alpha),
                    bet,  B * ka / (ka - bet),
                    gam, C * ka / (ka - gam)),
                    nrow=3, byrow=TRUE);
        expect_equal(e$g, g);
        expect_equal(e$oral, 1L);
        expect_equal(e$parameterization, 1L);
        expect_equal(e$ncmt, 3L);
        getLinDerivs(e);
        expect_true(all(ls(e)[regexpr("^d",ls(e))!= -1] %in% c("dCL", "dV", "dQ", "dV2", "dQ2", "dV3", "dKA")));
    }
})


context("3 Compartmentment: V K12 K21 K13 K31")
test_that("3 cmt V K12 K21 K13 K31", {
    e <- new.env(parent = emptyenv());
    for (i  in 1:25){
        e$K <- runif(1, 0.001, 1000)
        e$V <- runif(1, 0.001, 1000);
        e$K12 <- runif(1, 0.001, 1000);
        e$K21 <- runif(1, 0.001, 1000);
        e$K31 <- runif(1, 0.001, 1000);
        e$K13 <- runif(1, 0.001, 1000);
        V <- e$V;
        k <- e$K
        k12 <- e$K12
        k21 <- e$K21
        k13 <- e$K13
        k31 <- e$K31
        a0 <- k * k21 * k31;
        a1 <- k * k31 + k21 * k31 + k21 * k13 + k * k21 + k31 * k12;
        a2 <- k + k12 + k13 + k21 + k31;
        p <- a1 - a2 * a2 / 3.0;
        q <- 2.0 * a2 * a2 * a2 / 27.0 - a1 * a2 /3.0 + a0;
        r1 <- sqrt(-p * p * p / 27.0);
        r2 <- 2 *r1 ^ (1.0 / 3.0);
        theta <- acos(-q / (2.0 * r1)) / 3.0;
        alpha <- -(cos(theta) * r2 - a2 / 3.0);
        bet <- -(cos(theta + 2.0 / 3.0 * pi) * r2 - a2 / 3.0);
        gam <- -(cos(theta + 4.0 / 3.0 * pi) * r2 - a2 / 3.0);
        A <- (k21 - alpha) * (k31 - alpha) / (alpha - bet) / (alpha - gam) / V;
        B <- (k21 - bet) * (k31 - bet) / (bet - alpha) / (bet - gam) / V;
        C <- (k21 - gam) * (k31 - gam) / (gam - alpha) / (gam - bet) / V;
        getMacroConstants(e);
        expect_equal(e$A, A)
        expect_equal(e$alpha, alpha);
        expect_equal(e$B, B)
        expect_equal(e$beta, bet);
        expect_equal(e$C, C)
        expect_equal(e$gamma, gam);
        g <- matrix(c(alpha,  A,
                      bet,  B,
                      gam, C),
                    nrow=3, byrow=TRUE);
        expect_equal(e$g, g);
        expect_equal(e$oral, 0L);
        expect_equal(e$parameterization, 2L);
        expect_equal(e$ncmt, 3L);
        getLinDerivs(e);
        expect_true(all(ls(e)[regexpr("^d",ls(e))!= -1] %in% c("dV", "dK12", "dK21", "dK13", "dK31", "dK")));
    }
})
context("\tOral")
test_that("3 cmt V K12 K21 K13 K31", {
    e <- new.env(parent = emptyenv());
    for (i  in 1:25){
        e$K <- runif(1, 0.001, 1000)
        e$V <- runif(1, 0.001, 1000);
        e$K12 <- runif(1, 0.001, 1000);
        e$K21 <- runif(1, 0.001, 1000);
        e$K31 <- runif(1, 0.001, 1000);
        e$K13 <- runif(1, 0.001, 1000);
        e$KA <-  runif(1, 0.001, 1000);
        V <- e$V;
        k <- e$K
        k12 <- e$K12
        k21 <- e$K21
        k13 <- e$K13
        k31 <- e$K31
        ka <- e$KA
        a0 <- k * k21 * k31;
        a1 <- k * k31 + k21 * k31 + k21 * k13 + k * k21 + k31 * k12;
        a2 <- k + k12 + k13 + k21 + k31;
        p <- a1 - a2 * a2 / 3.0;
        q <- 2.0 * a2 * a2 * a2 / 27.0 - a1 * a2 /3.0 + a0;
        r1 <- sqrt(-p * p * p / 27.0);
        r2 <- 2 *r1 ^ (1.0 / 3.0);
        theta <- acos(-q / (2.0 * r1)) / 3.0;
        alpha <- -(cos(theta) * r2 - a2 / 3.0);
        bet <- -(cos(theta + 2.0 / 3.0 * pi) * r2 - a2 / 3.0);
        gam <- -(cos(theta + 4.0 / 3.0 * pi) * r2 - a2 / 3.0);
        A <- (k21 - alpha) * (k31 - alpha) / (alpha - bet) / (alpha - gam) / V;
        B <- (k21 - bet) * (k31 - bet) / (bet - alpha) / (bet - gam) / V;
        C <- (k21 - gam) * (k31 - gam) / (gam - alpha) / (gam - bet) / V;
        getMacroConstants(e);
        expect_equal(e$A, A)
        expect_equal(e$alpha, alpha);
        expect_equal(e$B, B)
        expect_equal(e$beta, bet);
        expect_equal(e$C, C)
        expect_equal(e$gamma, gam);
        g <- matrix(c(alpha,  A  * ka / (ka - alpha),
                      bet,  B * ka / (ka - bet),
                      gam, C * ka / (ka - gam)),
                    nrow=3, byrow=TRUE)
        expect_equal(e$g, g);
        expect_equal(e$oral, 1L);
        expect_equal(e$parameterization, 2L);
        expect_equal(e$ncmt, 3L);
        getLinDerivs(e);
        expect_true(all(ls(e)[regexpr("^d",ls(e))!= -1] %in% c("dV", "dK12", "dK21", "dK13", "dK31", "dK", "dKA")));
    }
})
