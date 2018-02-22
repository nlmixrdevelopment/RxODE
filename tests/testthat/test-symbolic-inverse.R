context("Test Symbolic inverse by sympy. (Needed for Almquist2015)")
rxPermissive({

    rmat <- function(dim){
        mat <- round(matrix(runif(dim ^ 2), dim, dim));
        mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)];
        diag(mat) <- 0;
        diag(mat) <- 1;
        return(mat);
    }

    mat <- diag(c(0.1, 0.1))
    symo <- rxSymInvCreate(mat);

    test_that("Creates the correct log.det.OmegaInv.5", {
        expect_equal(rxSymInv(symo, sqrt(c(0.1, 0.1)))$log.det.OMGAinv.5, 2.30258509299405)
    })

    mat <- diag(c(0.1, 0.1))
    symo <- rxSymInvCreate(mat, chol=TRUE);

    test_that("Creates the correct log.det.OmegaInv.5", {
        expect_equal(rxSymInv(symo, symo$th)$log.det.OMGAinv.5, 2.30258509299405)
    })

    mat2 <- matrix(c(1, 0.5, 0.5, 1), 2);
    symo <- rxSymInvCreate(mat2);

    test_that("Creates the correct matricies,", {
        expect_equal(digest::digest(as.list(symo %>% rxSymInv(1:3))),
                     "3fe4583876c9949f146ec3658d72d24c");
        expect_equal(symo %>% rxSymInv(1:3, pow=1), structure(c(1, 2, 2, 9), .Dim = c(2L, 2L)));
        expect_equal(symo %>% rxSymInv(1:3, pow=1, 1), structure(c(2, 0, 0, 0), .Dim = c(2L, 2L)));
        expect_equal(symo %>% rxSymInv(1:3, pow=1, 2), structure(c(0, 1, 1, 0), .Dim = c(2L, 2L)));
        expect_equal(symo %>% rxSymInv(1:3, pow=1, 3), structure(c(0, 0, 0, 6), .Dim = c(2L, 2L)));
        expect_error(symo %>% rxSymInv(1:3, pow=1, 4))
        expect_error(symo %>% rxSymInv(1:3, pow=1, 0.5))
        expect_error(symo %>% rxSymInv(1:3, pow=1, -1))
        expect_error(symo %>% rxSymInv(1:3, pow=0.5, 3))
        expect_error(symo %>% rxSymInv(1:2))
        expect_error(symo %>% rxSymInv(1:4))
        expect_equal(symo %>% rxSymInv(1:3, pow= -1), structure(c(1.8, -0.4, -0.4, 0.2), .Dim = c(2L, 2L)));
        expect_equal(symo %>% rxSymInv(c(1, 0.5, 2), -1, 0.5), 0.5 * log(det(symo %>% rxSymInv(c(1, 0.5, 2), -1))))
        ## f <- tempfile();
        ## sink(f)
        ## print(symo)
        ## sink();
        ## tmp <- readLines(f);
        ## unlink(f);
        ## expect_equal(tmp, c("Object to create Omega and Omega^-1 & derivitaves for a 2x2 matrix:", "     [,1]   [,2]  ", "[1,] \"t0^2\" \"t1\"  ", "[2,] \"t1\"   \"t2^2\"", "Use `rxSymInv' for the matrix."));
    })

    mat2 <- matrix(c(1, -0.5, -0.5, 1), 2);
    symo <- rxSymInvCreate(mat2);

    test_that("Creates the correct matrices,", {
        expect_equal(digest::digest(as.list(symo %>% rxSymInv(1:3))),
                     "3fe4583876c9949f146ec3658d72d24c");
        expect_equal(symo %>% rxSymInv(1:3, pow=1), structure(c(1, 2, 2, 9), .Dim = c(2L, 2L)));
        expect_equal(symo %>% rxSymInv(1:3, pow=1, 1), structure(c(2, 0, 0, 0), .Dim = c(2L, 2L)));
        expect_equal(symo %>% rxSymInv(1:3, pow=1, 2), structure(c(0, 1, 1, 0), .Dim = c(2L, 2L)));
        expect_equal(symo %>% rxSymInv(1:3, pow=1, 3), structure(c(0, 0, 0, 6), .Dim = c(2L, 2L)));
        expect_error(symo %>% rxSymInv(1:3, pow=1, 4))
        expect_error(symo %>% rxSymInv(1:3, pow=1, 0.5))
        expect_error(symo %>% rxSymInv(1:3, pow=1, -1))
        expect_error(symo %>% rxSymInv(1:3, pow=0.5, 3))
        expect_error(symo %>% rxSymInv(1:2))
        expect_error(symo %>% rxSymInv(1:4))
        expect_equal(symo %>% rxSymInv(1:3, pow= -1), structure(c(1.8, -0.4, -0.4, 0.2), .Dim = c(2L, 2L)));
        expect_equal(symo %>% rxSymInv(c(1, 0.5, 2), -1, 0.5), 0.5 * log(det(symo %>% rxSymInv(c(1, 0.5, 2), -1))))
        ## f <- tempfile();
        ## sink(f)
        ## print(symo)
        ## sink();
        ## tmp <- readLines(f);
        ## unlink(f);
        ## expect_equal(tmp, c("Object to create Omega and Omega^-1 & derivitaves for a 2x2 matrix:", "     [,1]   [,2]  ", "[1,] \"t0^2\" \"t1\"  ", "[2,] \"t1\"   \"t2^2\"", "Use `rxSymInv' for the matrix."));
    })

    test_that("Random inverses make sense", {
        for (i in 1:11){
            mat <- rmat(round(runif(1,0.5,3.5)));
            symo <- rxSymInvCreate(mat);
            nt <- rxSymInv(symo, "n");
            for (j in 1:10){
                par <- abs(rnorm(nt))
                m <- symo %>% rxSymInv(par, 1);
                m1 <- solve(m);
                m <- symo %>% rxSymInv(par, -1)
                dimnames(m) <- list(NULL, NULL);
                dimnames(m1) <- list(NULL, NULL);
                expect_equal(m1, m);
            }
        }
    })

    test_that("diag 3x3 random inverses", {
        mat <- diag(3);
        symo <- rxSymInvCreate(mat);
        nt <- rxSymInv(symo, "n");
        for (j in 1:10){
            par <- abs(rnorm(nt))
            m <- symo %>% rxSymInv(par, 1);
            m1 <- solve(m);
            m <- symo %>% rxSymInv(par, -1);
            dimnames(m) <- list(NULL, NULL);
            dimnames(m1) <- list(NULL, NULL);
            expect_equal(m1, m);
        }
    })

    mat2 <- matrix(rep(1, 4), 2)
    mat2[rxBlockZeros(mat2, 1)] <- 0
    mat2i <- rxSymInvCreate(mat2);

    mat3.1 <- matrix(rep(1, 3 * 3), 3)
    mat3.1[rxBlockZeros(mat3.1, 1)] <- 0
    mat3.1i <- rxSymInvCreate(mat3.1);


    mat3.2 <- matrix(rep(1, 3 * 3), 3)
    mat3.2[rxBlockZeros(mat3.2, 2)] <- 0
    mat3.2i <- rxSymInvCreate(mat3.2);

    mat3 <- diag(3)
    mat3 <- rxSymInvCreate(mat3);


    mat4.1 <- matrix(rep(1, 4 * 4), 4)
    mat4.1[rxBlockZeros(mat4.1, 1)] <- 0
    mat4.1i <- rxSymInvCreate(mat4.1);


    mat4.2 <- matrix(rep(1, 4 * 4), 4)
    mat4.2[rxBlockZeros(mat4.2, 2)] <- 0
    mat4.2i <- rxSymInvCreate(mat4.2);


    mat4.3 <- matrix(rep(1, 4 * 4), 4)
    mat4.3[rxBlockZeros(mat4.3, 3)] <- 0
    mat4.3i <- rxSymInvCreate(mat4.3);

    mat4 <- diag(4);
    mat4 <- rxSymInvCreate(mat4);


    test_that("Symbolc inverses produce block matricies", {
        expect_equal(class(mat2i), "rxSymInvBlock");
        expect_equal(class(mat3.1i), "rxSymInvBlock");
        expect_equal(class(mat3.2i), "rxSymInvBlock");
        expect_equal(class(mat3), "rxSymInvBlock");
        expect_equal(class(mat4.1i), "rxSymInvBlock");
        expect_equal(class(mat4.2i), "rxSymInvBlock");
        expect_equal(class(mat4.3i), "rxSymInvBlock");
        expect_equal(class(mat4), "rxSymInvBlock");
    })

}, silent=TRUE, on.validate=TRUE)
