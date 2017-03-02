context("Test Symbolic inverse by sympy. (Needed for Almquist2015)")
rxPermissive({
    rmat <- function(dim){
        mat <- round(matrix(runif(dim ^ 2), dim, dim));
        mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)];
        diag(mat) <- 0;
        diag(mat) <- 1;
        return(mat);
    }

    mat2 <- matrix(c(1, 0.5, 0.5, 1), 2);
    symo <- rxSymInvCreate(mat2);

    test_that("Creates the correct matricies,", {
        expect_equal(digest::digest(as.list(symo %>% rxSymInv(1:3))),
                     "8bdf60e2cb4048314b95158f77727885");
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
        f <- tempfile();
        sink(f)
        print(symo)
        sink();
        tmp <- readLines(f);
        unlink(f);
        expect_equal(tmp, c("Object to create Omega and Omega^-1 & derivitaves for a 2x2 matrix:", "     [,1]   [,2]  ", "[1,] \"t0^2\" \"t1\"  ", "[2,] \"t1\"   \"t2^2\"", "Use `rxSymInv' for the matrix."));
    })

    test_that("Random inverses make sense", {
        for (i in 1:10){
            mat <- rmat(round(runif(1,0.5,3.5)));
            symo <- rxSymInvCreate(mat);
            nt <- symo$fn(as.double(0), as.integer(0), as.integer(1));
            for (j in 1:10){
                par <- abs(rnorm(nt))
                m <- symo %>% rxSymInv(par, 1);
                m1 <- solve(m);
                expect_equal(m1, symo %>% rxSymInv(par, -1));
            }
        }
    })

}, silent=TRUE)
