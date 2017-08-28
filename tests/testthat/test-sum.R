rxPermissive({
    context("Test PythonFsum")
    test_that("Fsum", {
        expect_equal(rxPythonFsum(c()), 0);
        expect_equal(rxPythonFsum(c(0)), 0);
        expect_equal(rxPythonFsum(c(1e100, 1.0, -1e100, 1e-100, 1e50, -1.0, -1e50)), 1e-100);
        expect_equal(rxPythonFsum(c(1e100, 1.0, -1e100, 1e-100, 1e50, -1.0, -1e50)), 1e-100);
        expect_equal(rxPythonFsum(c(2.0^53, -0.5, -2.0^-54)), 2.0^53-1.0);
        expect_equal(rxPythonFsum(c(2.0 ^ 53, 1.0, 2.0 ^ -100)), 2.0 ^ 53+2.0)
        expect_equal(rxPythonFsum(c(2.0 ^ 53+10.0, 1.0, 2.0 ^ -100)), 2.0 ^ 53+12.0);
        expect_equal(rxPythonFsum(c(2.0 ^ 53-4.0, 0.5, 2.0 ^ -54)), 2.0 ^ 53-3.0);
        expect_equal(rxPythonFsum(c(1e16, 1., 1e-16)), 10000000000000002.0);
        expect_equal(rxPythonFsum(c(1e16-2., 1.-2. ^ -53, -(1e16-2.), -(1.-2. ^ -53))), 0.0)
        ## FIXME include these python tests...
        ## ([1./n for n in range(1, 1001)],
        ##     float.fromhex('0x1.df11f45f4e61ap+2')),
        ## ([(-1.)**n/n for n in range(1, 1001)],
        ##     float.fromhex('-0x1.62a2af1bd3624p-1')),
        ## ([1.7**(i+1)-1.7**i for i in range(1000)] + [-1.7**1000], -1.0),
                                        # exercise code for resizing partials array
        ## ([2.**n - 2.**(n+50) + 2.**(n+52) for n in range(-1074, 972, 2)] +
        ##  [-2.**1022],
        ##     float.fromhex('0x1.5555555555555p+970')),
    })
})
