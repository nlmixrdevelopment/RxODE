rxPermissive({
    context("nlmixr style matrix parsing");
    test_that("nlmixr matrix parsing", {
        expect_equal(rxMatrix({et2 + et3 + et4 ~ c(40,
                                                   0.1, 20,
                                                   0.1, 0.1, 30)}),
                     structure(c(40, 0.1, 0.1, 0.1, 20, 0.1, 0.1, 0.1, 30),
                               .Dim = c(3L, 3L),
                               .Dimnames = list(c("et2", "et3", "et4"),
                                                c("et2", "et3", "et4"))))

        expect_equal(rxMatrix(list(et2 + et3 + et4 ~ c(40,
                                                       0.1, 20,
                                                       0.1, 0.1, 30),
                                   matrix(1,dimnames=list("et5","et5")))),
                     structure(c(40, 0.1, 0.1, 0, 0.1, 20, 0.1, 0, 0.1, 0.1, 30, 0,
                                 0, 0, 0, 1),
                               .Dim = c(4L, 4L),
                               .Dimnames = list(c("et2", "et3", "et4", "et5"),
                                                c("et2", "et3", "et4", "et5"))))

        expect_equal(rxMatrix(list(et2 + et3 + et4 ~ c(40,
                                                       0.1, 20,
                                                       0.1, 0.1, 30),
                                   matrix(1,dimnames=list("et5","et5")))),
                     structure(c(40, 0.1, 0.1, 0, 0.1, 20, 0.1, 0, 0.1, 0.1, 30, 0,
                                 0, 0, 0, 1),
                               .Dim = c(4L, 4L),
                               .Dimnames = list(c("et2", "et3", "et4", "et5"),
                                                c("et2", "et3", "et4", "et5"))))

        expect_equal(rxMatrix({
            et2 + et3 + et4 ~ c(40,
                                0.1, 20,
                                0.1, 0.1, 30);
            et5 ~ 1;}),
            structure(c(40, 0.1, 0.1, 0, 0.1, 20, 0.1, 0, 0.1, 0.1, 30, 0,
                        0, 0, 0, 1),
                      .Dim = c(4L, 4L),
                      .Dimnames = list(c("et2", "et3", "et4", "et5"),
                                       c("et2", "et3", "et4", "et5"))))
        expect_equal(rxMatrix(et2 + et3 + et4 ~ c(40,
                                                  0.1, 20,
                                                  0.1, 0.1, 30),
                              et5 ~ 1),
                     structure(c(40, 0.1, 0.1, 0, 0.1, 20, 0.1, 0, 0.1, 0.1, 30, 0,
                                 0, 0, 0, 1),
                               .Dim = c(4L, 4L),
                               .Dimnames = list(c("et2", "et3", "et4", "et5"),
                                                c("et2", "et3", "et4", "et5"))))

        expect_equal(rxMatrix(et2 + et3 + et4 ~ c(40,
                                                  0.1, 20,
                                                  0.1, 0.1, 30),
                              list(et5 ~ 1, et6 ~ 3)),
                     structure(c(40, 0.1, 0.1, 0, 0, 0.1, 20, 0.1, 0, 0, 0.1, 0.1,
                                 30, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 3),
                               .Dim = c(5L, 5L),
                               .Dimnames = list(c("et2", "et3", "et4", "et5", "et6"),
                                                c("et2", "et3", "et4", "et5", "et6"))))
        expect_equal(rxMatrix(quote({et2 + et3 + et4 ~ c(40,
                                                         0.1, 20,
                                                         0.1, 0.1, 30)})),
                     structure(c(40, 0.1, 0.1, 0.1, 20, 0.1, 0.1, 0.1, 30),
                               .Dim = c(3L, 3L),
                               .Dimnames = list(c("et2", "et3", "et4"),
                                                c("et2", "et3", "et4"))))

        .mat  <- rxMatrix({et2 + et3 + et4 ~ c(40,
                                               0.1, 20,
                                               0.1, 0.1, 30)})
        ## Test for NSE issues
        expect_equal(.mat, rxMatrix(.mat))
        ## Test for NULL
        expect_equal(NULL, rxMatrix(NULL))
    })
})
