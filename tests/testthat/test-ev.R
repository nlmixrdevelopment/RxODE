context("Test event Table et(...)");
rxPermissive({

    et <- et()

    test_that("Empty event table check",{

        expect_equal(et$nobs, 0L)
        expect_equal(et$ndose, 0L)
        expect_equal(et$get.EventTable(), NULL);

        expect_equal(et$get.dosing(), NULL);
        expect_equal(et$get.sampling(), NULL);
    })

    et <- et(1:10);

    test_that("Observation only table check",{

        expect_equal(et$nobs, 10L);
        expect_equal(et$ndose, 0L);

        expect_equal(class(et$get.EventTable()), "data.frame");
        expect_true(is(et, "rxEt"));

        expect_false(et$.show["id"])
        expect_false(et$.show["cmt"])


        expect_equal(et$get.dosing(), NULL);
        expect_equal(et$get.sampling(), NULL);
    })

    et <- et(1:10, "matt");

    test_that("Observation only table check",{

        expect_equal(et$nobs, 10L);
        expect_equal(et$ndose, 0L);

        expect_equal(class(et$get.EventTable()), "data.frame");
        expect_true(is(et, "rxEt"));

        expect_false(et$.show["id"])
        expect_true(et$.show["cmt"])

        expect_equal(et$get.dosing(), NULL);
        expect_equal(et$get.sampling(), NULL);
    })

    et1 <- et(1:10,id=10)

    test_that("Observation only table check",{

        expect_equal(et$nobs, 100L);
        expect_equal(et$ndose, 0L);

        expect_equal(class(et$get.EventTable()), "data.frame");
        expect_true(is(et, "rxEt"));

        expect_true(et$.show["id"])
        expect_false(et$.show["cmt"])

        expect_equal(et$get.dosing(), NULL);
        expect_equal(et$get.sampling(), NULL);
    })

    ## now resize down
    et2 <- et1 %>% et(id=1)


    ## now resize back up
    et3 <- et2 %>% et(id=10)

})
