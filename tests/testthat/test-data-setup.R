rxPermissive({
    ## Test the behavior
    context("Test Data Setup (for RcppParallel-style for loop); 0 cov")
    library(dplyr);
    load(devtools::package_file("tests/testthat/test-data-setup.Rdata"))
    test_that("conversion without covariates", {
        convert1 <- rxDataSetup(dat);
        expect_equal(length(convert1$cov), 0);
        expect_equal(length(convert1$cov.names), 0);
        for (i in unique(dat$ID)){
            w <- seq(convert1$ids$posObs[i] + 1, convert1$ids$posObs[i] + convert1$ids$nObs[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID == 0,
                                                 c("DV", "TIME")])),
                         as.double(as.matrix(convert1$obs[w, ])))
            w <- seq(convert1$ids$posDose[i] + 1, convert1$ids$posDose[i] + convert1$ids$nDose[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID != 0,
                                                 c("EVID", "TIME", "AMT")])),
                         as.double(as.matrix(convert1$dose[w, ])))
            w <- seq(convert1$ids$posEvent[i] + 1, convert1$ids$posEvent[i] + convert1$ids$nEvent[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i, c("EVID", "TIME")])),
                         as.double(as.matrix(convert1$et[w, ])))
        }
        expect_false(convert1$missing.id)
        expect_false(convert1$missing.dv)
    })

    context("  matrix instead of data.frame")
    test_that("matrix", {
        convert1 <- rxDataSetup(as.matrix(dat));
        expect_equal(length(convert1$cov), 0);
        expect_equal(length(convert1$cov.names), 0);
        for (i in unique(dat$ID)){
            w <- seq(convert1$ids$posObs[i] + 1, convert1$ids$posObs[i] + convert1$ids$nObs[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID == 0,
                                                 c("DV", "TIME")])),
                         as.double(as.matrix(convert1$obs[w, ])))
            w <- seq(convert1$ids$posDose[i] + 1, convert1$ids$posDose[i] + convert1$ids$nDose[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID != 0,
                                                 c("EVID", "TIME", "AMT")])),
                         as.double(as.matrix(convert1$dose[w, ])))
            w <- seq(convert1$ids$posEvent[i] + 1, convert1$ids$posEvent[i] + convert1$ids$nEvent[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i, c("EVID", "TIME")])),
                         as.double(as.matrix(convert1$et[w, ])))
        }
        expect_false(convert1$missing.id)
        expect_false(convert1$missing.dv)
    })


    context("  tbl instead of data.frame")
    test_that("tbl", {
        convert1 <- rxDataSetup(as.tbl(dat));
        expect_equal(length(convert1$cov), 0);
        expect_equal(length(convert1$cov.names), 0);
        for (i in unique(dat$ID)){
            w <- seq(convert1$ids$posObs[i] + 1, convert1$ids$posObs[i] + convert1$ids$nObs[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID == 0,
                                                 c("DV", "TIME")])),
                         as.double(as.matrix(convert1$obs[w, ])))
            w <- seq(convert1$ids$posDose[i] + 1, convert1$ids$posDose[i] + convert1$ids$nDose[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID != 0,
                                                 c("EVID", "TIME", "AMT")])),
                         as.double(as.matrix(convert1$dose[w, ])))
            w <- seq(convert1$ids$posEvent[i] + 1, convert1$ids$posEvent[i] + convert1$ids$nEvent[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i, c("EVID", "TIME")])),
                         as.double(as.matrix(convert1$et[w, ])))
        }
        expect_false(convert1$missing.id)
        expect_false(convert1$missing.dv)
    })

    dat2 <- dat %>% mutate(id=ID, amt=AMT, time=TIME, evid=EVID, dv=DV) %>% select(-ID, -AMT, -TIME, -EVID, -DV);

    context("  lower case key column names")
    test_that("conversion without covariates; lower case names", {
        convert1 <- rxDataSetup(dat2);
        expect_equal(length(convert1$cov), 0);
        expect_equal(length(convert1$cov.names), 0);
        for (i in unique(dat$ID)){
            w <- seq(convert1$ids$posObs[i] + 1, convert1$ids$posObs[i] + convert1$ids$nObs[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID == 0,
                                                 c("DV", "TIME")])),
                         as.double(as.matrix(convert1$obs[w, ])))
            w <- seq(convert1$ids$posDose[i] + 1, convert1$ids$posDose[i] + convert1$ids$nDose[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID != 0,
                                                 c("EVID", "TIME", "AMT")])),
                         as.double(as.matrix(convert1$dose[w, ])))
            w <- seq(convert1$ids$posEvent[i] + 1, convert1$ids$posEvent[i] + convert1$ids$nEvent[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i, c("EVID", "TIME")])),
                         as.double(as.matrix(convert1$et[w, ])))
        }
        expect_false(convert1$missing.id)
        expect_false(convert1$missing.dv)
    })

    dat2 <- dat %>% mutate(Id=ID, Amt=AMT, Time=TIME, Evid=EVID, Dv=DV) %>% select(-ID, -AMT, -TIME, -EVID, -DV);
    context("  title case key column names")
    test_that("conversion without covariates; lower case names", {
        convert1 <- rxDataSetup(dat2);
        expect_equal(length(convert1$cov), 0);
        expect_equal(length(convert1$cov.names), 0);
        for (i in unique(dat$ID)){
            w <- seq(convert1$ids$posObs[i] + 1, convert1$ids$posObs[i] + convert1$ids$nObs[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID == 0,
                                                 c("DV", "TIME")])),
                         as.double(as.matrix(convert1$obs[w, ])))
            w <- seq(convert1$ids$posDose[i] + 1, convert1$ids$posDose[i] + convert1$ids$nDose[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID != 0,
                                                 c("EVID", "TIME", "AMT")])),
                         as.double(as.matrix(convert1$dose[w, ])))
            w <- seq(convert1$ids$posEvent[i] + 1, convert1$ids$posEvent[i] + convert1$ids$nEvent[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i, c("EVID", "TIME")])),
                         as.double(as.matrix(convert1$et[w, ])))
        }
        expect_false(convert1$missing.id)
        expect_false(convert1$missing.dv)
    })

    context(" STrange capitalization")
    dat2 <- dat %>% mutate(Id=ID, AMt=AMT, TIme=TIME, EVid=EVID, Dv=DV) %>% select(-ID, -AMT, -TIME, -EVID, -DV);
    test_that("bad setup", {
        expect_error(rxDataSetup(dat2))
    })

    dat2 <- dat %>% select(-DV)

    context("  missing DV")
    test_that("missing DV", {
        convert1 <- rxDataSetup(dat2);
        expect_equal(length(convert1$cov), 0);
        expect_equal(length(convert1$cov.names), 0);
        for (i in unique(dat$ID)){
            w <- seq(convert1$ids$posObs[i] + 1, convert1$ids$posObs[i] + convert1$ids$nObs[i])
            tmp <- dat[dat$ID == i & dat$EVID == 0,
                       c("DV", "TIME")]
            tmp$DV <- NA;
            expect_equal(as.double(as.matrix(tmp)),
                         as.double(as.matrix(convert1$obs[w, ])))
            w <- seq(convert1$ids$posDose[i] + 1, convert1$ids$posDose[i] + convert1$ids$nDose[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID != 0,
                                                 c("EVID", "TIME", "AMT")])),
                         as.double(as.matrix(convert1$dose[w, ])))
            w <- seq(convert1$ids$posEvent[i] + 1, convert1$ids$posEvent[i] + convert1$ids$nEvent[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i, c("EVID", "TIME")])),
                         as.double(as.matrix(convert1$et[w, ])))
        }
        expect_false(convert1$missing.id)
        expect_true(convert1$missing.dv)
    })

    dat2 <- dat %>% filter(ID == 1) %>% select(-ID)

    context("  missing ID")
    test_that("missing ID", {
        convert1 <- rxDataSetup(dat2);
        expect_equal(length(convert1$cov), 0);
        expect_equal(length(convert1$cov.names), 0);
        i <- 1;
        w <- seq(convert1$ids$posObs[i] + 1, convert1$ids$posObs[i] + convert1$ids$nObs[i])
        tmp <- dat2[ dat2$EVID == 0,
                   c("DV", "TIME")]
        expect_equal(as.double(as.matrix(tmp)),
                     as.double(as.matrix(convert1$obs[w, ])))
        w <- seq(convert1$ids$posDose[i] + 1, convert1$ids$posDose[i] + convert1$ids$nDose[i])
        expect_equal(as.double(as.matrix(dat2[dat2$EVID != 0,
                                             c("EVID", "TIME", "AMT")])),
                     as.double(as.matrix(convert1$dose[w, ])))
        w <- seq(convert1$ids$posEvent[i] + 1, convert1$ids$posEvent[i] + convert1$ids$nEvent[i])
        expect_equal(as.double(as.matrix(dat2[, c("EVID", "TIME")])),
                     as.double(as.matrix(convert1$et[w, ])))
        expect_true(convert1$missing.id)
        expect_false(convert1$missing.dv)
    })

    context("  missing DV/ID");
    dat2 <- dat %>% filter(ID == 1) %>% select(-ID, -DV)

    test_that("missing DV/ID", {
        convert1 <- rxDataSetup(dat2);
        expect_equal(length(convert1$cov), 0);
        expect_equal(length(convert1$cov.names), 0);
        i <- 1;
        w <- seq(convert1$ids$posObs[i] + 1, convert1$ids$posObs[i] + convert1$ids$nObs[i])
        dat2$DV <- NA
        tmp <- dat2[dat2$EVID == 0,
                    c("DV", "TIME")]
        tmp$DV <- NA;
        expect_equal(as.double(as.matrix(tmp)),
                     as.double(as.matrix(convert1$obs[w, ])))
        w <- seq(convert1$ids$posDose[i] + 1, convert1$ids$posDose[i] + convert1$ids$nDose[i])
        expect_equal(as.double(as.matrix(dat2[dat2$EVID != 0,
                                              c("EVID", "TIME", "AMT")])),
                     as.double(as.matrix(convert1$dose[w, ])))
        w <- seq(convert1$ids$posEvent[i] + 1, convert1$ids$posEvent[i] + convert1$ids$nEvent[i])
        expect_equal(as.double(as.matrix(dat2[, c("EVID", "TIME")])),
                     as.double(as.matrix(convert1$et[w, ])))
        expect_true(convert1$missing.id)
        expect_true(convert1$missing.dv)
    })

    dat2 <- dat %>% select(-ID);

    context("  unsorted ID/TIME throw error.")
    test_that("error on unsorted data", {
        expect_error(rxDataSetup(dat2));
    })

    cn <- c("V")
    convert2 <- rxDataSetup(dat, cn);
    context("Test Data Setup (for RcppParallel-style for loop); 1 cov")
    test_that("conversion with 1 covariate", {
        expect_equal(length(convert2$cov.names), 1);
        for (i in unique(dat$ID)){
            w <- seq(convert2$ids$posObs[i] + 1,
                     convert2$ids$posObs[i] + convert2$ids$nObs[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID == 0,
                                                 c("DV", "TIME")])),
                         as.double(as.matrix(convert2$obs[w, ])))
            w <- seq(convert2$ids$posDose[i] + 1,
                     convert2$ids$posDose[i] + convert2$ids$nDose[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID != 0,
                                                 c("EVID", "TIME", "AMT")])),
                         as.double(as.matrix(convert2$dose[w, ])));
            w <- seq(convert2$ids$posCov[i] + 1,
                     convert2$ids$posCov[i] + convert2$ids$nCov[i])
            expect_equal(as.double(convert2$cov[w]),
                         as.double(as.matrix(dat[dat$ID == i & dat$EVID == 0,cn])));
            w <- seq(convert2$ids$posEvent[i] + 1, convert2$ids$posEvent[i] + convert2$ids$nEvent[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i, c("EVID", "TIME")])),
                         as.double(as.matrix(convert2$et[w, ])))
        }
    })

    cn <- c("V", "CL")
    convert2 <- rxDataSetup(dat, cn);
    context("Test Data Setup (for RcppParallel-style for loop); 2 cov")
    test_that("conversion with 2 covariate", {
        expect_equal(length(convert2$cov.names), 2);
        for (i in unique(dat$ID)){
            w <- seq(convert2$ids$posObs[i] + 1,
                     convert2$ids$posObs[i] + convert2$ids$nObs[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID == 0,
                                                 c("DV", "TIME")])),
                         as.double(as.matrix(convert2$obs[w, ])))
            w <- seq(convert2$ids$posDose[i] + 1,
                     convert2$ids$posDose[i] + convert2$ids$nDose[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID != 0,
                                                 c("EVID", "TIME", "AMT")])),
                         as.double(as.matrix(convert2$dose[w, ])));
            w <- seq(convert2$ids$posCov[i] + 1,
                     convert2$ids$posCov[i] + convert2$ids$nCov[i])
            expect_equal(as.double(convert2$cov[w]),
                         as.double(as.matrix(dat[dat$ID == i & dat$EVID == 0,cn])));
            w <- seq(convert2$ids$posEvent[i] + 1, convert2$ids$posEvent[i] + convert2$ids$nEvent[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i, c("EVID", "TIME")])),
                         as.double(as.matrix(convert2$et[w, ])))
        }
    })

    cn <- c("V", "CL", "DOSE")
    convert2 <- rxDataSetup(dat, cn);
    context("Test Data Setup (for RcppParallel-style for loop); 3 cov")
    test_that("conversion with 3 covariate", {
        expect_equal(length(convert2$cov.names), 3);
        for (i in unique(dat$ID)){
            w <- seq(convert2$ids$posObs[i] + 1,
                     convert2$ids$posObs[i] + convert2$ids$nObs[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID == 0,
                                                 c("DV", "TIME")])),
                         as.double(as.matrix(convert2$obs[w, ])))
            w <- seq(convert2$ids$posDose[i] + 1,
                     convert2$ids$posDose[i] + convert2$ids$nDose[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i & dat$EVID != 0,
                                                 c("EVID", "TIME", "AMT")])),
                         as.double(as.matrix(convert2$dose[w, ])));
            w <- seq(convert2$ids$posCov[i] + 1,
                     convert2$ids$posCov[i] + convert2$ids$nCov[i])
            expect_equal(as.double(convert2$cov[w]),
                         as.double(as.matrix(dat[dat$ID == i & dat$EVID == 0,cn])));
            w <- seq(convert2$ids$posEvent[i] + 1, convert2$ids$posEvent[i] + convert2$ids$nEvent[i])
            expect_equal(as.double(as.matrix(dat[dat$ID == i, c("EVID", "TIME")])),
                         as.double(as.matrix(convert2$et[w, ])))
        }
    })

    context("Expand event table with covariate information")
    test_that("Setup Event table", {

        et <- eventTable()   # default time units
        et$add.sampling(seq(from=0, to=100, by=0.01))

        cov <- data.frame(c=et$get.sampling()$time+1, d=et$get.sampling()$time+1);

        tmp1 <- rxDataSetup(et, as.matrix(cov));
        tmp2 <- rxDataSetup(et, cov)

        expect_equal(tmp1, tmp2)

        cov2 <- data.frame(c=et$get.sampling()$time[-1]+1, d=et$get.sampling()$time[-1]+1);

        expect_error(rxDataSetup(et, as.matrix(cov2)));
        expect_error(rxDataSetup(et, cov2));

        cov2 <- data.frame(c=c(et$get.sampling()$time, 1)+1, d=c(et$get.sampling()$time, 1)+1);

    })

})
