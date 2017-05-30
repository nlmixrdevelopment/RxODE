context("Test Printing and summary functions");
rxPermissive({

    ode <- "
   C2 = centr/V2;
   C3 = peri/V3;
   d/dt(depot) =-KA*depot;
   d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
   d/dt(peri)  =                    Q*C2 - Q*C3;
   d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
"
    m1 <- RxODE(model = ode)

    et1 <- eventTable(amount.units="ug", time.units = "hours")
    et1$add.dosing(dose=10000, nbr.doses=5, dosing.interval = 24)
    et1$add.sampling(0:24)
    et1$add.sampling(24);
    et1$add.sampling(seq(from = 24+8, to = 5*24, by = 8))

    pred <- predict(m1, params = c(KA=.291, CL=18.6, V2=40.2, Q=10.5, V3=297.0,
                                   Kin=1.0, Kout=1.0, EC50=200.0),
                    events = et1,
                    inits = c(0, 0, 0, 1));

    ode2 <- "
   KA = .291
   CL = 18.6
   V2 = 40.2
   Q  = 10.5
   V3 = 297.0
   Kin = 1.0
   Kout = 1.0
   EC50 = 200.0
   eff(0) = 1
   C2 = centr/V2;
   C3 = peri/V3;
   d/dt(depot) =-KA*depot;
   d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
   d/dt(peri)  =                    Q*C2 - Q*C3;
   d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
"
    m2 <- RxODE(model = ode2)
    pred2 <- predict(m2, params = c(KA=.3, CL=18.6, V2=40.2, Q=10.5, V3=297.0,
                                    Kin=1.0, Kout=1.0, EC50=200.0),
                     events = et1,
                     inits = c(0, 0, 0, 1));

    tmpsink <- function(){
        while (sink.number() != 0){
            sink();
        }
        tmpfile <- tempfile();
        sink(tmpfile);
        return(tmpfile);
    }

    test_that("Print for RxODE m1 works correctly",{
        ##
        tmpfile <- tmpsink();
        print(m1)
        sink();
        p1 <- readLines(tmpfile);
        unlink(tmpfile);
        expect_equal(p1,sprintf("RxODE model named \"%s\" (ready to run)",basename(getwd())));
    });

    test_that("Print m1$cmpMgr", {

        tmpfile <- tempfile()
        tmpfile <- tmpsink();
        print(m1$cmpMgr)
        sink();
        p1 <- readLines(tmpfile);
        unlink(tmpfile);
        expect_equal(p1,sprintf("RxCompilationManager for RxODE model '%s'",basename(getwd())));
    })

    test_that("Print m1$cmpMgr$rxDll", {
        tmpfile <- tmpsink();
        print(m1$cmpMgr$rxDll())
        sink();
        p1 <- readLines(tmpfile);
        unlink(tmpfile);
        expect_equal(p1,sprintf("RxODE dll named \"%s\" is loaded and ready to use.",basename(rxDll(m1))));
    })

    ##
    cf <- c("","User Supplied Parameters:",
            "  V2   V3   KA   CL    Q  Kin Kout EC50 ",
            "  NA   NA   NA   NA   NA   NA   NA   NA ",
            "",
            "User Initial Conditions:",
            "depot centr  peri   eff ",
            "    0     0     0     0 ",
            "",
            "Compartents:",
            "  cmt=1   cmt=2   cmt=3   cmt=4 ",
            "\"depot\" \"centr\"  \"peri\"   \"eff\" ");


    test_that("Print coef(m1)", {
        ##
        tmpfile <- tmpsink();
        print(coef(m1));
        sink();
        p1 <- readLines(tmpfile);
        unlink(tmpfile);
        expect_equal(p1,cf);
    })
    cf2 <- c("","User Supplied Parameters:",
             "     KA      CL      V2       Q      V3     Kin    Kout    EC50 ",
             "  0.291  18.600  40.200  10.500 297.000   1.000   1.000 200.000 ",
             "",
             "User Initial Conditions:",
             "depot centr  peri   eff ",
             "    0     0     0     1 ",
             "",
             "Compartents:",
             "  cmt=1   cmt=2   cmt=3   cmt=4 ",
             "\"depot\" \"centr\"  \"peri\"   \"eff\" ")

    test_that("Print coef(m2)", {
        tmpfile <- tmpsink();
        print(coef(m2));
        sink();
        p1 <- readLines(tmpfile);
        unlink(tmpfile);
    })

    test_that("Print coef(m1$cmpMgr)", {
        tmpfile <- tmpsink();
        print(coef(m1$cmpMgr));
        sink();
        p1 <- readLines(tmpfile);
        unlink(tmpfile);
        expect_equal(p1,cf)
    })

    cfp <- c("","User Supplied Parameters ($params):",
             "     V2      V3      KA      CL       Q     Kin    Kout    EC50 ",
             " 40.200 297.000   0.291  18.600  10.500   1.000   1.000 200.000 ",
             "",
             "User Initial Conditions ($state):",
             "depot centr  peri   eff ",
             "    0     0     0     1 ",
             "",
             "Compartents:",
             "  cmt=1   cmt=2   cmt=3   cmt=4 ",
             "\"depot\" \"centr\"  \"peri\"   \"eff\" ");

    test_that("Print coef(pred)", {
        tmpfile <- tmpsink();
        print(coef(pred));
        sink();
        p1 <- readLines(tmpfile);
        unlink(tmpfile);
        expect_equal(p1,cfp);
    })

    cfp2 <- c("","User Supplied Parameters ($params):",
              "   KA    CL    V2     Q    V3   Kin  Kout  EC50 ",
              "  0.3  18.6  40.2  10.5 297.0   1.0   1.0 200.0 ",
              "",
              "User Initial Conditions ($state):",
              "depot centr  peri   eff ",
              "    0     0     0     1 ",
              "",
              "Default parameter values:",
              "     KA      CL      V2       Q      V3     Kin    Kout    EC50     eff ",
              "  0.291  18.600  40.200  10.500 297.000   1.000   1.000 200.000   1.000 ",
              "",
              "Compartents:",
              "  cmt=1   cmt=2   cmt=3   cmt=4 ",
              "\"depot\" \"centr\"  \"peri\"   \"eff\" ");

    test_that("Print coef(pred2)", {
        tmpfile <- tmpsink();
        print(coef(pred2));
        sink();
        p1 <- readLines(tmpfile);
        unlink(tmpfile);
        expect_equal(p1,cfp2)
    })
    ##
    s.base <- c(sprintf("dll: %s",rxDll(m1)),
                cf,
                "",
                "Calculated Variables:",
                "[1] \"C2\" \"C3\"",
                "",
                "Model:",
                "",
                "   C2 = centr/V2;",
                "   C3 = peri/V3;",
                "   d/dt(depot) =-KA*depot;",
                "   d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;",
                "   d/dt(peri)  =                    Q*C2 - Q*C3;",
                "   d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;",
                "",
                "");
    test_that("Print summary(m1)", {
        tmpfile <- tmpsink();
        summary(m1);
        sink();
        p1 <- readLines(tmpfile);
        unlink(tmpfile);
        expect_equal(p1,c(sprintf("RxODE model named \"%s\" (ready to run)",basename(getwd())),
                          s.base))
    })

    test_that("Print summary(m1$cmpMgr)", {
        tmpfile <- tmpsink();
        summary(m1$cmpMgr);
        sink();
        p1 <- readLines(tmpfile);
        unlink(tmpfile);
        expect_equal(p1,c(sprintf("RxCompilationManager for RxODE model '%s'",basename(getwd())),
                          s.base))
    })
    ##
    t1 <- strsplit(sprintf("Solved RxODE object
Dll: %s

Parameters:
     V2      V3      KA      CL       Q     Kin    Kout    EC50
 40.200 297.000   0.291  18.600  10.500   1.000   1.000 200.000


Initial Conditions:
depot centr  peri   eff
    0     0     0     1


First part of data:
# A tibble: 38 x 7
   time     depot    centr      peri      eff       C2        C3
  <dbl>     <dbl>    <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
1     0 10000.000    0.000    0.0000 1.000000  0.00000 0.0000000
2     1  7475.157 1768.532  270.6751 1.083968 43.99334 0.9113641
3     2  5587.797 2191.248  787.3677 1.179529 54.50866 2.6510696
4     3  4176.966 2076.396 1314.0348 1.227523 51.65163 4.4243597
5     4  3122.347 1783.880 1765.1486 1.233503 44.37513 5.9432612
6     5  2334.004 1465.845 2120.2772 1.214084 36.46382 7.1389804
# ... with 32 more rows", rxDll(m1)),"\n")[[1]];

    t1 <- gsub("(EC50|200.000|eff|1) *$", "\\1", t1);

    options(RxODE.display.tbl = TRUE);
    rxSyncOptions();
    test_that("print(pred); RxODE.display.tbl = TRUE", {
        tmpfile <- tmpsink();
        print(pred);
        sink();
        p1 <- readLines(tmpfile);
        unlink(tmpfile);
        p1 <- gsub("(EC50|200.000|eff|1) *$", "\\1", p1);
        expect_equal(gsub("(# A tibble: [0-9]+ ).*( [0-9]+)","\\1x\\2",p1),t1);
    })
    options(RxODE.display.tbl = FALSE)
    rxSyncOptions();
    t1 <- strsplit(sprintf("Solved RxODE object
Dll: %s

Parameters:
     V2      V3      KA      CL       Q     Kin    Kout    EC50
 40.200 297.000   0.291  18.600  10.500   1.000   1.000 200.000


Initial Conditions:
depot centr  peri   eff
    0     0     0     1


First part of data:
     time     depot    centr      peri      eff       C2        C3
[1,]    0 10000.000    0.000    0.0000 1.000000  0.00000 0.0000000
[2,]    1  7475.157 1768.532  270.6751 1.083968 43.99334 0.9113641
[3,]    2  5587.797 2191.248  787.3677 1.179529 54.50866 2.6510696
[4,]    3  4176.966 2076.396 1314.0348 1.227523 51.65163 4.4243597
[5,]    4  3122.347 1783.880 1765.1486 1.233503 44.37513 5.9432612
[6,]    5  2334.004 1465.845 2120.2772 1.214084 36.46382 7.1389804",rxDll(m1)),"\n")[[1]]
    t1 <- gsub("(EC50|200.000|eff|1) *$", "\\1", t1);


    test_that("print(pred); RxODE.display.tbl = FALSE", {
        tmpfile <- tmpsink();
        print(pred);
        sink();
        p1 <- readLines(tmpfile);
        unlink(tmpfile);
        p1 <- gsub("(EC50|200.000|eff|1) *$", "\\1", p1);
        expect_equal(p1,t1);
    })
    ################################################################################
    ## Delete dll, see what happens.
    rxDelete(m1);

    test_that("print(m1), no dll", {
        tmpfile <- tmpsink();
        print(m1);
        sink();
        p1 <- readLines(tmpfile);
        unlink(tmpfile);
        ## print(p1);
        expect_equal(p1,sprintf("RxODE model named \"%s\" (invalid object, needs to be re-created)",basename(getwd())))
    })
}, silent=TRUE)
