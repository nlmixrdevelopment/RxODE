rxPermissive({
    context("Working Directory + Model Name (Issue #4)")
    model_text <- "
  C2 = centr/V2;
  C3 = peri/V3;
  d/dt(depot) =-KA*depot;
  d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
  d/dt(peri)  =                    Q*C2 - Q*C3;
  d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
  "
    modName <- paste(sample(LETTERS, 5), collapse="")
    model   <- RxODE::RxODE(model = model_text, modName=modName, wd='mytempdir')

    test_that("Model has the correct states", {
        expect_equal(rxState(model),
                     c("depot", "centr", "peri", "eff"));
        expect_equal(rxParams(model),
                     c("V2", "V3", "KA", "CL", "Q", "Kin", "Kout", "EC50"));
    })

    rxDelete(model);
    unlink("mytempdir", recursive=TRUE);

}, silent=TRUE, cran=TRUE)

rxPermissive({
    context("Working Directory + Model Name inside package (Issue #4)")
    tf <- tempfile();
    sink(tf)
    authors_at_r <- paste0(
        "'",
        person(
            "Matthew",
            "Fidler",
            email = "matthew.fidler@gmail.com",
            role  = c("aut", "cre")),
        "'"
    )
    options(devtools.desc.author=authors_at_r)
    options(
        repos = c(CRAN = "https://cran.rstudio.com/"),
        devtools.name="Matt Fidler",
        devtools.desc.license="GPL (>= 2)"#,
    )

    ## Create a package to do testing on...
    if (dir.exists("tmp"))
        unlink("tmp", recursive=TRUE)
    devtools::create("tmp");
    devtools::use_testthat("tmp");
    l <- writeLines(c("Sys.setenv(R_TESTS = \"\")", readLines("tmp/tests/testthat.R")),
                    "tmp/tests/testthat.R");
    writeLines(c("test_that(\"inside_package rxode check\", {", "    model <- RxODE::RxODE({",  "        C2 = centr/V2", "        C3 = peri/V3", "        d/dt(depot) = -KA * depot",  "        d/dt(centr) = KA * depot - CL * C2 - Q * C2 + Q * C3",  "        d/dt(peri) = Q * C2 - Q * C3", "        d/dt(eff) = Kin - Kout * (1 - C2/(EC50 + C2)) * eff",  "    }, modName = paste(sample(LETTERS, 5), collapse = \"\"), wd = tempfile())",  "    expect_equal(RxODE::rxState(model), c(\"depot\", \"centr\", \"peri\", ",  "        \"eff\"))", "    expect_equal(RxODE::rxParams(model), c(\"V2\", \"V3\", \"KA\", ",  "        \"CL\", \"Q\", \"Kin\", \"Kout\", \"EC50\"))", "})"),
               "tmp/tests/testthat/test-inside-package.R")
    tmp <- devtools::check("tmp", document=FALSE, quiet=TRUE);
    sink()
    unlink(tf);

    test_that("Inside a different package works.", {
        expect_equal(length(tmp$errors), 0)
    })
    unlink("tmp", recursive=TRUE)
}, silent=TRUE, cran=FALSE)
