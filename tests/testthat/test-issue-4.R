if (Sys.getenv("covr") !="true"){
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
        fun <- function(){
            tmpDir <- tempfile();
            if (dir.exists(tmpDir))
                unlink(tmpDir, recursive=TRUE)
            dir.create(tmpDir, recursive=TRUE)
            owd <- getwd();
            on.exit({setwd(owd);
                if (dir.exists(tmpDir))
                    unlink(tmpDir, recursive=TRUE)});
            setwd(tmpDir)
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
                devtools.name="Matt Fidler"
            )

            ## Create a package to do testing on...
            ## devtools::create("tmp");
            usethis::create_package(tmpDir);
            usethis::use_testthat();
            test <- file.path(tmpDir, "tests/testthat.R")
            test1 <- file.path(tmpDir, "tests/testthat/test-inside-package.R");
            l <- writeLines(c("Sys.setenv(R_TESTS = \"\")", readLines(test)),
                            test);
            writeLines(c("test_that(\"inside_package rxode check\", {", "    model <- RxODE::RxODE({",  "        C2 = centr/V2", "        C3 = peri/V3", "        d/dt(depot) = -KA * depot",  "        d/dt(centr) = KA * depot - CL * C2 - Q * C2 + Q * C3",  "        d/dt(peri) = Q * C2 - Q * C3", "        d/dt(eff) = Kin - Kout * (1 - C2/(EC50 + C2)) * eff",  "    }, modName = paste(sample(LETTERS, 5), collapse = \"\"), wd = tempfile())",  "    expect_equal(RxODE::rxState(model), c(\"depot\", \"centr\", \"peri\", ",  "        \"eff\"))", "    expect_equal(RxODE::rxParams(model), c(\"V2\", \"V3\", \"KA\", ",  "        \"CL\", \"Q\", \"Kin\", \"Kout\", \"EC50\"))", "})"),
                       test1)
            tmp <- devtools::check(tmpDir, document=FALSE, quiet=TRUE);
            sink()
            unlink(tf);

            test_that("Inside a different package works.", {
                expect_equal(length(tmp$errors), 0)
            })

        }
        fun();
    }, silent=TRUE, cran=FALSE)
}
