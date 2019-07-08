rxPermissive({

    expand.grid.jc <- function(seq1,seq2) {
      cbind(Var1 = rep.int(seq1, length(seq2)),
       Var2 = rep.int(seq2, rep.int(length(seq1),length(seq2))))
     }

    x <- microbenchmark::microbenchmark(rxExpandGrid(letters, letters), expand.grid.jc(letters, letters))

    library(dplyr)

    x <- as.data.frame(x) %>%
        group_by(expr) %>%
        summarize(x=median(time)) %>%
        arrange(expr)

    context("Faster Expand Grid")
    test_that("rxExpandGrid is faster that expand.grid.jc", {
        expect_true(x$x[1] < x$x[2]);
    })

}, on.validate=TRUE)
