library(digest);
context("Checking output of grammars.")
df <- NULL
for (file in list.files(pattern=".*\\.test\\.g$")){
    mkdparse(file);
    out <- sprintf("%s.d_parser.c",file);
    df <- rbind(df,data.frame(basename(file),basename(out),digest(file),digest(out)));
    unlink(out);
    unlink(gsub("\\.c$",".h",out));
}
test_that("Grammars headers are produced correctly.",{
    expect_equal(digest(df),"774f2969d245acc70124ce4dec4bb712")
})
