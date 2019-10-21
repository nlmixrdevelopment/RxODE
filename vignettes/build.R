library(RxODE)

.files <- list.files(pattern=".*[.]Rmd")
.files <- .files[!(.files %in% c("RxODE-intro.Rmd", "RxODE-syntax.Rmd"))]
for (.f in .files){rmarkdown::render(.f)}
