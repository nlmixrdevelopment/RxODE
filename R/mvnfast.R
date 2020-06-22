.mvnfast <- NULL
.rmvn <- function(n, mu, sigma, ncores = 1, isChol = FALSE, A = NULL,
                  kpnames = FALSE) {
  if (is.null(.mvnfast)) {
    rxReq("mvnfast")
  }
  .mvnfast$rmvn(n, mu, sigma, ncores, isChol, A, kpnames)
}

.rmvt <- function(...) {
  if (is.null(.mvnfast)) {
    rxReq("mvnfast")
  }
  .mvnfast$rmvn(...)
}
