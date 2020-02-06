##' Expand the omega matrix to include the IOV matrix
##'
##' .. content for \details{} ..
##' @param omega Original (named) omega matrix
##' @param iov Original (named) iov matrix
##' @param iovNames IOV names that are used in the model
##' @param omegaMatNames Omega matrix names (if not in omega matrix)
##' @param iovMatNames IOV matrix names (if not in iov matrix)
##' @return List with three elements; `m` is Named matrix for internal
##'     R simulation. `lower`/`upper` are the upper and lower
##'     boundaries of the truncated multi-variate normal.
##' @author Matthew Fidler
##' @examples
##'
##' omega <- lotri(eta.cl+eta.v ~c(1, 0.5, 1))
##' iov <- lotri(iov.cl+iov.v ~c(1, 0.5, 1))
##' .rxExpandOmegaIov(omega, iov, paste0("ETA[",1:4, "]"))
##'
##' .rxExpandOmegaIov(omega)
##' @noRd
.rxExpandOmegaIov <- function(omega, iov=NULL, iovNames=NULL,
                              omegaLower=-Inf, omegaUpper=Inf,
                              iovLower=-Inf, iovUpper=Inf,
                              omegaMatNames=NULL, iovMatNames=NULL){
    if (!is.null(omegaMatNames)) {
        .dim1 <- omegaMatNames
    } else {
        if (is.null(dimnames(omega))) stop("'omega' needs to be named")
        .dim1 <- dimnames(omega)[[1]]
    }
    .len0 <- length(.dim1)
    if (length(omegaLower) == 1) omegaLower <- rep(omegaLower, .len0)
    else if (length(omegaLower) != .len0)
        stop("'omega' and 'omegaLower' need compatible dimensions")
    if (length(omegaUpper) == 1) omegaUpper <- rep(omegaUpper, .len0)
    else if (length(omegaUpper) != .len0)
        stop("'omega' and 'omegaUpper' need compatible dimensions")
    if (is.null(iov) && is.null(iovNames)){
        return(list(m=omega, lower=omegaLower,
                    upper=omegaUpper))
    } else {
        if (!is.null(iovMatNames)) {
            .dim2 <- iovMatNames
        } else {
            if (is.null(dimnames(iov))) stop("'iov' needs to be named")
            .dim2 <- dimnames(iov)[[1]]
        }
        if (is.null(iovNames)) stop("'iovNames' needs to have model-based terms")
        .len1 <- length(.dim2)
        if (length(iovLower) == 1) iovLower <- rep(iovLower, .len1)
        else if (length(iovLower) != .len0)
            stop("'iov' and 'iovLower' need compatible dimensions")
        if (length(iovUpper) == 1) iovUpper <- rep(iovUpper, .len1)
        else if (length(iovUpper) != .len0)
            stop("'iov' and 'iovUpper' need compatible dimensions")
        .len2 <- length(iovNames)
        if (.len2 %% .len1 != 0)
            stop('iov and full iov expansion does not match');
        .len3 <- .len2 / .len1
        iovLower <- rep(iovLower, .len3)
        iovUpper <- rep(iovUpper, .len3)
        .m <- lapply(seq(1, .len3), function(...) iov)
        .m <- as.matrix(Matrix::bdiag(c(list(omega), .m)));
        .dn <- c(.dim1, iovNames)
        dimnames(.m) <- list(.dn, .dn)
        return(list(m=.m,
                    lower=c(omegaLower, iovLower),
                    upper=c(omegaUpper, iovUpper)))
    }
}
