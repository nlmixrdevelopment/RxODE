## Refactor Jacobian calculation
##' Faster expand.grid
##'
##' Only support x and y as characters right now
##'
##' @param x first element (must be character)
##' @param y second element (must be character)
##' @param type Internal type=0L is traditional expand gridn and
##'     type=1L is jacobian expand grid (adds symbols)
##' @return Expand grid
##' @author Matthew Fidler
##' @keywords internal
##' @examples
##'
##' ##
##' rxExpandGrid(letters, letters)
##'
##' ## Another fast method; See
##' ## https://stackoverflow.com/questions/10405637/use-outer-instead-of-expand-grid
##'
##' expand.grid.jc <- function(seq1,seq2) {
##'  cbind(Var1 = rep.int(seq1, length(seq2)),
##'   Var2 = rep.int(seq2, rep.int(length(seq1),length(seq2))))
##' }
##'
##' \dontrun{
##'  microbenchmark::microbenchmark(rxExpandGrid(letters, letters), expand.grid.jc(letters, letters))
##' }
##' @export
rxExpandGrid <- function(x, y, type=0L, symengine=FALSE){
    rxExpandGrid_(x, y, type, symengine)
}

## Assumes model is loaded.
.rxJacobian <- function(model, vars=TRUE){
    .symengine <- rxIs(model, "rxS");
    if (rxIs(vars,"logical")){
        if (vars){
            .pars  <- .rxParams(model, FALSE);
            if (any(.pars=="ETA[1]")){
                .pars  <- .pars[regexpr(rex::rex(start,"ETA[",any_numbers,"]"), .pars) != -1]
            }
            .jac <- rxExpandGrid(rxState(model),
                                 c(rxState(model), .pars),
                                 1L, .symengine);
        } else  {
            .jac <- rxExpandGrid(rxState(model),
                                 rxState(model),
                                 1L, .symengine);
        }
    } else if (rxIs(vars,"character")){
        .jac <- rxExpandGrid(rxState(model),
                             c(rxState(model), vars),
                             1L, .symengine)
    }
    rxCat("Calculate Jacobian\n");
    rxProgress(dim(.jac)[1]);
    on.exit({rxProgressAbort()});
    if (.symengine){
        .ret <- apply(.jac, 1, function(x){
            .l <- x["line"];
            .l <- eval(parse(text=.l));
            rxTick();
            paste0(x["rx"], "=", rxFromSE(.l));
        })
        assign(".jacobian", .ret, envir=model);
    } else {
        .ret <- apply(.jac, 1, function(x){
            .l <- x["line"];
            .l <- rxSymPy(.l)
            rxTick();
            paste0(x["rx"], "=", rxFromSymPy(.l));
        })
    }
    rxProgressStop();
    return(.ret)
}

## Assumes .rxJacobian called on model c(state,vars)
.rxSens <- function(model, vars){
    .state <- rxState(model);
    .grd <- rxExpandSens_(.state, vars);
    rxCat("Calculate Sensitivites\n");
    rxProgress(dim(.grd)[1]);
    on.exit({rxProgressAbort()});
    rxSymPyVars(.grd$ddtS);
    rxSymPyVars(.grd$ddS2);
    .ret <- apply(.grd, 1, function(x){
        .l <- x["line"];
        .l <- rxSymPy(.l)
        .ret <- paste0(x["ddt"], "=", rxFromSymPy(.l));
        if (any(rxSymPy.vars == x["s0"])){
            .l <- x["s0D"];
            .l <- rxSymPy(.l);
            if (.l != "0"){
                .ret <- paste0(.ret, "\n", x["s0r"], "=", rxFromSymPy(.l), "+0.0");
            }
        }
        rxTick();
        return(.ret);
    })
    rxProgressStop();
    return(.ret)
}



## norm <- RxODE("
## d/dt(y)  = dy
## d/dt(dy) = mu*(1-y^2)*dy - y
## ## Initial conditions
## y(0) = 2
## dy(0) = 0
## ## mu
## mu = 1 ## nonstiff; 10 moderately stiff; 1000 stiff
## ")
