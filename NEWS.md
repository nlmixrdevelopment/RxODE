# RxODE xxxx
## Breaking changes

* `lag(cmt) = ` is no longer supported because `lag()` is now allowed
  to use lagged covariates and observations.

* RxODE can only use supported functions (could be breaking); You may
  add your own functions with `rxFun` and their derivatives with `rxD`

* RxODE now uses its own internal truncated multivariate normal
  simulations based on the threefry sitmo library.  Therefore random
  numbers generated within `RxODE` like providing
  `rxSolve(...,omega=)` will have different results with this new
  random number generator.  This was done to allow internal resampling
  of sigmas/etas with thread-safe random number generators (calling R
  through `mvnfast` or R's simulation engines are not thread safe).
  
  * A backward-compatible simulation can be used by `rxSolve(...,
    mvnfast=TRUE)` (but of course requires the `mvnfast` package
    installed);  NEED TO CHECK IF THIS IS STILL TRUE;  Is this worthwhile?

* `RxODE` now moved the precise sum/product type options for `sum()`
  and `prod()` to `rxSolve` o `rxControl`
  
* `cvPost` now will returned a named list of matrices if the input matrix was named

## New features

* Completion for all api elements of `rxSolve()` objects, and `et()`
  objects have been added (accessed through `$`)

* `et()` can dose with `length(amt) > 1` as long as the other
  arguments can create a event table.

* Rstudio notebook output makes more sense

* Printing upgraded to cli 2.0
  
* Caching of internal C data setup is now supported increasing speed
  of `optim` code when:
  - Event Table doesn't change
  - The size of the parameters doesn't change
  - `inits` do not change (though you can specify them as `cmt(0)=...`
    in the model and change them by parameters)
  - See Issue #109

* Allow accessing different time-varying components of an input dataset for each indivdiual with:
  - `lag(var, #)`
  - `lead(var, #)`
  - `first(var)`
  - `last(var)`
  - `diff(var)`
Each of these are similar to the R `lag`, `lead`, `first`, `last` and
`diff`.  However when undefined, it returns `NA`

* Allow sticky left-handed side of the equation; This means for an
  observation the left handed values are saved for the next
  observations and then reassigned to the last calculated value.
  
  This allows nonmem-style of calculating parameters like tad:
  
```r
mod1 <-RxODE({
    KA=2.94E-01;
    CL=1.86E+01;
    V2=4.02E+01;
    Q=1.05E+01;
    V3=2.97E+02;
    Kin=1;
    Kout=1;
    EC50=200;
    C2 = centr/V2;
    C3 = peri/V3;
    d/dt(depot) =-KA*depot;
    d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
    d/dt(peri)  =                    Q*C2 - Q*C3;
    d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
    if (!is.na(amt)){
        tdose <- time
    } else {
        tad <- time - tdose
    }
})
```

It is still simpler to use:

```r
mod1 <-RxODE({
    KA=2.94E-01;
    CL=1.86E+01;
    V2=4.02E+01;
    Q=1.05E+01;
    V3=2.97E+02;
    Kin=1;
    Kout=1;
    EC50=200;
    C2 = centr/V2;
    C3 = peri/V3;
    d/dt(depot) =-KA*depot;
    d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
    d/dt(peri)  =                    Q*C2 - Q*C3;
    d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
    tad <- time - tlast
})
```
If the `lhs` parameters haven't been defined yet, they are `NA`

* Now the NONMEM-style `newind` flag can be used to initialize `lhs`
  parameters.

* Added "advan" style `linCmt()` solutions, to allow correct solutions
  of time-varying covariate values with solved systems

* Added sensitivity auto-differentiation of `linCmt()` solutions (via
  stan's math headers).  This allows sensitivities of `linCmt()`
  solutions and enables `nlmixr` focei to support solved systems.
  - As such, `RxODE` now requires `C++14` support.
  
* Added ability to prune branching with `rxPrune`. This converts
  `if`/`else` or `ifelse` to single line statements without any
  `if`/`then` branching within them.
  
* Added ability to take more complex conditional expressions, including:
  - `ifelse(expr, yes, no)`
  - `x = (x==1)*1 + (!(x==1))*2`
  - `if (logic){ expr} else if (logic) {expr} else {}`.  The preferred
    syntax is still only `if`/`else` and the corresponding parsed code
    reflects this preference.
    - Note `ifelse` is not allowed as an ODE compartment or a variable.
	
* Switched to `symengine` instead of using `sympy`
  - Remove dependence on python.
  - Since symengine is C-based and doesn't require the python
    interface it is much faster than `sympy`, though some functions in
    `sympy` are no longer accessible.
  - Also symengine requires R 3.6, so now RxODE requires R 3.6
	
* Added new ODE solving method "indLin", or inductive linearization.
  When the full model is a linear ODE system this becomes simply the
  matrix exponential solution.  Currently this requires a different
  setup.

* Added arbitrary function definition to RxODE using `rxFun`
  - Requires function, arguments and corresponding C-code
  - Derivatives (if required) can be added to the derivative table
    `rxD`.  When taking deviates without a derivative function, RxODE
    will use numerical differences.

* Will error if RxODE does not know of a function that you are trying
  to use; This could be a breaking change.  Currently:
  - C's functions from `math.h` are supported
  - R's function returning and taking doubles are supported
  - Other functions can be added using `rxFun` and `rxD`

* Added `NA`, `NaN`, `Inf` and `+Inf` handling to a RxODE model.  Can
  be useful to diagnose problems in models and provide alternate
  solutions. In addition, added R-like functions `is.nan`, `is.na`,
  `is.finite` and `is.infinite` which can be called within the RxODE
  block.

* Allowed the following data variables to be accessed (but not
  assigned or used as a state):
  - `cmt`
  - `dvid`
  - `addl`
  - `ss`
  - `amt`
  - `rate`

* Kept `evid` and `ii` as restricted items since they are not part of
  the covariate table and are restricted in use.
  
* Added the following random number generators; They are thread safe
  (based on `threefry` `sitmo` and c++11) and your simulations with
  them will depend on the number of cores used in your simulation (Be
  careful about reproduciblility with large number of threads; Also
  use parallel-solve type of RxODE simulations to avoid the [birthday
  problem](https://www.johndcook.com/blog/2016/01/29/random-number-generator-seed-mistakes/)).
  
  
  During ODE solving, the values of these are `0`, but while
  calculating the final output the variable is randomized at least for
  every output. These are:
  - `rxnorm()` and `rxnormV()` (low discrepancy normal)
  - `rxcauchy()`
  - `rxchisq()`
  - `rxexp()`
  - `rxf()`
  - `rxgamma()`
  - `rxbeta()`
  - `rxgeom()`
  - `rxpois()`
  - `rxt()`
  - `rxunif()`
  - `rxweibull()`
  
* Now have a method to use `lotri` to simulate between occasion
  variability and other levels of nesting.
  
* Updated to [TimSort 1.2.0](https://github.com/timsort/cpp-TimSort/releases/tag/v1.2.0)

* More descriptive errors when types of input are different than expected
  
## Bug fixes:
 - Occasionally RxODE misidentified dual `lhs`/`param` values.  An
   additional check is performed so that this does not happen.
 - Transit compartment models using the `transit` ODE or variable are
   now allowed.  Also check for more internally parsed items (see
   Issue #145).

# RxODE v0.9.2-0
* New plotting engine
* Various bug fixes for upcoming R 4.0 release:
  - Dropped some imports for 21 imports restriction
  - Fixed incompatibility with new `ggplot2` 3.3.0
  - Fixed allowing `NA`s in RxODE dataset
  - Fixed setting all compartment default values for bioavailability, rate, etc.
  - Added additional protection against floating point -> NaN for power functions
  
# RxODE v0.9.1-9
* Minor namespace/documentation changes for R 4.0 compatibility

# RxODE v0.9.1-8
* Added the ability to have an input parameter to be assigned to a new
  value (Issue #135)
* Added LINPACK authors as contributors
* Added a `NEWS.md` file to track changes to the package
