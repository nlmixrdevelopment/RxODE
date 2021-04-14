# RxODE development version

# RxODE 1.0.8

* Fix issue #399


# RxODE 1.0.7 

* Change syntax vignette to use markdown option
  `screenshot.force=FALSE`.  This should get rid of the `webshot`
  error

* Change to depend on dparser 1.3.0, which has some memory fixes

# RxODE 1.0.6

* RxODE imports but does not link to `checkmate` any longer.  This change
   should make recompilation of RxODE to work with different releases
   of `checkmate` unnecessary.

* Default Solaris solver changed back to "lsoda"

* Fix Bug #393, where in certain circumstances `rxSolve(...,theta=)`
  did not solve for all subjects.

* Will not ignore NEWS and README when building the package so that
  they will show up on CRAN.  You can also access the news by
  `news(package="RxODE")`

* Changed `ODR` model names from time id to `_rx` followed by the
  `md5` hash id and a per-session counter id; For packages the id is
  `_rxp` followed by the `md5` hash and a per-session counter id.

* Changed `qs` to be more conservative in hash creation. Add a check
  hash as well as NOT using altrep stringfish representation.

# RxODE 1.0.5

* Maintainance release -- use `std::floor` and cast variables to
  `double` for internal C functions.  This should allow a successful
  compile on Solaris CRAN.

* Changed `units` from an Imports to a Suggests to allow testing on
  Solaris rhub

* Changed `ODR` model names from time id to `_rx` followed by the
  `md5` hash id; For packages the id is `_rxp` followed by the `md5`
  hash.

* Removed AD linear compartment solutions for Windows R 3.6, though
  they still work for Windows R 4.0 (You can get them back for Windows
  R 3.6 if you install `BH` 1.66.0-1 and then recompile from source).

   - This will cause `nlmixr` to fail with solved systems on Windows 3.6.
     Currently the Stan Headers do not compile on this system so they are
     disabled at this time.

 * RxODE imports but does not link to `qs` any longer; This change
   should make recompilation of RxODE to work with different releases
   of `qs` unnecessary.

 * RxODE now checks for binary compatibility for `Rcpp`, `dparser`,
   `checkmate`, and `PreciseSums`

# RxODE 1.0.4
## Breaking changes

* RxODE can only use supported functions (could be breaking); You may
  add your own functions with `rxFun` and their derivatives with `rxD`

* RxODE now uses its own internal truncated multivariate normal
  simulations based on the threefry sitmo library.  Therefore random
  numbers generated within `RxODE` like providing
  `rxSolve(...,omega=)` will have different results with this new
  random number generator.  This was done to allow internal re-sampling
  of sigmas/etas with thread-safe random number generators (calling R
  through `mvnfast` or R's simulation engines are not thread safe).

* `RxODE` now moved the precise sum/product type options for `sum()`
  and `prod()` to `rxSolve` or `rxControl`

* `cvPost` now will returned a named list of matrices if the input
  matrix was named

* `rxSolve` will now return an integer `id` instead of a factor `id`
  when `id` is integer or integerish (as defined by checkmate).
  Otherwise a factor will be returned.

* When mixing ODEs and `linCmt()` models, the `linCmt()` compartments
  are 1 and possibly 2 instead of right after the last internal ODE.
  This is more aligned with how PK/PD models are typically defined.

* `EVID=3` and `EVID=4` now (possibly) reset time as well.  This
  occurs when the input dataset is sorted before solving.

* When `EVID=2` is present, an `evid` column is output to distinguish
  `evid=0` and `evid=2`

## New features

* Add the ability to order input parameters with the `param()`
  pseudo-function

* Add the ability to resample covariates with `resample=TRUE` or
  `resample=c("SEX", "CRCL")`.  You can resample all the covariates by
  `ID` with `resampleID=TRUE` or resample the covariates without
  respect to `ID` with `resampleID=FALSE`

* Comparison of factors/strings is now supported in `RxODE`; Therefore
  ID=="Study-1" is now allowed.

* Completion for elements of `rxSolve()` objects, and `et()`
  objects have been added (accessed through `$`)

* Completion of `rxSolve()` arguments are now included since they are
  part of the main method

* Allow simulation with zero matrices, that provide the simulation
  without variability.  This affects `rxSolve` as well as `rxMvnrnd` and
  `cvPost` (which will give a zero matrix whenever one is specified)

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

* Allow `while(logical)` statements with ability to break out if them
  by `break`. The while has an escape valve controlled by `maxwhere`
  which by default is 10000 iterations. It can be change with
  `rxSolve(..., maxwhere = NNN)`

* Allow accessing different time-varying components of an input
  dataset for each individual with:

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

  This allows NONMEM-style of calculating parameters like tad:

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

* Added `tad()`, `tad(cmt)` functions for time since last dose and time
  since last dose for a compartment; Also added time after first dose
  and time after first dose for a compartment `tafd()`, `tafd(cmt)`;
  time of last dose `tlast()`, `tlast(cmt)` and dose number
  `dosenum()` (currently not for each compartment)

* Changed linear solved systems to use "advan" style `linCmt()`
  solutions, to allow correct solutions of time-varying covariates
  values with solved systems; As such, the solutions may be slightly
  different.  Infusions to the depot compartment are now supported.


* Added sensitivity auto-differentiation of `linCmt()` solutions.
  This allows sensitivities of `linCmt()` solutions and enables
  `nlmixr` focei to support solved systems.
  - One solution is to use Stan's auto-differentiation which requires
    `C++14`

* When calculating the empirical Bayesian estimates for with `rxInner`
  (used for nlmixr's 'focei') ignore any variable beginning with `rx_`
  and `nlmixr_` to hide internal variables from table output.  This
  also added `tad=tad()` and `dosenum=dosenum()` to the `ebe` output
  allowing grouping by id, dose number and use TAD for individual plot
  stratification.

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

* Allowed the following data variables can be accessed (but not
  assigned or used as a state):
  - `cmt`
  - `dvid`
  - `addl`
  - `ss`
  - `amt`
  - `rate`
  - `id` which requires calling the id as factor `ID=="1"` for
    instance.

* Kept `evid` and `ii` as restricted items since they are not part of
  the covariate table and are restricted in use.

* Added the following random number generators; They are thread safe
  (based on `threefry` `sitmo` and c++11) and your simulations with
  them will depend on the number of cores used in your simulation (Be
  careful about reproducibility with large number of threads; Also
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

  In addition, while initializing the system, the following values are
  simulated and retained for each individual:

  - `rinorm()` and `rinormV()` (low discrepancy normal)
  - `ricauchy()`
  - `richisq()`
  - `riexp()`
  - `rif()`
  - `rigamma()`
  - `ribeta()`
  - `rigeom()`
  - `ripois()`
  - `rit()`
  - `riunif()`
  - `riweibull()`

* Added `simeta()` which simulates a new `eta` when called based
  on the possibly truncated normal `omega` specified by the original
  simulation.  This simulation occurs at the same time as the ODE is
  initialized or when an ODE is missing, before calculating the final
  output values.  The `omega` will reflect whatever study is being simulated.

*  Added `simeps()` which simulates a new `eps` from the possibly
  truncated normal `sigma` at the same time as calculating the final
  output values. Before this time, the `sigma` variables are zero.

  All these change the solving to single thread by default to make sure the
  simulation is reproducible. With high loads/difficult problems the
  random number generator may be on a different thread and give a
  different number than another computer/try.

  Also please note that the `clang` and `gcc` compiler use different
  methods to create the more complex random numbers.  Therefore
  `MacOS` random numbers will be different than `Linux`/`Windows` at
  this time (with the exception of uniform numbers).

  These numbers are still non-correlated random numbers (based on the
  sitmo test) with the exception of the vandercorput distributions, so
  if you increase the number of threads (cores=...) the results should
  still be valid, though maybe harder to reproduce.  The faster the
  random number generation, the more likely these results will be
  reproduced across platforms.

* Added the ability to integrate standard deviations/errors of omega
  diagonals and sigma diagonals.  This is done by specifying the omega
  diagonals in the theta matrix and having them represent the
  variabilities or standard deviations. Then these standard deviations
  are simulated along with the correlations using the IJK correlation
  matrix (omega dimension < 10) or a correlation matrix or Inverse
  Wishart-based correlation matrix (omega dimension > 10).  The
  information about how to simulate this is in the variability
  simulation vignette.

* Now have a method to use `lotri` to simulate between occasion
  variability and other levels of nesting.

* Added lower gamma functions See Issue #185

* Upgraded comparison sort to timsort 2.0.1

* Changed in-place sort to a modified radix sort from
  `data.table`.  The radix search was modified to:
 - Work directly with `RxODE` internal solved structures
 - Assume no infinite values or `NA`/`NaN` values of time
 - Always sort time in ascending order
 - Changed sorting to run in a single thread instead of taking over
   all the threads like data.table

* Changed method for setting/getting number of threads based on
  `data.table`'s method

* Added function `rxDerived` which will calculate derived parameters
  for 1, 2, and 3 compartment models

* More descriptive errors when types of input are different than expected

## Engine changes

* Moved many C functions to C++.  CRAN OpenMP support requires C++
  only when C and C++ are mixed.  See:

  https://stackoverflow.com/questions/54056594/cran-acceptable-way-of-linking-to-openmp-some-c-code-called-from-rcpp

* No longer produces C code that create the model variables. Instead,
  use `qs` to serialize, compress and encode in base91 and then write
  the string into the C file. The `qs` package then decodes all of
  that into the model variables.  This also increases the compilation
  speed for models in RxODE.

* Pre-compile RxODE headers once (if cache is enabled), which
  increases compilation speed for models in RxODE

* `RxODE`'s translation from the mini-language to C has been refactored

## Bug fixes:
 - Occasionally RxODE misidentified dual `lhs`/`param` values.  An
   additional check is performed so that this does not happen.

 - For solved matrices with similar names (like "tadd" and "tad")
   RxODE will now prefer exact matches instead of the first match
   found when accessing the items with `$tad`.

 - A fix where all ID information is kept with `keep=c(""..."")`

 - Transit compartment models using the `transit` ODE or variable are
   now allowed.  Also check for more internally parsed items (see
   Issue #145).

 - Bug fix for `etSeq` and `etRep` where greater than 2 items were
   mis-calculated

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

<!--  LocalWords:  resample covariates Rstudio NONMEM advan focei
 -->
<!--  LocalWords:  nlmixr's symengine linearization RxODE
 -->
<!--  LocalWords:  reproducibility
 -->
