<!--
---
output: github_document
---
-->
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- https://blog.r-hub.io/2019/12/03/readmes/ -->



# RxODE

<!-- badges: start -->
[![Build Status](https://travis-ci.org/nlmixrdevelopment/RxODE.svg?branch=master)](https://travis-ci.org/nlmixrdevelopment/RxODE)
[![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/8vv1e3hncve9tnva?svg=true)](https://ci.appveyor.com/project/mattfidler/rxode)
[![codecov.io](https://codecov.io/github/nlmixrdevelopment/RxODE/coverage.svg)](https://codecov.io/github/nlmixrdevelopment/RxODE?branch=master)
[![CRAN version](http://www.r-pkg.org/badges/version/RxODE)](https://cran.r-project.org/package=RxODE)
[![CRAN total downloads](https://cranlogs.r-pkg.org/badges/grand-total/RxODE)](https://cran.r-project.org/package=RxODE)
[![CRAN total downloads](https://cranlogs.r-pkg.org/badges/RxODE)](https://cran.r-project.org/package=RxODE)

<!-- badges: end -->

## Overview

**RxODE** is an R package for solving and simulating from ode-based
models. These models are convert the RxODE mini-language to C and
create a compiled dll for fast solving. ODE solving using RxODE has a
few key parts:

 - `RxODE()` which creates the C code for fast ODE solving based on a
   [simple syntax](https://nlmixrdevelopment.github.io/RxODE/articles/RxODE-syntax.html) related to Leibnitz notation.
 - The event data, which can be:
   - a `NONMEM` or `deSolve` [compatible data frame](https://nlmixrdevelopment.github.io/RxODE/articles/RxODE-event-types.html), or
   - created with `et()` or `EventTable()` for [easy simulation of events](https://nlmixrdevelopment.github.io/RxODE/articles/RxODE-event-table.html)
   - The data frame can be augmented by adding
     [time-varying](https://nlmixrdevelopment.github.io/RxODE/articles/RxODE-covariates.html#time-varying-covariates)
     or adding [individual covariates](file:///home/matt/src/RxODE/docs/articles/RxODE-covariates.html#individual-covariates) (`iCov=` as needed)
 - `rxSolve()` which solves the system of equations using initial
   conditions and parameters to make predictions
   - With multiple subject data, [this may be
     parallelized](https://nlmixrdevelopment.github.io/RxODE/articles/RxODE-speed.html).
   - With single subject the [output data frame is adaptive](https://nlmixrdevelopment.github.io/RxODE/articles/RxODE-data-frame.html)
   - Covariances and other metrics of uncertanty can be used to
     [simulate while solving](https://nlmixrdevelopment.github.io/RxODE/articles/RxODE-sim-var.html)

## Installation

You can install the released version of RxODE from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("RxODE")
```

To run RxODE, you need a working c compiler.  To use parallel threaded
solving in RxODE, this c compiler needs to support open-mp.

You can check to see if R has working c compiler you can check with:

```r
## install.packages("pkgbuild")
pkgbuild::has_build_tools(debug = TRUE)
```

If you do not have the toolchain, you can set it up as described by
the platform information below:

### Windows

In windows you may simply use installr to install rtools:

```r
install.packages("installr")
library(installr)
install.rtools()
```

Alternatively you can
[download](https://cran.r-project.org/bin/windows/Rtools/) and install
rtools directly.

### Mac OSX

Installation on a mac is much similar to RxODE installation under
windows.  To enable open mp on R and RxODE, you will need to install
the gfortran and clang compilers located at
https://cran.r-project.org/bin/macosx/tools/

### Linux

To install on linux make sure you install gcc (with openmp support)
and gfortran using your distribution's package manager.


## Development Version

Since the development version of RxODE uses StanHeaders, you will need
to make sure your compiler is setup to support C++14, as described in
the [rstan setup page](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started#configuration-of-the-c-toolchain)

Once the C++ toolchain is setup appropriately, you can install the
development version from
[GitHub](https://github.com/nlmixrdevelopment/RxODE) with:

``` r
# install.packages("devtools")
devtools::install_github("nlmixrdevelopment/RxODE")
```

# Related R Packages

This is a brief comparison of pharmacometric ODE solving R packages to
`RxODE`.

There are several [R packages for differential
equations](https://cran.r-project.org/web/views/DifferentialEquations.html).
The most popular is
[deSolve](https://cran.r-project.org/package=deSolve).  

However for pharmacometrics-specific ODE solving, there are only 2
packages other than [RxODE](https://CRAN.R-project.org/package=RxODE)
released on CRAN.  Each uses compiled code to have faster ODE solving.

- [mrgsolve](https://CRAN.R-project.org/package=mrgsolve), which uses
  a thread-safe C++ lsoda solver to solve ODE systems.  The user is
  required to write hybrid R/C++ code to create a mrgsolve
  model.
  
  In contrast, `RxODE` has a R-like mini-language that is parsed into
  C code that solves the ODE system.
  
  Both `mrgsolve` and `RxODE` implement thread-safe ODE solving. In
  `RxODE` threading for each individual ODE solved is turned on by
  default and the supported RxODE functions are all thread-safe.
  
  Unlike `RxODE`, `mrgsolve` does not currently support symbolic
  manipulation of ODE systems, like automatic Jacobian calculation or
  forward sensitivity calculation (`RxODE` currently supports this and
  this is the basis of
  [nlmixr](https://cran.r-project.org/package=nlmixr)'s FOCEi
  algorithm)

- [dMod](https://cran.r-project.org/package=dMod), which uses a unique
  syntax to create "reactions".  These reactions are not low-level ODE
  equations, but create them and then created c code for a compiled
  deSolve model.
  
  In contrast `RxODE` defines ODE systems at a lower level.  `RxODE`'s
  parsing of the mini-language comes from C, whereas `dMod`'s parsing
  comes from R.
  
  Like `RxODE`, `dMod` supports symbolic manipulation of ODE systems
  and calculates forward sensitivities and adjoint sensitivities of
  systems.
  
  Unlike `mrgsolve` and `RxODE`, `dMod` is not thread-safe since
  `deSolve` is not yet thread-safe.

And there is one package that is not released on CRAN:

- [PKPDsim](https://github.com/InsightRX/PKPDsim) which defines models
  in an R-like syntax and converts the system to compiled code. 
  
  Like `mrgsolve`, `PKPDsim` does not currently support symbolic
  manipulation of ODE systems.
  
  `PKPDsim` is not thread-safe.

The open pharmacometrics open source community is fairly friendly, and
the RxODE maintainers has had positive interactions with all the
maintainer all of the ODE-solving pharmacometric projects listed.
