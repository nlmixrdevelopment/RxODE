# CRAN Comments

* Please add \value to .Rd files regarding exported methods and
  explain the functions results in the documentation. Please write
  about the structure of the output (class) and also what the output
  means. (If a function does not return a value, please document that
  too, e.g.  \value{No return value, called for side effects} or
  similar)

```
Missing Rd-tags:
      as.data.table.rxEt.Rd: \value
      ...
```

 - Completed

* Please do not modify the .GlobalEnv. This is not allowed by the CRAN
  policies.

- The global environment is not modified by the RxODE functions,
  though there is one test that modifies the global environment by
  using `assign()`.  The authors believed that since the other tests
  modified the global environment while running, this was acceptable
  and within CRAN policies.  That test has been removed for this CRAN
  package submission.  If you think this argument makes sense we can
  add it back for the next release.
  
- Other references to `globalenv()` in `RxODE` are used in
  printouts to show immediately where to access certain properties
  without consulting the documentation.  Currently this searches the
  parent environment until it reaches the global environment (or there
  are no more environments to search).  For example the print-out
  would say `solveDf$params` instead of `$params` which may be confusing
  to a novice R user. While this keyword is used in the RxODE package
  it is not used to change the global environment.

- The last reference to the `globalenv()` calls `eval(parse(text =
  sprintf(".Call(\"%s\")", .modVars)), envir = .GlobalEnv)`.  This is
  simply to get the proper environmental `.Call` for the user-built
  `dll` and does not modify the environment.  Rather it allows all the
  platforms tested to grab the ODE model variables from the DLL
  correctly.  Without this evaluation certain platforms do not get the
  model variables and will error on a simple `.Call`.

- In the RxODE github repository there are building scripts that
  modify the global environment which are not included in the package
  submitted to CRAN.

## Test environments
* local R installation, R 4.0.3
* ubuntu 16.04 (on travis-ci), R 4.0.3
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
