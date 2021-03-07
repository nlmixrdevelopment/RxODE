# CRAN Comments

The following platforms were tested with a succesful check: MacOS release, oldrel, Windows
  release, oldrel, ubuntu devel, ubuntu release and ubuntu oldrel:

https://github.com/nlmixrdevelopment/RxODE/actions/runs/628714043

Additionally solaris was tested on our local virtual machines as well
as on rhub with a successful check.


* Solaris issue [here](https://www.r-project.org/nosvn/R.check/r-patched-solaris-x86/RxODE-00install.html);

  By removing the suggests ‘devtools’ ‘symengine’ ‘units’ ‘usethis’,
  the rhub check ran with minimal errors (related to removing the
  packages from the build). The removed packages currently do not
  successfully build with Solaris rhub.

  - https://builder.r-hub.io/status/RxODE_1.0.5.tar.gz-94aedb1ff28849e8b651c6bfff3da783

* Old Release Mac Issue [here](https://www.r-project.org/nosvn/R.check/r-oldrel-macos-x86_64/RxODE-00check.html)

  We believe this was caused by using unique symbols in ODE models
  based on time-stamps.  For fast compiling/fast parsing models the
  time stamps could be the same.  This caused a ODR-violation crash.

  The unique symbols were changed to model-based md5 hashes which
  should be unique for each model that is produced.  If the model has
  the same md5 instead of recompiling the model, it is reloaded into
  the R session.

  The check for Old Release Mac on our systems are clean:

  https://github.com/nlmixrdevelopment/RxODE/runs/1983444522?check_suite_focus=true

  We also verified the ASAN version on our machines

* New Release Mac Issue [here](https://www.r-project.org/nosvn/R.check/r-release-macos-x86_64/RxODE-00check.html)

  This check cannot find the csv files in the vignette file; We moved
  this to a pkgdown-only

* Old Release Windows Issue [here](https://www.r-project.org/nosvn/R.check/r-oldrel-windows-ix86+x86_64/RxODE-00check.html)

  This is caused by `BH` being updated to support C++14, which isn't
  supported by the old 3.5 Rtools.  We turned of stan AD support for
  Windows R 3.6 unless they are using `BH` `1.66.0-1` or earlier.  If
  they try to use some of the stan-based functionality the user is
  suggested to upgrade R and Rtools or temporarily downgrade `BH`.

  With this new change, the R build for the old release works well:
