# RxODE v0.9.2-0
* New plotting engine
* Various bug fixes for upcoming R 4.0 release:
  - Dropped some imports for 21 imports restriction
  - Fixed incompatibility with new `ggplot2` 3.3.0
  - Fixed allowing `NA`s in RxODE dataset
  - Fixed setting all compartment default values for bioavailability, rate, etc.
  - Added additional protection against floating point -> NaN for power functions

# RxODE v0.9.1-8
* Added the ability to have an input parameter to be assigned to a new
  value (Issue #135)
* Added LINPACK authors as contributors
* Added a `NEWS.md` file to track changes to the package
