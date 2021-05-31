Update Version: 1.0.1

Update addresses 2 notes and 1 error discovered by CRAN package check run 2021-05-30 23:49:52 CEST. Solutions to points brought up are listed below:

Check: installed package size
Result: NOTE
     installed size is 6.3Mb
     sub-directories of 1Mb or more:
     libs 5.6Mb
     
* We further cleaned the package to remove any garbage. Unfortunately, the package is still this note but we believe that any further removals from the package will hamper usability.

Check: dependencies in R code
Result: NOTE
    Namespace in Imports field not imported from: ‘RcppProgress’
     All declared Imports should be used.
     
* We moved ‘RcppProgress’ from "Imports" to "Suggests" within the DESCRIPTION.

Check: whether package can be installed
Result: ERROR
    Installation failed.
Flavor: r-patched-solaris-x86

* This error is produced due to problems between Solaris and testthat (r-lib/testthat#1257). To prevent this error we have exluded the unit tests via compiler macros.

---

# spectre v1.0.0

## Resubmission
This is a resubmission. Below is a list of changes made to address each comment:

Found the following (possibly) invalid URLs: URL: https://www.tidyverse.org/lifecycle/#stable (moved to https://lifecycle.r-lib.org/articles/stages.html)
From: README.md
Status: 200
Message: OK

Please change http --> https, add trailing slashes, or follow moved content as appropriate.

* The link was updated

Is there some reference about the method you can add in the Description field in the form Authors (year) <doi:.....>?

* We added a reference about the method in the Description field

## Test environments
* Local
  * macOS BigSur, R 4.0.5
  * Windows 10 Home, R 4.0.3
  * Ubuntu 20.04 LTS, R 4.0.5
* GitHub Actions 
  * windows-latest, R: 'release'
  * macOS-latest, R: 'release'
  * ubuntu-20.04, R: 'release'
* win-builder (release, devel)

## R CMD check results

0 errors | 0 warnings | 0 note

* This is a new release.

## Reverse dependencies
There are currently no reverse dependencies.